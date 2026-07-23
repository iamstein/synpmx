`%||%` <- function(x, y) if (is.null(x)) y else x

.is_zero <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) || is.integer(x)) {
    !is.na(x) & x == 0
  } else {
    !is.na(x) & trimws(as.character(x)) %in% c("0", "0.0")
  }
}

.observation_rows <- function(data, roles, require_present = FALSE) {
  out <- .is_zero(data[[roles$evid]])
  if (!is.null(roles$mdv)) out <- out & .is_zero(data[[roles$mdv]])
  if (require_present) out <- out & !is.na(data[[roles$dv]])
  out
}

.event_rows <- function(data, roles) {
  !is.na(data[[roles$evid]]) & !.is_zero(data[[roles$evid]])
}

.dose_rows <- function(data, roles) {
  out <- .event_rows(data, roles)
  if (!is.null(roles$amt)) {
    amount <- suppressWarnings(as.numeric(data[[roles$amt]]))
    positive <- out & is.finite(amount) & amount > 0
    if (any(positive)) out <- positive
  }
  out
}

.endpoint <- function(data, roles) {
  if (is.null(roles$dvid)) rep("DV", nrow(data)) else {
    out <- as.character(data[[roles$dvid]])
    out[is.na(out)] <- "<missing>"
    out
  }
}

.unique_in_order <- function(x) x[!duplicated(x)]

.with_local_seed <- function(seed, code) {
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)) {
    stop("`seed` must be one integer from 0 to `.Machine$integer.max`.",
         call. = FALSE)
  }
  seed <- as.integer(seed)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.ar1_noise <- function(n, phi = 0.65, sd = 1) {
  if (!n || sd == 0) return(numeric(n))
  out <- numeric(n)
  out[1L] <- stats::rnorm(1L, sd = sd)
  if (n > 1L) {
    innovation_sd <- sd * sqrt(max(0, 1 - phi^2))
    for (i in 2:n) {
      out[i] <- phi * out[i - 1L] + stats::rnorm(1L, sd = innovation_sd)
    }
  }
  out
}

.schema_names <- function(schema) {
  vapply(schema$columns, `[[`, character(1), "name")
}

.schema_column <- function(schema, name) {
  index <- match(name, .schema_names(schema))
  if (is.na(index)) NULL else schema$columns[[index]]
}

.typed_missing <- function(column, n) {
  class <- column$class
  if ("factor" %in% class) {
    return(factor(rep(NA_character_, n), levels = column$levels,
                  ordered = isTRUE(column$ordered)))
  }
  if ("integer" %in% class || identical(column$typeof, "integer")) {
    return(rep(NA_integer_, n))
  }
  if ("numeric" %in% class || identical(column$typeof, "double")) {
    return(rep(NA_real_, n))
  }
  if ("logical" %in% class) return(rep(NA, n))
  if ("character" %in% class) return(rep(NA_character_, n))
  rep(NA, n)
}

.cast_public_column <- function(x, column, is_id = FALSE) {
  class <- column$class
  if ("factor" %in% class) {
    values <- as.character(x)
    levels <- column$levels
    if (is_id) levels <- unique(c(levels, values[!is.na(values)]))
    return(factor(values, levels = levels, ordered = isTRUE(column$ordered)))
  }
  if ("integer" %in% class || identical(column$typeof, "integer")) {
    return(as.integer(round(as.numeric(x))))
  }
  if ("numeric" %in% class || identical(column$typeof, "double")) {
    return(as.double(x))
  }
  if ("logical" %in% class) return(as.logical(x))
  if ("character" %in% class) return(as.character(x))
  x
}

.restore_public_schema <- function(result, schema, roles) {
  names <- .schema_names(schema)
  missing <- setdiff(names, names(result))
  for (name in missing) {
    result[[name]] <- .typed_missing(.schema_column(schema, name), nrow(result))
  }
  result <- result[, names, drop = FALSE]
  for (name in names) {
    result[[name]] <- .cast_public_column(
      result[[name]], .schema_column(schema, name),
      is_id = identical(name, roles$id)
    )
  }
  data_class <- schema$data_class
  if (!identical(data_class, "data.frame")) class(result) <- data_class
  result
}

.new_public_ids <- function(schema, id_name, n) {
  column <- .schema_column(schema, id_name)
  if (is.null(column)) stop("The public schema omits the ID role.", call. = FALSE)
  if ("factor" %in% column$class || "character" %in% column$class) {
    width <- max(3L, nchar(as.character(n)))
    values <- sprintf(paste0("syn_%0", width, "d"), seq_len(n))
  } else if ("integer" %in% column$class || column$typeof == "integer") {
    values <- as.integer(100000000L + seq_len(n))
  } else if ("numeric" %in% column$class || column$typeof == "double") {
    values <- as.numeric(100000000 + seq_len(n))
  } else {
    stop("The public ID schema must be integer, numeric, character, or factor.",
         call. = FALSE)
  }
  .cast_public_column(values, column, is_id = TRUE)
}

.schema_matches <- function(data, schema, exclude = NULL) {
  expected <- .schema_names(schema)
  actual <- setdiff(names(data), exclude)
  if (!identical(actual, expected)) {
    stop("Source columns after exclusion must exactly match the declared public schema.",
         call. = FALSE)
  }
  for (name in expected) {
    declared <- .schema_column(schema, name)
    if (!identical(class(data[[name]]), declared$class)) {
      stop("Source class for `", name,
           "` does not match the declared public schema.", call. = FALSE)
    }
    if (is.factor(data[[name]]) &&
        !identical(levels(data[[name]]), declared$levels)) {
      stop("Source factor levels for `", name,
           "` do not match the declared public schema.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Every column a role points at, other than the ones `exclude` names for
# removal. This is the allowlist: what AVATAR retains, and what the
# identifier-name check treats as consciously handled.
.retained_role_columns <- function(roles) {
  referenced <- unlist(roles[setdiff(names(roles), "exclude")],
                       use.names = FALSE)
  unique(referenced[!is.na(referenced) & nzchar(referenced)])
}

.direct_identifier_names <- function(names) {
  lowered <- tolower(names)
  normalized <- gsub("[^a-z0-9]", "", lowered)
  exact <- c(
    "name", "firstname", "lastname", "fullname", "email", "emailaddress",
    "ssn", "socialsecuritynumber", "mrn", "medicalrecordnumber",
    "address", "streetaddress", "phone", "telephone", "dob", "birthdate",
    "dateofbirth", "passport", "patientname", "subjectname",
    "contactemail", "phonenumber", "telephonenumber"
  )
  tokens <- strsplit(gsub("[^a-z0-9]+", "_", lowered), "_", fixed = TRUE)
  token_patterns <- c(
    "name", "firstname", "lastname", "fullname", "email", "ssn", "mrn",
    "address", "phone", "telephone", "dob", "birthdate", "passport"
  )
  token_match <- vapply(tokens, function(x) any(x %in% token_patterns),
                        logical(1))
  names[normalized %in% exact | token_match]
}

.release_id <- function() {
  paste0("pmx-", format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC"), "-",
         Sys.getpid())
}

.recursive_names <- function(x) {
  out <- names(x) %||% character()
  if (is.list(x)) {
    for (value in x) out <- c(out, .recursive_names(value))
  }
  out
}

.jitter_times <- function(nominal, fraction) {
  if (!length(nominal) || fraction == 0) return(as.numeric(nominal))
  unique_grid <- sort(unique(nominal))
  spacing <- diff(unique_grid)
  spacing <- spacing[spacing > 0]
  sd <- if (length(spacing)) min(spacing) * fraction else fraction
  jitter <- stats::rnorm(length(unique_grid), sd = sd)
  as.numeric(nominal) + jitter[match(nominal, unique_grid)]
}

# AVATAR (Version 4) helpers ------------------------------------------------
#
# Restored from the Version 1 engine. Used by synpmx_avatar() in synthesis.R and
# profiles.R for template events, schema restoration, and diagnostics.
.aligned_time <- function(data, roles) {
  time <- as.numeric(data[[roles$time]])
  id <- data[[roles$id]]
  aligned <- time
  for (subject in .unique_in_order(id)) {
    rows <- !is.na(id) & id == subject
    event <- rows & !.is_zero(data[[roles$evid]])
    start <- event
    if (!is.null(roles$amt)) {
      amount <- data[[roles$amt]]
      positive <- event & !is.na(amount) & amount > 0
      if (any(positive)) start <- positive
    }
    if (any(start)) {
      origin <- min(time[start], na.rm = TRUE)
      aligned[rows] <- time[rows] - origin
    }
  }
  aligned
}

.first_present <- function(x) {
  present <- which(!is.na(x))
  if (!length(present)) return(x[NA_integer_][1L])
  x[present[1L]]
}

.representative_values <- function(data, roles, subjects = NULL) {
  if (is.null(subjects)) subjects <- .unique_in_order(data[[roles$id]])
  covariates <- roles$covariates
  result <- stats::setNames(vector("list", length(covariates)), covariates)
  for (covariate in covariates) {
    result[[covariate]] <- lapply(subjects, function(subject) {
      rows <- !is.na(data[[roles$id]]) & data[[roles$id]] == subject
      .first_present(data[[covariate]][rows])
    })
  }
  result
}

.new_ids <- function(source_ids, n) {
  template <- source_ids
  if (is.factor(template)) {
    width <- max(3L, nchar(as.character(n)))
    labels <- sprintf(paste0("syn_%0", width, "d"), seq_len(n))
    return(factor(labels, levels = c(levels(template), labels),
                  ordered = is.ordered(template)))
  }
  if (is.integer(template)) {
    start <- if (all(is.na(template))) 1L else max(template, na.rm = TRUE) + 1L
    return(as.integer(start + seq_len(n) - 1L))
  }
  if (is.numeric(template)) {
    start <- if (all(is.na(template))) 1 else max(template, na.rm = TRUE) + 1
    return(as.numeric(start + seq_len(n) - 1))
  }
  if (is.character(template)) {
    width <- max(3L, nchar(as.character(n)))
    labels <- sprintf(paste0("syn_%0", width, "d"), seq_len(n))
    while (any(labels %in% template)) labels <- paste0("new_", labels)
    return(labels)
  }
  stop("ID columns must be integer, numeric, character, or factor.",
       call. = FALSE)
}

.restore_column <- function(x, template, is_id = FALSE) {
  if (is.factor(template)) {
    values <- as.character(x)
    lev <- levels(template)
    if (is_id) lev <- unique(c(lev, values[!is.na(values)]))
    return(factor(values, levels = lev, ordered = is.ordered(template)))
  }
  if (inherits(template, "Date")) return(as.Date(x))
  if (inherits(template, "POSIXct")) {
    return(as.POSIXct(x, origin = "1970-01-01", tz = attr(template, "tzone")))
  }
  if (is.integer(template)) return(as.integer(round(x)))
  if (is.double(template)) return(as.double(x))
  if (is.logical(template)) return(as.logical(x))
  if (is.character(template)) return(as.character(x))
  x
}

.restore_schema <- function(result, source, roles) {
  result <- result[, names(source), drop = FALSE]
  for (column in names(source)) {
    result[[column]] <- .restore_column(
      result[[column]], source[[column]], is_id = identical(column, roles$id)
    )
  }
  source_class <- class(source)
  if (!identical(source_class, "data.frame")) class(result) <- source_class
  result
}

.warning_collector <- function() {
  env <- new.env(parent = emptyenv())
  env$messages <- character()
  env$add <- function(message) env$messages <- unique(c(env$messages, message))
  env
}

.weighted_available <- function(values, weights) {
  okay <- is.finite(values) & is.finite(weights) & weights >= 0
  if (!any(okay)) return(NA_real_)
  available_weights <- weights[okay]
  total <- sum(available_weights)
  if (!is.finite(total) || total <= 0) {
    available_weights <- rep(1 / sum(okay), sum(okay))
  } else {
    available_weights <- available_weights / total
  }
  sum(values[okay] * available_weights)
}

