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
    values <- sprintf(paste0("mock_%0", width, "d"), seq_len(n))
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
