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
  if (!is.null(roles$mdv)) {
    out <- out & .is_zero(data[[roles$mdv]])
  }
  if (require_present) out <- out & !is.na(data[[roles$dv]])
  out
}

.endpoint <- function(data, roles) {
  if (is.null(roles$dvid)) {
    rep("DV", nrow(data))
  } else {
    out <- as.character(data[[roles$dvid]])
    out[is.na(out)] <- "<missing>"
    out
  }
}

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

.unique_in_order <- function(x) x[!duplicated(x)]

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

.new_ids <- function(source_ids, n) {
  template <- source_ids
  if (is.factor(template)) {
    width <- max(3L, nchar(as.character(n)))
    labels <- sprintf(paste0("mock_%0", width, "d"), seq_len(n))
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
    labels <- sprintf(paste0("mock_%0", width, "d"), seq_len(n))
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

.ar1_noise <- function(n, phi, sd) {
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
