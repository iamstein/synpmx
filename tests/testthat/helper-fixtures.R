pmx_fixture <- function(n = 7L, id_type = c("integer", "character"),
                        duplicated_profiles = FALSE, missing_dv = FALSE) {
  id_type <- match.arg(id_type)
  pieces <- lapply(seq_len(n), function(subject) {
    time <- c(0, 0, 1, 2, 12, 12, 13, 14)
    evid <- c(1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L)
    base <- if (duplicated_profiles) 0 else subject * 0.15
    dv <- c(0, 0.2 + base, 4 + base, 2 + base,
            0, 0.3 + base, 5 + base, 2.5 + base)
    mdv <- ifelse(evid == 0L, 0L, 1L)
    if (missing_dv && subject == 1L) {
      dv[4L] <- NA_real_
      mdv[4L] <- 1L
    }
    data.frame(
      ID = if (id_type == "integer") as.integer(subject) else
        paste0("S", subject),
      TIME = time,
      DV = dv,
      AMT = ifelse(evid == 1L, 100 + subject, 0),
      RATE = 0,
      EVID = evid,
      CMT = ifelse(evid == 1L, 1L, 2L),
      MDV = mdv,
      WT = if (duplicated_profiles) 70 else 55 + subject * 4,
      AGE = if (duplicated_profiles) 40L else as.integer(25 + subject),
      SEX = if (duplicated_profiles) "female" else
        ifelse(subject %% 2L, "female", "male"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  out$SEX <- factor(out$SEX, levels = c("female", "male"))
  rownames(out) <- NULL
  out
}

fixture_roles <- function() {
  pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", mdv = "MDV", rate = "RATE",
    covariates = c("WT", "AGE", "SEX")
  )
}

multidvid_fixture <- function(n = 7L) {
  pieces <- lapply(seq_len(n), function(subject) {
    data.frame(
      id = as.integer(subject),
      time = c(0, 0, 0, 1, 1, 2, 2, 4, 4),
      amt = c(100 + subject, rep(0, 8)),
      dv = c(0, 0.2 + subject / 10, 80 - subject,
             4 + subject / 10, 60 - subject,
             2 + subject / 10, 45 - subject,
             0.8 + subject / 10, 35 - subject),
      dvid = factor(c("cp", "cp", "pd", "cp", "pd", "cp", "pd", "cp", "pd"),
                    levels = c("cp", "pd")),
      evid = c(1L, rep(0L, 8)),
      wt = rep(55 + 3 * subject, 9),
      sex = factor(rep(ifelse(subject %% 2L, "female", "male"), 9),
                   levels = c("female", "male"))
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

multidvid_roles <- function() {
  pmx_roles(
    id = "id", time = "time", dv = "dv", amt = "amt",
    evid = "evid", dvid = "dvid", covariates = c("wt", "sex")
  )
}
