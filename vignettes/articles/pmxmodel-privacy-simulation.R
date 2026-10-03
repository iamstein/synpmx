# Membership inference against two simplified models of the release
# `synpmx_model()` makes. Read by pmxmodel-privacy-simulation.Rmd, which works
# through each step, and by pmxmodel-privacy.Rmd, which quotes the results.
#
# A member is a patient in the simulated study and a non-member a new patient
# from the same population. The attacker holds the candidate's own values and
# reference values for the population, and scores the candidate by how far the
# release has moved from the reference toward them. Each protection is scored
# on the same simulated studies with and without it, so two columns differ only
# by what the protection does.

# The probability that a member scores above a non-member, ties counted half:
# the area under the receiver operating characteristic curve.
auc <- function(members, others) {
  ranks <- rank(c(members, others))
  n <- length(members)
  (sum(ranks[seq_len(n)]) - n * (n + 1) / 2) / (n * length(others))
}

# Attendance, one row per patient and one column per scheduled visit: 1 where
# the patient was sampled there and 0 where the sample was missed. The chance of
# a sample differs by visit, and `correlation` is shared between one patient's
# visits: a patient who misses one visit is more likely to miss the next.
draw_attendance <- function(patients, rate, correlation) {
  latent <- sqrt(correlation) * rnorm(patients) +
    sqrt(1 - correlation) * matrix(rnorm(patients * length(rate)), patients)
  1 * sweep(latent, 2, stats::qnorm(rate), "<")
}

# The share of patients sampled at each visit, as the visit model releases it.
# With `floor = 3`, a share with one or two patients on either side is rounded
# to 0 or 1, as the release does; `floor = 1` releases every share as it is.
release_shares <- function(attendance, floor) {
  patients <- nrow(attendance)
  held <- colSums(attendance)
  released <- held / patients
  thin <- pmin(held, patients - held) > 0 &
    pmin(held, patients - held) < floor
  released[thin] <- as.numeric(released[thin] >= 0.5)
  released
}

# The attacker's score for one candidate: at each visit, whether the candidate
# was sampled minus the visit's average attendance, times the released share
# minus that average, summed over visits. It is high when the visits the
# candidate missed are the ones whose released attendance came out below
# average, and the ones they kept came out above.
attendance_score <- function(candidate, rate, released) {
  sum((candidate - rate) * (released - rate))
}

# One member's and one non-member's score from each of `draws` simulated
# studies: `patients` patients behind each share, scheduled at `visits` visits,
# with an average attendance between 0.2 and 0.8 at each, and the attacker's
# reference the true average. The release pools a share over the arms that
# have the visit, so the patients behind it are usually the whole cohort.
attendance_scores <- function(patients, visits, correlation, floor,
                              draws = 1500) {
  scores <- vapply(seq_len(draws), function(draw) {
    rate <- runif(visits, 0.2, 0.8)
    study <- draw_attendance(patients, rate, correlation)
    outsider <- draw_attendance(1, rate, correlation)[1, ]
    released <- release_shares(study, floor)
    c(attendance_score(study[1, ], rate, released),
      attendance_score(outsider, rate, released))
  }, numeric(2))
  data.frame(member = scores[1, ], non_member = scores[2, ])
}

# The same studies scored with every share released and with the floor.
attendance_attack <- function(patients, visits, correlation, draws = 1500) {
  scores <- vapply(seq_len(draws), function(draw) {
    rate <- runif(visits, 0.2, 0.8)
    study <- draw_attendance(patients, rate, correlation)
    outsider <- draw_attendance(1, rate, correlation)[1, ]
    whole <- release_shares(study, 1)
    floored <- release_shares(study, 3)
    c(attendance_score(study[1, ], rate, whole),
      attendance_score(outsider, rate, whole),
      attendance_score(study[1, ], rate, floored),
      attendance_score(outsider, rate, floored))
  }, numeric(4))
  c(without_floor = auc(scores[1, ], scores[2, ]),
    released = auc(scores[3, ], scores[4, ]))
}

# The population model: five parameters of a two-compartment oral model, each
# released as a typical value and a between-subject variance. On the log scale
# a typical value is close to the mean of the patients' random effects and the
# variance close to the mean of their squares. The attacker holds the
# candidate's own random effects, and `off_by` is the SD, on the log scale, of
# the error in the attacker's reference values. The same studies are scored
# with the release unrounded and at two significant figures.
population_attack <- function(patients, off_by, draws = 2000) {
  typical <- c(cl = 4.3, v = 31, q = 2.1, v2 = 55, ka = 1.2)
  omega <- c(cl = 0.3, v = 0.25, q = 0.4, v2 = 0.35, ka = 0.6)^2
  scores <- vapply(seq_len(draws), function(draw) {
    eta <- vapply(omega, function(w) rnorm(patients, 0, sqrt(w)),
                  numeric(patients))
    outsider <- vapply(omega, function(w) rnorm(1, 0, sqrt(w)), numeric(1))
    reference_typical <- typical * exp(rnorm(5, 0, off_by))
    reference_omega <- omega * exp(rnorm(5, 0, 2 * off_by))
    score <- function(e, released_typical, released_omega) {
      sum(e * log(released_typical / reference_typical) / reference_omega) +
        sum((e^2 / reference_omega - 1) *
              (released_omega / reference_omega - 1)) / 2
    }
    exact_typical <- typical * exp(colMeans(eta))
    exact_omega <- colMeans(eta^2)
    rounded_typical <- signif(exact_typical, 2)
    rounded_omega <- signif(exact_omega, 2)
    c(score(eta[1, ], exact_typical, exact_omega),
      score(outsider, exact_typical, exact_omega),
      score(eta[1, ], rounded_typical, rounded_omega),
      score(outsider, rounded_typical, rounded_omega))
  }, numeric(4))
  c(unrounded = auc(scores[1, ], scores[2, ]),
    released = auc(scores[3, ], scores[4, ]))
}

# Both tables, from one seed, so that every article quoting them quotes the
# same numbers.
membership_tables <- function() {
  set.seed(20261002)
  population <- expand.grid(patients = c(12, 32, 60), off_by = c(0, 0.1, 0.2))
  population <- cbind(population, t(mapply(
    population_attack, population$patients, population$off_by)))
  attendance <- expand.grid(visits = c(5, 20, 60), patients = c(10, 30, 100),
                            correlation = c(0, 0.6))
  attendance <- cbind(attendance, t(mapply(
    attendance_attack, attendance$patients, attendance$visits,
    attendance$correlation)))
  list(population = population, attendance = attendance)
}
