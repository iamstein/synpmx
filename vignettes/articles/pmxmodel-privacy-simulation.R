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

# Yes/no outcomes, one row per patient and one column per visit. The rate of a
# yes differs by visit, and `correlation` is shared between one patient's
# visits: a patient with a yes at one visit is more likely to have one at the
# next.
draw_outcomes <- function(patients, rate, correlation) {
  latent <- sqrt(correlation) * rnorm(patients) +
    sqrt(1 - correlation) * matrix(rnorm(patients * length(rate)), patients)
  1 * sweep(latent, 2, stats::qnorm(rate), "<")
}

# The arm's share of patients with a yes at each visit. With `floor = 3`, a
# share with one or two patients on either side is rounded to 0 or 1, as the
# release does; `floor = 1` releases every share as it is.
release_shares <- function(outcomes, floor) {
  patients <- nrow(outcomes)
  held <- colSums(outcomes)
  released <- held / patients
  thin <- pmin(held, patients - held) > 0 &
    pmin(held, patients - held) < floor
  released[thin] <- as.numeric(released[thin] >= 0.5)
  released
}

# The attacker's score for one candidate: at each visit, the candidate's
# outcome minus the reference rate, times the released share minus the
# reference rate, summed over visits.
frequency_score <- function(candidate, rate, released) {
  sum((candidate - rate) * (released - rate))
}

# One member's and one non-member's score from each of `draws` simulated
# studies: `patients` patients behind each share, measured at `visits` visits,
# with a rate between 0.2 and 0.8 at each visit, and the attacker's reference
# the true rate. The release pools a share over the arms that have the visit,
# so the patients behind it are usually the whole cohort.
frequency_scores <- function(patients, visits, correlation, floor,
                             draws = 1500) {
  scores <- vapply(seq_len(draws), function(draw) {
    rate <- runif(visits, 0.2, 0.8)
    arm <- draw_outcomes(patients, rate, correlation)
    outsider <- draw_outcomes(1, rate, correlation)[1, ]
    released <- release_shares(arm, floor)
    c(frequency_score(arm[1, ], rate, released),
      frequency_score(outsider, rate, released))
  }, numeric(2))
  data.frame(member = scores[1, ], non_member = scores[2, ])
}

# The same studies scored with every share released and with the floor.
frequency_attack <- function(patients, visits, correlation, draws = 1500) {
  scores <- vapply(seq_len(draws), function(draw) {
    rate <- runif(visits, 0.2, 0.8)
    arm <- draw_outcomes(patients, rate, correlation)
    outsider <- draw_outcomes(1, rate, correlation)[1, ]
    whole <- release_shares(arm, 1)
    floored <- release_shares(arm, 3)
    c(frequency_score(arm[1, ], rate, whole),
      frequency_score(outsider, rate, whole),
      frequency_score(arm[1, ], rate, floored),
      frequency_score(outsider, rate, floored))
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
  frequency <- expand.grid(visits = c(5, 20, 60), patients = c(10, 30, 100),
                           correlation = c(0, 0.6))
  frequency <- cbind(frequency, t(mapply(
    frequency_attack, frequency$patients, frequency$visits,
    frequency$correlation)))
  list(population = population, frequency = frequency)
}
