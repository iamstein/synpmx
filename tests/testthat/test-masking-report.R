test_that("pmx_masking_report(section =) keeps whole blocks and nothing else", {
  data <- pmx_simulated_fixture(30)
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "DVID", covariates = "WT"
  )
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))

  full <- pmx_masking_report(synthetic, data, roles)
  dose <- pmx_masking_report(synthetic, data, roles, section = "dose_schedules")

  # A section is a strict subset of the full report, header row included, and
  # the rows keep the order and the wording they have in it.
  expect_lt(nrow(dose), nrow(full))
  expect_true(all(dose$Quantity %in% full$Quantity))
  expect_equal(dose$Quantity[[1L]],
               "**Dose schedules: WHEN each patient was dosed**")
  expect_identical(dose, full[match(dose$Quantity, full$Quantity), ],
                   ignore_attr = "row.names")

  # Asking for every section is asking for the whole report.
  every <- pmx_masking_report(synthetic, data, roles,
                              section = names(synpmx:::.masking_sections))
  expect_equal(nrow(every), nrow(full))

  # Two sections come back in report order, not in the order they were asked
  # for, so a caller cannot silently reorder the reading.
  pair <- pmx_masking_report(synthetic, data, roles,
                             section = c("dose_schedules", "anchors"))
  expect_equal(pair$Quantity[[1L]], "**Who was available to build on**")

  expect_error(pmx_masking_report(synthetic, data, roles, section = "doses"),
               "Unknown `section`")
})
