test_that("simulation evaluation registry covers every demonstration dataset", {
  expect_setequal(
    names(sim_eval_registry(include_optional = TRUE)),
    c(
      "censoring",
      if (requireNamespace("nlmixr2data", quietly = TRUE)) {
        c("theo_md", "warfarin", "wbcSim", "nimoData", "mavoglurant")
      }
    )
  )
})

test_that("censoring simulation evaluation passes every hard gate", {
  case <- sim_eval_case("censoring")
  model <- sim_eval_fit(case)
  synthetic <- generate_pmx(model, seed = 401)
  gates <- sim_eval_gate_results(case, model, synthetic)
  expect_true(all(gates$pass), info = paste(
    gates$gate[!gates$pass], gates$detail[!gates$pass], collapse = "; "
  ))
  metrics <- sim_eval_metric_rows(case, model, synthetic, 401, 5, "public")
  expect_true(all(metrics$gates_failed == 0L))
})

test_that("named-data simulations pass the shared evaluation gates", {
  skip_if_not_installed("nlmixr2data")
  for (id in c(
    "theo_md", "warfarin", "wbcSim", "nimoData", "mavoglurant"
  )) {
    case <- sim_eval_case(id)
    model <- sim_eval_fit(case)
    synthetic <- generate_pmx(model, seed = 401)
    gates <- sim_eval_gate_results(case, model, synthetic)
    expect_true(all(gates$pass), info = paste(
      id, paste(gates$gate[!gates$pass], collapse = ", ")
    ))
    expect_true(sim_eval_public_overrides_absent(model$public$design))
  }
})

test_that("evaluation plot data exclude events and preserve panel semantics", {
  skip_if_not_installed("ggplot2")
  case <- sim_eval_case("censoring")
  model <- sim_eval_fit(case)
  synthetic <- generate_pmx(model, seed = 402)
  plotted <- sim_eval_plot_data(
    synthetic, case$roles, case$endpoints, "Synthetic",
    time_bounds = case$bounds$time
  )
  expect_equal(
    nrow(plotted), sum(sim_eval_observation_rows(synthetic, case$roles))
  )
  expect_identical(levels(plotted$dataset), c("Source", "Synthetic"))
  expect_setequal(unique(plotted$endpoint), names(case$endpoints))
  plot <- sim_eval_plot(case, synthetic)
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$facet$params$free$y, TRUE)
  y_scale <- plot$scales$get_scales("y")
  expect_true(is.null(y_scale) || identical(y_scale$trans$name, "identity"))
})
