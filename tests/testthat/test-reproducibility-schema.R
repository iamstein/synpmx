test_that("the complete result is reproducible without changing global RNG", {
  source <- pmx_fixture()
  roles <- fixture_roles()

  set.seed(987)
  state <- .Random.seed
  first <- mock_pmx(source, roles, n_subjects = 4, seed = 123)
  expect_identical(.Random.seed, state)
  second <- mock_pmx(source, roles, n_subjects = 4, seed = 123)
  different <- mock_pmx(source, roles, n_subjects = 4, seed = 124)

  expect_identical(first, second)
  expect_false(identical(first, different))
})

test_that("requested subjects and genuinely new typed IDs are produced", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 11, seed = 8)
  ids <- unique(mock$ID)
  expect_length(ids, 11L)
  expect_type(ids, "integer")
  expect_false(any(ids %in% source$ID))

  character_source <- pmx_fixture(id_type = "character")
  character_mock <- mock_pmx(character_source, roles, n_subjects = 4, seed = 8)
  expect_type(character_mock$ID, "character")
  expect_false(any(unique(character_mock$ID) %in% character_source$ID))
})

test_that("schema, classes, factors, and constant covariates are retained", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 5, seed = 22)

  expect_identical(names(mock), names(source))
  expect_identical(vapply(mock, class, character(1)),
                   vapply(source, class, character(1)))
  expect_identical(levels(mock$SEX), levels(source$SEX))
  for (column in roles$covariates) {
    by_subject <- split(mock[[column]], mock$ID)
    expect_true(all(vapply(by_subject, function(x) length(unique(x)) == 1L,
                           logical(1))))
  }
  expect_true(is.list(attr(mock, "pmx_settings")))
  expect_false(inherits(mock, "pmxmock"))
})
