# A scorecard as a coloured HTML table

Displays a
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
as an interactive table with the verdicts coloured: `"FAIL"` in bold red
on a light red, `"review"` in bold orange on a light orange,
`"unavailable"` in muted grey, and `"pass"` left as ordinary text. A
card is five verdicts among thirty-odd rows of prose, and the rows that
need reading are the ones that have to be findable without reading all
of it.

## Usage

``` r
synpmx_scorecard_datatable(x, ...)
```

## Arguments

- x:

  A
  [`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).

- ...:

  Passed to
  [`DT::datatable()`](https://rdrr.io/pkg/DT/man/datatable.html). Paging
  is off and row numbers are suppressed by default, since the whole card
  is meant to be read at once and `check` already names each row.

## Value

An
[`htmltools::tagList`](https://rstudio.github.io/htmltools/reference/tagList.html)
holding the coloured card and the notes that knitting one carries.
Without `DT` installed, `x` invisibly, having printed it.

## Details

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
computes the card and this displays it, so the object is unchanged and
can be subset, saved or printed as usual. The colouring is the only
thing added.

What is emitted is what knitting the card itself emits, with the
colouring added: the card, then the B5 rare-level detail where a study
has any.

`DT` is a suggested package rather than a required one – `synpmx` has no
hard dependencies – so without it installed this says so and prints the
card in the console form instead. The verdicts are all still there; only
the colour is missing.

The same restriction applies as to the card itself. Rows reading
`"source"` or `"both"` were computed from real patient data, so an HTML
file holding the whole table belongs in the environment the source lives
in.

## See also

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).

## Examples

``` r
data <- pmx_simulated_fixture(30)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT"
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 9 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE, SEX.
#>   Declare a column in `keep` to carry it through verbatim.
synpmx_scorecard_datatable(synpmx_scorecard(data, synthetic, roles))
#> <div class="datatables html-widget html-fill-item" id="htmlwidget-ac96cb3ee4656e2e9ec3" style="width:100%;height:auto;"></div>
#> <script type="application/json" data-for="htmlwidget-ac96cb3ee4656e2e9ec3">{"x":{"filter":"none","vertical":false,"data":[["A1","A2","A3","A4","A5a","A5b","A6","B1a","B1b","B2","B3","B4a","B4b","B5","C1","C2","C3","D1"],["Synthetic table is a legal PMX dataset","Source is legal under the declared roles","Every endpoint survived","Cohort size survived","Observations per patient","Doses per patient","Discrete endpoints keeping their source scale","Avatars with a visit set nobody else shares","Avatars with a dose schedule nobody else shares","Synthetic patients unusual within their stratum","Adversarial accuracy inside its null interval","Generated time vectors copying an exposed real one","Generated DV vectors copying an exposed real one","Rare source levels copied into the output","Strata keeping their source size","Distinct dose-time schedules represented","Arms keeping their source endpoints","Values landing in the same range"],["synthetic","source","both","both","both","both","both","run settings","run settings","synthetic","both","both","both","both","both","run settings","both","both"],["TRUE","TRUE","2 of 2","30 -&gt; 30","14 -&gt; 14","2 -&gt; 2","no discrete endpoint","0","0","0 of 30","0.767 in [0.248, 0.692]","0","0","no categorical covariate or stratum","no strata declared","1 of 1","no strata declared","sd x1.4 on pd (furthest of 3)"],["pass","pass","pass","pass","pass","pass","pass","pass","pass","pass","review","pass","pass","pass","pass","pass","pass","review"],["validate_pmx(synthetic, roles)","validate_pmx(source, roles, strict = FALSE)","compare_pmx_distributions(source, synthetic, roles)","pmx_masking_report(synthetic, source, roles, section = \"anchors\")","compare_pmx_distributions(source, synthetic, roles)","pmx_masking_report(synthetic, source, roles, section = \"dose_schedules\")","pmx_endpoint_types(source, roles)","unmaskable_strata(source, roles)","unmaskable_strata(source, roles)","flag_identifiable_subjects(synthetic, roles)","compare_pmx_proximity(source, synthetic, roles)","skeleton_uniqueness(source, roles, coarsen_time = TRUE)","compare_pmx_proximity(source, synthetic, roles)","pmx_roles(strata = , covariates = )","pmx_roles(strata = )","pmx_masking_report(synthetic, source, roles, section = \"dose_schedules\")","pmx_roles(strata = )","compare_pmx_distributions(source, synthetic, roles, output = \"tables\")"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>check<\/th>\n      <th>question<\/th>\n      <th>reads<\/th>\n      <th>result<\/th>\n      <th>verdict<\/th>\n      <th>explore<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"paging":false,"columnDefs":[{"name":"check","targets":0},{"name":"question","targets":1},{"name":"reads","targets":2},{"name":"result","targets":3},{"name":"verdict","targets":4},{"name":"explore","targets":5}],"order":[],"autoWidth":false,"orderClasses":false,"rowCallback":"function(row, data, displayNum, displayIndex, dataIndex) {\nvar value=data[4]; $(this.api().cell(row, 4).node()).css({'font-weight':value == \"FAIL\" ? \"bold\" : value == \"review\" ? \"bold\" : \"normal\",'color':value == \"FAIL\" ? \"#B00020\" : value == \"review\" ? \"#B45309\" : value == \"unavailable\" ? \"#6C757D\" : \"inherit\",'background-color':value == \"FAIL\" ? \"#FDECEA\" : value == \"review\" ? \"#FFF4E5\" : \"transparent\"});\n}"}},"evals":["options.rowCallback"],"jsHooks":[]}</script>
```
