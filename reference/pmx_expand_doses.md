# Expand or compress repeated dose records

`pmx_expand_doses()` writes out the doses that `ADDL` and `II` imply,
one row each, so that every dose the patient received is a row.
`pmx_compress_doses()` does the reverse: a run of dose records that
agree in every column except time and sit at one interval becomes one
record carrying `ADDL` and `II`.

## Usage

``` r
pmx_expand_doses(data, roles)

pmx_compress_doses(data, roles)
```

## Arguments

- data:

  A data frame of pharmacometric events.

- roles:

  A `pmx_roles` object declaring `addl` and `ii`.

## Value

`data`, with the repeated doses written out or folded back, and the
`addl` and `ii` columns set accordingly. Rows are ordered by subject and
time, with a dose before an observation at the same time.

## Details

Every generator calls both itself – expansion on the source it reads,
compression on the data it returns – so the synthetic study comes back
in the encoding its source used. They are exported for the caller who
wants to see what the generator saw, or to hand `nlmixr2` the compact
form directly.

Neither does anything unless
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
declares `addl` and `ii`. A source that writes every dose out is left
exactly as it is.

Expansion places dose `k` at `time + k * ii`, and moves a declared
`nominal_time` by the same amount. Compression requires the recorded and
the nominal gaps to agree, so a schedule whose recorded times drift off
the plan is not compressed – `ADDL` can only place a dose on an exact
interval, and writing one would change the study.

## Examples

``` r
roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                   evid = "EVID", addl = "ADDL", ii = "II")
compact <- data.frame(
  ID = 1, TIME = c(0, 12, 36), DV = c(NA, 3.1, 2.4),
  AMT = c(100, 0, 0), EVID = c(1, 0, 0), ADDL = c(2, 0, 0), II = c(24, 0, 0)
)
expanded <- pmx_expand_doses(compact, roles)
nrow(expanded)                        # three doses, two observations
#> [1] 5
identical(pmx_compress_doses(expanded, roles)$ADDL, compact$ADDL)
#> [1] TRUE
```
