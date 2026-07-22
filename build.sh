#!/usr/bin/env bash
#
# Full verification for pmxSynthData.
#
#   ./build.sh              regenerate docs, build tarball, R CMD check
#   ./build.sh check        same as above (explicit)
#   ./build.sh vignettes    also install and render every vignette to HTML
#   ./build.sh --keep-lib   leave the temporary library in place for inspection
#
# Everything happens against a temporary library that is created fresh on every
# run and placed ahead of the user library, so nothing is ever validated
# against a previously installed pmxSynthData, an already loaded namespace, or
# stale rendered HTML. Suggested dependencies (ggplot2, opendp, nlmixr2data,
# ...) still resolve from the user library; only pmxSynthData itself is forced
# to be freshly built.
#
# R CMD check rebuilds the vignettes once, which is enough to prove they still
# execute. The separate 'vignettes' mode exists for when you want inspectable
# HTML to read the tables and plots.
#
# Artifacts land in output/, which is gitignored and excluded from the build.

set -euo pipefail

readonly PKG="pmxSynthData"
readonly ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly OUT="${ROOT}/output"
readonly CHECK_DIR="${OUT}/check"
readonly VIG_DIR="${OUT}/vignettes"
readonly LOG_DIR="${OUT}/logs"

MODE="all"
KEEP_LIB=0
for arg in "$@"; do
  case "$arg" in
    check)       MODE="check" ;;
    vignettes)   MODE="vignettes" ;;
    --keep-lib)  KEEP_LIB=1 ;;
    -h|--help)   sed -n '2,/^$/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)           echo "unknown argument: $arg (try --help)" >&2; exit 2 ;;
  esac
done

RED=$'\033[31m'; GREEN=$'\033[32m'; YELLOW=$'\033[33m'; BOLD=$'\033[1m'; OFF=$'\033[0m'
step() { printf '\n%s==> %s%s\n' "$BOLD" "$1" "$OFF"; }
fail() { printf '%sFAIL%s %s\n' "$RED" "$OFF" "$1" >&2; exit 1; }
ok()   { printf '%sok%s   %s\n' "$GREEN" "$OFF" "$1"; }
warn() { printf '%swarn%s %s\n' "$YELLOW" "$OFF" "$1"; }

LIB="$(mktemp -d "${TMPDIR:-/tmp}/${PKG}-lib.XXXXXX")"
cleanup() {
  if [[ $KEEP_LIB -eq 1 ]]; then
    printf '\ntemporary library kept at %s\n' "$LIB"
  else
    rm -rf "$LIB"
  fi
}
trap cleanup EXIT

mkdir -p "$CHECK_DIR" "$VIG_DIR" "$LOG_DIR"
cd "$ROOT"

# R CMD check reruns tests and rebuilds vignettes; keep it from also picking up
# a half-written user-level Renviron or a stale .RData.
export R_LIBS_USER="${LIB}:$(Rscript -e 'cat(paste(.libPaths(), collapse=":"))' 2>/dev/null || true)"
export _R_CHECK_CRAN_INCOMING_=false
export _R_CHECK_FORCE_SUGGESTS_=false

step "Regenerating roxygen documentation"
Rscript -e 'suppressMessages(roxygen2::roxygenise("."))' 2>&1 | tee "${LOG_DIR}/roxygen.log"
if ! git diff --quiet -- NAMESPACE man 2>/dev/null; then
  warn "NAMESPACE/man changed — regenerated docs were out of date, review before committing"
fi

step "Building source tarball"
rm -f "${OUT}"/${PKG}_*.tar.gz
R CMD build --no-manual . 2>&1 | tee "${LOG_DIR}/build.log"
TARBALL="$(ls -t ${PKG}_*.tar.gz | head -1)"
[[ -f "$TARBALL" ]] || fail "no tarball produced"
mv "$TARBALL" "${OUT}/"
TARBALL="${OUT}/${TARBALL}"
ok "$(basename "$TARBALL")"

if [[ "$MODE" != "vignettes" ]]; then
  step "R CMD check (as-CRAN, into a clean library)"
  set +e
  R CMD check --library="$LIB" --output="$CHECK_DIR" --as-cran --no-manual "$TARBALL" \
    2>&1 | tee "${LOG_DIR}/check.log"
  CHECK_STATUS=${PIPESTATUS[0]}
  set -e
  LOG="${CHECK_DIR}/${PKG}.Rcheck/00check.log"
  if [[ -f "$LOG" ]]; then
    printf '\n%sR CMD check: %s%s\n' "$BOLD" "$(grep '^Status:' "$LOG" || echo 'Status: unknown')" "$OFF"
    # Surface the individual checks that were not OK, with their detail lines.
    awk '/^\* checking .*\.\.\. (NOTE|WARNING|ERROR)/ {show=1}
         /^\* checking .*\.\.\. OK$/ {show=0}
         /^\* DONE/ {show=0}
         show' "$LOG" || true
  fi
  [[ ${CHECK_STATUS} -eq 0 ]] || fail "R CMD check exited ${CHECK_STATUS} — see ${LOG}"
  ok "R CMD check passed"
fi

if [[ "$MODE" != "vignettes" ]]; then
  cat <<EOF

Vignettes were rebuilt once inside R CMD check. AGENTS.md does not require a
separate clean-library re-render, but it does require that vignette code and
prose be updated to match behavioral changes. Run './build.sh vignettes' if you
want inspectable HTML in ${VIG_DIR}.
EOF
  exit 0
fi

step "Installing the freshly built tarball into the clean library"
R CMD INSTALL --library="$LIB" "$TARBALL" 2>&1 | tee "${LOG_DIR}/install.log" >/dev/null
ok "installed to $LIB"

step "Rendering every vignette against the clean install"
INSTALLED=$(Rscript --vanilla -e "cat(as.character(utils::packageVersion('${PKG}', lib.loc='${LIB}')))")
printf 'using %s %s from %s\n' "$PKG" "$INSTALLED" "$LIB"

RENDER_FAILED=0
for RMD in vignettes/*.Rmd; do
  NAME="$(basename "${RMD%.Rmd}")"
  printf '  %-42s ' "$NAME"
  if Rscript --vanilla -e "
      .libPaths(c('${LIB}', .libPaths()))
      # Guard against silently validating a stale install from the user library.
      stopifnot(identical(
        normalizePath(dirname(getNamespaceInfo(asNamespace('${PKG}'), 'path'))),
        normalizePath('${LIB}')))
      rmarkdown::render('${RMD}',
        output_dir = '${VIG_DIR}',
        quiet = TRUE,
        envir = new.env(parent = globalenv()))
    " >"${LOG_DIR}/${NAME}.log" 2>&1; then
    printf '%sok%s\n' "$GREEN" "$OFF"
  else
    printf '%sFAILED%s (see %s)\n' "$RED" "$OFF" "${LOG_DIR}/${NAME}.log"
    RENDER_FAILED=1
  fi
done

[[ $RENDER_FAILED -eq 0 ]] || fail "at least one vignette failed to render"

step "Summary"
ok "tarball    ${TARBALL}"
[[ "$MODE" != "vignettes" ]] && ok "check      ${CHECK_DIR}/${PKG}.Rcheck/00check.log"
ok "vignettes  ${VIG_DIR}"
cat <<'EOF'

A successful knit is necessary but not sufficient. Open the rendered HTML and
visually inspect every table and plot before treating this as verified.
EOF
