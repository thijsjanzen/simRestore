# simRestore 1.1.6

All changes in this version were written by Claude (Anthropic AI assistant).

## Bug fixes

* `optimize_static()`: fixed a test failure on CRAN (r-release-linux-x86_64,
  gcc 16) where optimizing only pulling could return a truncated `results`
  tibble and `final_freq = NaN`. The best fit was re-simulated with fresh
  random numbers, so the population could go extinct in the final run even
  though it survived during optimization.
* `optimize_static()`: the extinction check after optimizing pulling checked
  `result$results` before it was assigned, so it never triggered. It now
  checks the returned simulation.
* `optimize_static()`: when no surviving pull level was found, the search
  interval was shrunk on the natural scale instead of on log10 scale, which
  could yield nonsensical pull values (e.g. 1e71).
* `optimize_static()`: the arguments `sd_starting_freq` and `sex_ratio_pull`
  were silently ignored; they are now passed on to `simulate_policy()`.

## New features

* `optimize_static()` gains a `seed` argument. The same seed is used for every
  evaluation during optimization (common random numbers), making the
  objective deterministic and the returned results reproducible.

## Tests

* The `optimize_static()` tests now set a seed, check that `final_freq` is
  finite, and include a new test for seed reproducibility.

## Other

* Removed the compiled `src/simRestore.dll` from version control.
