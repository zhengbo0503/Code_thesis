# Chapter 4 — Accurate Eigenvalues and Eigenvectors

Code for `chaps/eigenvalue_eigenvector.tex` and the two eigenvalue figures of
Appendix A, `chaps/complex_experiments.tex`. Originally from
[`Code_htwz25`](https://github.com/zhengbo0503/Code_htwz25) (eigenvalues, MATLAB and
Julia, [arXiv:2501.03742](https://arxiv.org/abs/2501.03742)) and
[`Code_twz26b`](https://github.com/zhengbo0503/Code_twz26b) (eigenvectors).
`complex/` is new work with no upstream repo.

## Experiments

| Thesis           | Label                                      | Script                                       | Writes                                                                    |
| ---------------- | ------------------------------------------ | -------------------------------------------- | ------------------------------------------------------------------------- |
| Figure 4.1       | `fig.ev-maxferrk-varykappa`                | `experiments/evals_varying_cond.m`           | `output/varykappa_mode1–5.csv`                                            |
| Figure 4.2       | `fig.ev-maxferrk-varydim`                  | `experiments/evals_varying_dim.m`            | `output/varydim_mode1–5.csv`                                              |
| Figure 4.3       | `fig.ev-special-mtx`                       | `experiments/evals_special_matrix.m`         | `output/random_matrix.csv`, `whiskycorr.csv`, `covariance.csv`, `lau.csv` |
| Figures 4.4, 4.5 | `fig.ev-timing-union`, `fig.ev-timing-all` | `experiments/evals_timing.m`                 | `output/timing_union.csv`, `timings_2.csv`                                |
| Figure 4.6       | `fig.evec-vary-kappa`                      | `experiments/evecs_varying_cond.m`           | `data/evecs_varying_cond.csv`                                             |
| Figure 4.7       | `fig.evec-vary-n`                          | `experiments/evecs_varying_dim.m`            | `data/evecs_varying_dim.csv`                                              |
| Figure 4.8       | `fig.special-matrix-evec`                  | `experiments/evecs_special_matrix.m`         | `data/evecs_special_matrix_{kms,lehmer}.csv`                              |
| Figure A.1       | `fig.app-complex-ev-kappa`                 | `complex/test_complex_evals_varying_kappa.m` | `complex/data/complex_evals_varying_kappa_full.{mat,csv}`                 |
| Figure A.2       | `fig.app-complex-ev-size`                  | `complex/test_complex_evals_varying_n.m`     | `complex/data/complex_evals_varying_n_full.{mat,csv}`                     |

## Running the tests

```bash
tools/run_smoke.sh                    # smoke_test.m + julia/smoke_test.jl
```

```matlab
cd ch4_eigenvalue_eigenvector/experiments   % scripts assume this working directory
evals_varying_cond                          % etc.
```

- Advanpix is required and not on the path by default:
  `addpath('/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac')`.
- Score headless runs by the printed `N/M passed` line, not the exit code.
- `evecs_varying_cond.m` and `evecs_varying_dim.m` call `pause()`. Add
  `tools/batch_pause_guard` to the path for unattended runs.
- `evals_timing.m` only plots. Run `cd julia && julia test.jl` first to write
  `julia/result/timing_{order,cnd}.csv`. Install the vendored package with
  `Pkg.develop`, not `Pkg.add`.
- `evals_special_matrix.m` reads `data/superconductivty+data/train.csv`.
- Complex checks are not in `run_smoke.sh`: run `validate_complex_evd` by hand.

## Renames from the original repos

| Was                                               | Now                                  |
| ------------------------------------------------- | ------------------------------------ |
| `Code_htwz25/MATLAB/test/test_varying_cnd.m`      | `experiments/evals_varying_cond.m`   |
| `Code_htwz25/MATLAB/test/test_varying_dim.m`      | `experiments/evals_varying_dim.m`    |
| `Code_htwz25/MATLAB/test/test_special_matrix.m`   | `experiments/evals_special_matrix.m` |
| `Code_htwz25/MATLAB/test/test_timing.m`           | `experiments/evals_timing.m`         |
| `Code_htwz25/MATLAB/src/`, `Code_twz26b/evdalgs/` | `shared/`                            |
| `Code_htwz25/Julia/`                              | `julia/`                             |
| `Code_twz26b/test3.m`                             | `experiments/evecs_varying_cond.m`   |
| `Code_twz26b/test4.m`                             | `experiments/evecs_varying_dim.m`    |
| `Code_twz26b/test6.m`                             | `experiments/evecs_special_matrix.m` |

`compute_error` is split into `compute_error_evals` and `compute_error_evecs`, which
were different functions sharing a name. `scond` takes the scaling type
(`'D'`/`'C'`) second, not the norm type.
