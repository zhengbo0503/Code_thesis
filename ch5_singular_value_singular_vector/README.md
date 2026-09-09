# Chapter 5 — Accurate Singular Values and Singular Vectors

Code for `chaps/singular_value_singular_vector.tex` and the three singular-value
figures of Appendix A, `chaps/complex_experiments.tex`. Originally from
[`Code_twz26a`](https://github.com/zhengbo0503/Code_twz26a) (singular values,
[arXiv:2602.18134](https://arxiv.org/abs/2602.18134)) and
[`Code_twz26b`](https://github.com/zhengbo0503/Code_twz26b) (singular vectors).
`complex/` is new work with no upstream repo.

## Experiments

| Thesis     | Label                           | Script                                  | Writes                                                          |
| ---------- | ------------------------------- | --------------------------------------- | --------------------------------------------------------------- |
| Figure 5.1 | `fig.fwderr_diff_cond`          | `experiments/svals_varying_cond.m`      | `data/svals_varying_cond.{mat,csv}`                             |
| Figure 5.2 | `fig.fwderr_diff_size`          | `experiments/svals_varying_cols.m`      | `data/svals_varying_cols.{mat,csv}`                             |
| Figure 5.3 | `fig.special-matrix`            | `experiments/svals_special_matrix.m`    | `data/svals_special_matrix_{blur,kahan,whiskycorr,lauchli}.csv` |
| Figure 5.4 | `fig.compare-timing`            | `experiments/svals_timing.m`            | `data/svals_timing.{mat,csv}`                                   |
| Figure 5.5 | `fig.compare-timing-ssd`        | `experiments/svals_ssd.m`               | `data/svals_ssd.{mat,csv}`                                      |
| Figure 5.6 | `fig.svec-vary-kappa`           | `experiments/svecs_varying_cond.m`      | `data/svecs_varying_cond.{mat,csv}`                             |
| Figure 5.7 | `fig.svec-vary-n`               | `experiments/svecs_varying_cols.m`      | `data/svecs_varying_cols.{mat,csv}`                             |
| Figure 5.8 | `fig.special-matrix-svec`       | `experiments/svecs_special_matrix.m`    | `data/svecs_special_matrix_{kms,lehmer}.csv`                    |
| Figure 5.9 | `fig.svec-vary-kappa-left-svec` | `experiments/svecs_left_varying_cond.m` | `data/svecs_left_varying_cond.{mat,csv}`                        |
| Figure A.3 | `fig.app-complex-kappa`         | `complex/test_complex_varying_kappa.m`  | `complex/data/complex_varying_kappa_full.mat`                   |
| Figure A.4 | `fig.app-complex-size`          | `complex/test_complex_varying_n.m`      | `complex/data/complex_varying_n_full.mat`                       |
| Figure A.5 | `fig.app-complex-ssd`           | `complex/test_complex_ssd_timing.m`     | `complex/data/complex_ssd_timing_full.mat`                      |

## Running the tests

```bash
tools/run_smoke.sh      # smoke_test_mex.m, then smoke_test.m, then the complex preflight
```

```matlab
cd ch5_singular_value_singular_vector/experiments   % scripts assume this working directory
svals_varying_cond                                  % etc.
```

```matlab
cd ch5_singular_value_singular_vector/complex       % full complex run, unattended
setup_complex_paths
run_all_complex_full
```

- Run `smoke_test_mex` first; everything else depends on the MEX wrappers loading.
- Advanpix is required and not on the path by default:
  `addpath('/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac')`.
- Score headless runs by the printed `N/M passed` line, not the exit code.
- The shipped `shared/*.mexmaca64` work under my setup as-is and I believe it is not possible to run under other setups. To rebuild, use
  `mex/get_*/build_*_mex.m` with OpenBLAS ≤ 0.3.30 and check with `mex/test_xgesvj.m`.
- `svecs_varying_cond.m`, `svecs_varying_cols.m` and `svecs_left_varying_cond.m` call
  `pause()`. Add `tools/batch_pause_guard` to the path for unattended runs.
- `svals_special_matrix.m` needs Anymatrix v1.4.

## Renames from the original repos

| Was                                                    | Now                                                                                                                   |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------- |
| `Code_twz26a/test1.m` – `test5.m`                      | `experiments/svals_varying_cond.m`, `svals_varying_cols.m`, `svals_special_matrix.m`, `svals_timing.m`, `svals_ssd.m` |
| `Code_twz26b/test1.m`, `test2.m`, `test5.m`, `test7.m` | `experiments/svecs_varying_cond.m`, `svecs_varying_cols.m`, `svecs_special_matrix.m`, `svecs_left_varying_cond.m`     |
| `Code_twz26a/get_*/`                                   | `mex/get_*/`                                                                                                          |
| `Code_twz26a/`, `Code_twz26b/svdalgs/`                 | `shared/`                                                                                                             |

Names that collided across the two repos, now disambiguated:

- `mposj` and `mposj_svec` are **two different algorithms**, not two versions of one.
  `mposj_svec` is the QR + transposed variant used by the `svecs_*` scripts.
- `compute_error` is split into `compute_error_svals` and `compute_error_svecs`.
- `get_testmatrix` is split into `get_testmatrix_anymatrix` (blur, kahan, whiskycorr,
  lauchli) and `get_testmatrix_kms_lehmer` (KMS, Lehmer).
- `scond` takes the scaling type (`'D'`/`'C'`) second, not the norm type.
