# Chapter 3 — Preconditioning the Jacobi Algorithm

Code for `chaps/preconditioning.tex`, section *Numerical experiments*.
Originally from [`zhengbo0503/Code_htwz25`](https://github.com/zhengbo0503/Code_htwz25),
directory `MATLAB/` ([arXiv:2501.03742](https://arxiv.org/abs/2501.03742)).

## Experiments

| Thesis     | Label                       | Script                                 | Writes                           |
| ---------- | --------------------------- | -------------------------------------- | -------------------------------- |
| Figure 3.1 | `fig.prec-offAt`            | `experiments/test_off_quantity.m`      | `output/offAt_mode1–5.csv`       |
| Figure 3.2 | `fig.prec-reduct-cond`      | `experiments/test_scond_reduction.m`   | `output/reducekappa_mode1–5.csv` |
| Table 3.1  | `tab.prec-anymatrix-reduct` | none — values hard-coded in the thesis | —                                |

## Running the tests

```bash
tools/run_smoke.sh                    # smoke_test.m, R2026a by default
```

```matlab
cd ch3_preconditioning/experiments    % scripts assume this working directory
test_off_quantity
test_scond_reduction
```

- Advanpix is required and not on the path by default:
  `addpath('/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac')`.
- Score headless runs by the printed `N/M passed` line, not the exit code.

## Renames from the original repo

| Was                            | Now                                  |
| ------------------------------ | ------------------------------------ |
| `MATLAB/test/test_offA.m`      | `experiments/test_off_quantity.m`    |
| `MATLAB/test/test_reduction.m` | `experiments/test_scond_reduction.m` |
| `MATLAB/src/`, `src/private/`  | `shared/` (flattened)                |

`scond(A, type, nrmtype)` now takes the **scaling type** (`'D'`/`'C'`) second, not
the norm type.
