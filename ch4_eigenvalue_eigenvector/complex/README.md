# Complex two-sided Jacobi accuracy experiments

These are complex Hermitian adaptations of the real two-sided functions, with experiment structure copied from the one-sided complex tests. There are two accuracy experiments, without timing measurements.

## Run

From the Code_thesis root in MATLAB:

```matlab
addpath('/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac')
addpath('ch4_eigenvalue_eigenvector/complex')
resultsK = test_complex_evals_varying_kappa();
resultsN = test_complex_evals_varying_n();
```

The top-level `run_all_experiments` includes both when `IncludeComplex` is true (the default). To select just these two:

```matlab
run_all_experiments('Only', ...
    {'test_complex_evals_varying_kappa','test_complex_evals_varying_n'})
```

Both functions print `[START]`, progress after each matrix, and `[END]` after output generation. They do not clear the Command Window. Full sweeps can take several hours; no complete runtime has been measured. Each matrix prints a starting message before computation and a completion message afterward, so a long individual case can be identified. Results are saved after the whole grid; there is no per-case checkpoint/restart mechanism.

## Experiments and correspondence

| New function | One-sided template | Full parameters |
| --- | --- | --- |
| `test_complex_evals_varying_kappa` | `ch5_singular_value_singular_vector/complex/test_complex_varying_kappa.m` | Square order n=800; `logspace(3,15,20)`; modes 1:5; 100 matrices |
| `test_complex_evals_varying_n` | `ch5_singular_value_singular_vector/complex/test_complex_varying_n.m` | `unique(round(logspace(1,3,15)))`; target kappa=1e8; modes 1:5; 75 matrices |

The one-sided tests' mode-major loops, seed formulas (100000/200000 + 1000*mode + case index), precision settings, metric arrays, MAT output, and five-panel plots are retained. The SVD-specific rectangular dimension, second unitary factor, QR/direct branch and LAPACK-specific outputs are replaced by the corresponding EVD quantities. Method order is MP3Jacobi, Jacobi, MP2Jacobi, MATLAB eig. Each Jacobi method records INFO, sweeps and rotations. Both mixed-precision methods use the same preconditioner construction as their real counterparts.

`complex_randhpd` follows the one-sided `complex_randsvd` generator: the same five spectra, private random stream and phase-normalised Gaussian QR. For a Hermitian eigenproblem it uses a single square Q and forms `Q*diag(lambda)*Q'`, then stores `(A+A')/2`. It does not reproduce the one-sided matrix entries for the same seed. Mode 5's prescribed condition number can be smaller than the requested kappa. No diagonal shift is applied if roundoff destroys positive definiteness; the stored-input reference check stops the experiment instead.

## Minimal changes to the real functions

| New function | Real two-sided source | Arithmetic changes |
| --- | --- | --- |
| `jacobi_pair_complex` | `../shared/jacobi_pair.m` | Compute real tau using `abs(A(p,q))`, then attach the phase `A(p,q)/abs(A(p,q))` to s. Keep the original t and c formulas. |
| `cjacobi_complex` | `../shared/cjacobi.m` | Check Hermitian input; use `J=[c,s;-conj(s),c]`; make updated pivot diagonals real. Keep the row/column updates, cyclic order, current-diagonal stopping test, 30-sweep cap and sorting. |
| `mp_pjacobi_complex` | `../shared/mp_pjacobi.m` | Call the complex kernel and preserve Hermitian structure of the high-precision product. Keep low-precision eig, double QR, MP2/MP3 multiplication, scaling and reconstruction. Expose existing INFO as a seventh output; unavailable MP2 BOUND/SCOND are NaN. |

All original solver arithmetic and the one-sided functions remain unchanged. The real and complex two-sided solver docstrings now describe accumulated orthogonality/unitarity error without claiming exact equality to unit roundoff; SCOND is explicitly described as diagonally scaled conditioning. The complex wrapper uses the existing `float_params` and **explicit diagonal scaling** `scond(Athcomp,'D')`; it does not call the real-only column-scaling branch.

MP3 multiplication uses object-local 34-digit Advanpix values, and references use `eig(mp(A,71))` on the **actual stored double matrix**. The prescribed eigenvalues are not used as ground truth. The global multiprecision setting is not modified. The scaled-condition diagnostic uses the wrapper's 34-digit preconditioned matrix; the original real experiment separately recomputes that diagnostic at 71 digits, so the formula agrees but the diagnostic precision differs.

The retained wrapper BOUND formula is `7*n*u*SCOND + sqrt(n)*u`. The plots and `reference_at` use its leading term `7*n*u*SCOND`, matching the real eigenvalue plots. This is labelled a **real-case reference** in the complex experiment, not asserted as a proved complex-rounding bound. A curve exceeding it does not trigger a theoretical-bound failure. MP2/MP3 nonconvergence still fails the experiment. Plain Jacobi reaching its 30-sweep cap produces a warning; its final-iterate errors and INFO=-1 are retained so the comparison baseline cannot discard the complete sweep. Such points receive black x overlays with a separate legend entry in the PNG/PDF plots. Future thesis captions using these data must identify those points as results that did not satisfy the stopping criterion.

## Files and exact CSV schema

Each successful full sweep writes beneath this folder:

| Experiment | CSV | Other files with the same stem |
| --- | --- | --- |
| Varying condition number | `data/complex_evals_varying_kappa_full.csv` | `data/*.mat`, `plots/*.png`, `plots/*.pdf` |
| Varying order | `data/complex_evals_varying_n_full.csv` | `data/*.mat`, `plots/*.png`, `plots/*.pdf` |

Each CSV has one header row and 29 numeric columns. Rows are ordered by mode, then case index: 100 rows for varying kappa, 75 for varying n. Values are comma-delimited with `%.17g` precision; replacement uses a temporary file and rename. The CSVs are native outputs of these new experiments, separate from the 43 files currently consumed by the thesis.

```text
mode,seed,n,kappa_target,kappa_prescribed,kappa_actual,forward_mp3,forward_jacobi,forward_mp2,forward_eig,residual_mp3,residual_jacobi,residual_mp2,residual_eig,unitarity_mp3,unitarity_jacobi,unitarity_mp2,unitarity_eig,info_mp3,info_jacobi,info_mp2,sweeps_mp3,sweeps_jacobi,sweeps_mp2,rotations_mp3,rotations_jacobi,rotations_mp2,scond_at,reference_at
```

| Columns | Meaning |
| --- | --- |
| `mode,seed,n` | Spectrum mode, private generator seed, square matrix order |
| `kappa_target` | Generator input kappa |
| `kappa_prescribed` | max(lambda)/min(lambda) before storing the matrix |
| `kappa_actual` | max(dref)/min(dref) from the 71-digit stored-input reference |
| `forward_*` | max_i abs(dref_i - d_i)/abs(dref_i), descending eigenvalues |
| `residual_*` | norm(A - V*diag(d)*V','fro')/norm(A,'fro') |
| `unitarity_*` | norm(V'*V - eye(n),inf) |
| `info_*` | Jacobi convergence status: 0 successful, -1 at the retained sweep limit; plain Jacobi may be -1 in exported results; MP2/MP3 must be 0 |
| `sweeps_*,rotations_*` | Counts from each Jacobi kernel; no invented counts for MATLAB eig |
| `scond_at` | MP3 preconditioned matrix's diagonally scaled condition number, computed at high precision and exported as double |
| `reference_at` | 7*n*(eps('double')/2)*scond_at; empirical complex reference |

The MAT structure additionally preserves methods, grids, precision settings, completion flags and MATLAB release. `passed` means the experiment collected its required results, allowing a flagged plain-Jacobi final iterate; `allJacobiConverged` separately states whether all three Jacobi methods reported INFO=0 for every case. A saved MAT may precede CSV/plot export, so verify the final END marker and output files as well. The plotting routine reads this structure as in the one-sided template. The CSV includes every plotted coordinate/error/reference for subsequent TikZ use.

## Small validation

```matlab
addpath('/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac')
addpath('ch4_eigenvalue_eigenvector/complex')
report = validate_complex_evd();
```

This checks arbitrary rotation phases (including equal diagonals), diagonal/scalar matrices, all five spectra at small order through target kappa=1e15, real-input agreement with the original functions, the original scaling paths, private random-state preservation, global precision preservation, and both reduced experiment/output paths. It saves temporary MAT/CSV/PNG/PDF outputs and reports their directory. CSV numeric fields are checked against the in-memory arrays after re-reading. For moderate-condition cases the validator checks MP3 relative accuracy against the 71-digit reference. At kappa=1e15 it checks convergence, residual and unitarity, and records relative accuracy and scaled condition in `report.endpoint`: a fixed near-roundoff forward-error threshold is inappropriate when the preconditioned matrix remains ill-conditioned.

For just one reduced experiment:

```matlab
results = test_complex_evals_varying_kappa('smoke',tempname);
results = test_complex_evals_varying_n('smoke',tempname);
```

The smoke profile retains all five modes, with n=8 and kappa=[1e3,1e8] for the condition sweep, and n=[6,10] at kappa=1e8 for the size sweep. Smoke output names end in `_smoke` and cannot overwrite full-result filenames.

## Large-grid convergence probe

The smoke validator cannot establish behavior at n=800. Before the full rerun, also use:

```matlab
large = probe_complex_evd_convergence();
```

This runs only plain Jacobi on the actual condition-grid endpoint n=800, kappa=1e15, mode=3, caseIndex=20, seed=103020. It records INFO, rotations, sweeps, relative Frobenius residual and infinity-norm unitarity loss in a temporary MAT report. It accepts INFO=-1 as a recorded baseline outcome, while checking finite diagnostics and residual/unitarity below `100*n*u`. It does not compute a high-precision eigenvalue reference, certify forward accuracy, or establish that more sweeps could never converge. The core still uses its original 30-sweep cap. A passed probe means the diagnostic checks passed; inspect `converged` and `info` separately.
