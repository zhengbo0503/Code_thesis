# Complex-arithmetic experiments

## 1. Description

This folder contains the MATLAB codes for the complex-arithmetic experiments
of our mixed-precision one-sided Jacobi SVD algorithm. We present three
numerical experiments to assess the relative forward accuracy and the speed
of the algorithms on complex matrices.

The first experiment uses matrices with varying condition numbers. The second
experiment uses matrices with varying numbers of columns. Finally, we perform
an accuracy and timing experiment for the `(single, single, double)` version
of our algorithm. 

## 2. Requirements

The codes were tested in the following configuration:

```text
OS: macOS 26.6.2 25G83 arm64
CPU: Apple M3 Pro
Memory: 36 GB
MATLAB version: 25.2.0.3042426 (R2025b) Update 1
Advanpix version: 5.4.4 Build 16174
```

We use the Advanpix Multiprecision Computing Toolbox to simulate quadruple
and octuple precision. The four complex LAPACK routines `ZGESVJ`, `ZGEJSV`,
`CGESVJ`, and `CGEJSV` are called through MEX files linked with the
`mwlapack` library supplied by MATLAB R2025b.

## 3. Codes description

- `mposj_complex.m`: complex implementation of MP3JacobiSVD with the
  `(single, double, quadruple)` setting.
- `mposj_ssd_complex.m`: complex implementation with the
  `(single, single, double)` setting.
- `complex_randsvd.m`: generates complex random matrices for the experiments.
- `test_complex_varying_kappa.m`: accuracy test for varying condition
  numbers.
- `test_complex_varying_n.m`: accuracy test for varying numbers of columns.
- `test_complex_ssd_timing.m`: accuracy and timing test for the SSD version.
- `plot_complex_results.m`: produces the figures from the computed results.
- `run_all_complex_full.m`: runs all three experiments and saves the results.

## 4. Build the MEX files

Open the `complex_experiments` folder in MATLAB R2025b and run:

```matlab
setup_complex_paths;
build_all_complex_mex;
validate_complex_stack;
```

The last command performs a short check of the four MEX files and the two
versions of our algorithm. If all checks are passed, MATLAB prints:

```text
Complex preflight validation passed under MATLAB 2025b.
```

## 5. Quick example

To run a quick example of the complex MP3JacobiSVD algorithm, use the
following code in the MATLAB command window:

```matlab
clc;
m = 100; n = 50;
A = complex_randsvd(m,n,1e8,3,'double',1);
[U,S,V] = mposj_complex(A);

fprintf('Backward error: %e\n', ...
    norm(A-U*S*V','fro')/norm(A,'fro'));
fprintf('Orthogonality of U: %e\n', ...
    norm(U'*U-eye(n),'fro'));
fprintf('Orthogonality of V: %e\n', ...
    norm(V'*V-eye(n),'fro'));
```

## 6. Run all experiments

For a single run on macOS, open `run_all_complex_full.m` and click **Run**, or
enter the following command:

```matlab
run_all_complex_full
```

This file first checks the MEX files and the algorithms, and then runs the
three experiments in sequence. It also keeps the computer awake during the
calculation. The complete run took about 2 hours and 36 minutes on the test
computer described above.

The three experiments can also be run separately:

```matlab
test_complex_varying_kappa;
test_complex_varying_n;
test_complex_ssd_timing;
```

## 7. Experiment settings

To ensure reproducibility, we use fixed random seeds and save them in the
result files. The five positive `MODE` values in `complex_randsvd` give the
same singular value distributions as MATLAB `gallery('randsvd')`. The left
and right singular-vector factors are generated from complex Gaussian
matrices using QR factorizations. Therefore, the singular value distributions
are the same, but the generated matrices are not exactly the same as those
from `gallery('randsvd')`. For `MODE = 5`, the input condition number is an
upper limit. The actual condition number is saved and used in the figure.

For the first two experiments, we test the following four algorithms:

- MP3JacobiSVD;
- `ZGESVJ`, the LAPACK one-sided Jacobi algorithm;
- `ZGEJSV`, the LAPACK preconditioned one-sided Jacobi algorithm;
- MATLAB `svd`.

The working precision is double precision. All the reference singular values
are computed using MATLAB `svd` at 71 decimal digits with Advanpix. We use
`sqrt(m*n)*u*scond(At)` as a practical reference line, where `At` is the
preconditioned matrix computed before conversion to the working precision.

The settings of the three experiments are:

1. **Varying matrix condition number.** We generate complex matrices with
   `m = 1000` and `n = 800`. The condition number takes 20 logarithmically
   spaced values from `1e3` to `1e15`, and `MODE = 1,...,5`.
2. **Varying matrix size.** We fix `m = 1000` and the condition number to
   `1e8`. The number of columns takes 15 logarithmically spaced values from
   `10` to `1000`, and `MODE = 1,...,5`.
3. **SSD accuracy and timing.** In a separate experiment, we use complex
   single-precision matrices with `m = 1000`, condition number `1e6`, and
   `MODE = 3`. The number of columns takes 10 logarithmically spaced values
   from `100` to `1000`. We compare the SSD version of our algorithm with
   `CGESVJ` and `CGEJSV`. For each matrix size, the first timing run is
   discarded, and the average of the second and third runs is reported. The
   reference singular values and `scond(At)` are computed outside the timing.

As in Figure 5 of the paper, the SSD figure shows both the maximum relative
forward error and the runtime. All the figures use logarithmic scales.

## 8. Output files

The numerical results are saved in the `data` folder:

- `complex_varying_kappa_full.mat`;
- `complex_varying_n_full.mat`;
- `complex_ssd_timing_full.mat`.

Each experiment also saves its figure in both PNG and PDF formats in the
`plots` folder. For example, a saved result can be plotted again using:

```matlab
load('data/complex_varying_kappa_full.mat','results');
plot_complex_results(results,'varying_kappa');
```

## 9. Generative AI disclosure

The code in this folder was generated 
using OpenAI GPT-5.6 Sol through Codex, 
based on my own original code.
I reviewed and tested the generated code.
