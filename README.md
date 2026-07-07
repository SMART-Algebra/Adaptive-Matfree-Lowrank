# Adaptive, Matrix-Free Low-Rank Approximation — Companion Code

MATLAB code reproducing the numerical experiments in the manuscript
*Adaptive, Matrix-Free Low-Rank Approximation*.

Requirements: MATLAB (R2016b or later — the code uses implicit array expansion).

## Repository layout

```
randQB_MF_Fro.m            proposed method, Algorithm 3.1 (Frobenius norm)
randQB_AF_Fro.m            proposed method, Algorithm 3.3 (adjoint-free)
Prune_Rank_Fro.m           rank pruning, Algorithm 3.2 (used by the two above)
randQB_MF_Spec.m           proposed method, Algorithm 4.1 (spectral norm)
randQB_AF_Spec.m           proposed method, Section 4.3 (adjoint-free)
Prune_Rank_Spec.m          rank pruning, Algorithm 4.2 (used by the two spectral methods)
spec_norm_power.m          power-iteration spectral-norm estimate (shared by the spectral methods)
Generate_Test_Matrix.m     shared test-matrix generator (used by all experiments)
Fro_Norm_Experiments/      Frobenius-norm experiment drivers & baselines  (Section 5.2)
Spec_Norm_Experiments/     spectral-norm experiment drivers & baselines   (Section 5.3)
```

The proposed methods live at the repository root so they can be called directly;
the experiment scripts in the subfolders use them.

The repository root must be on the MATLAB path so that `Generate_Test_Matrix`
is visible from the experiment subfolders (e.g. run `addpath('..')` from within
a subfolder, or run from the repository root).

## Test matrices

`A = Generate_Test_Matrix(m, n, name)` returns `A = U*S*V'` with random
orthonormal `U`, `V` and a prescribed singular-value decay selected by the
string `name`:

| `name`       | singular values            | decay |
|--------------|----------------------------|-------|
| `'exp'`      | `sigma_i = exp(-i/2)`      | fast / exponential |
| `'poly2'`    | `sigma_i = 1/i^2`          | polynomial |
| `'poly4'`    | `sigma_i = 1/i^4`          | polynomial |
| `'s_shape'`  | piecewise (linear then polynomial, in log10) | S-shaped |

In the paper these are Matrix 1 (`'exp'`), Matrix 2 (`'poly4'` for the
Frobenius norm, `'poly2'` for the spectral norm), and Matrix 3 (`'s_shape'`).

---

## Proposed methods (repository root)

| File | Paper | Description |
|------|-------|-------------|
| `randQB_MF_Fro.m`  | Algorithm 3.1 | Matrix-free, adaptive low-rank QB approximation in the Frobenius norm. |
| `randQB_AF_Fro.m`  | Algorithm 3.3 | Adjoint-free variant (uses only the forward operator `A*x`). |
| `Prune_Rank_Fro.m` | Algorithm 3.2 | Rank pruning of the final block, called by the two algorithms above. |
| `randQB_MF_Spec.m` | Algorithm 4.1 | Matrix-free, adaptive low-rank QB approximation in the spectral norm. |
| `randQB_AF_Spec.m` | Section 4.3   | Adjoint-free spectral variant (uses only the forward operator `A*x`). |
| `Prune_Rank_Spec.m`| Algorithm 4.2 | Spectral-norm rank pruning of the final block, called by `randQB_MF_Spec` and `randQB_AF_Spec`. |
| `spec_norm_power.m`| Remark (§4.1) | Power-iteration estimate of `\|M\|_2`, shared by `randQB_MF_Spec` and `randQB_AF_Spec`. |

```matlab
[Q, B, Err_ind, Exact_err] = randQB_MF_Fro(A, b, epsilon, tol_type, compute_exact)
[Q,    Err_ind, Exact_err] = randQB_AF_Fro(A, b, epsilon, tol_type, compute_exact)
s = Prune_Rank_Fro(B, alpha, tau)

[Q, B, Err_ind, Exact_err] = randQB_MF_Spec(A, b, epsilon, tol_type, compute_exact)
[Q,    Err_ind, Exact_err] = randQB_AF_Spec(A, b, epsilon, tol_type, compute_exact)
[U_r, SV_r]     = Prune_Rank_Spec(B, tau)   % SV_r = S_r*V_r'; pruned block is U_r*SV_r
```

The Frobenius and spectral methods share the same `tol_type`
(`'relative'`/`'absolute'`) and `compute_exact` options.  The spectral
methods need no power-iteration parameter: `randQB_MF_Spec` estimates
`||A||_2` from the first block via `||B_1||_2` (Remark 4.1).

- `b` — block size; `epsilon` — tolerance (`epsilon > 0`).
- `tol_type` (optional, default `'relative'`): `'relative'` stops at
  `||A - Q*B||_F <= epsilon*||A||_F`; `'absolute'` stops at
  `||A - Q*B||_F <= epsilon`.
- `compute_exact` (optional, default `false`): if `true`, also returns the
  exact residual norm per iteration. This forms the dense residual and is
  **not** matrix-free — use for diagnostics only.
- The computed rank is `size(Q,2)`.

## Frobenius-norm experiments (`Fro_Norm_Experiments/`)

### Reference methods (baselines)

| File | Description |
|------|-------------|
| `ref_svd_Fro.m`   | Truncated SVD (the optimal rank for a given tolerance). |
| `ref_randQB_EI.m` | `randQB_EI` — an independent implementation of the algorithm of Yu, Gu & Li (*SIMAX* 2018), used as a comparison baseline. |

```matlab
[A_k, k]          = ref_svd_Fro(A, tol, tol_type)
[Q, B, k, EE, fail, neg] = ref_randQB_EI(A, tol, b, P, tol_type)   % P = power-iteration parameter; neg = indicator went negative
```

Both accept the same `'relative'`/`'absolute'` `tol_type` option.

### Comparison engine

```matlab
[ERR_IND, RANKS, REL_ERR, EI_NEG] = compare_Fro_methods(A, b, TOL)
```

Runs all four methods (t-SVD, `randQB_EI`, `randQB_MF_Fro`, `randQB_AF_Fro`)
over a list of tolerances `TOL`. With `ntol = numel(TOL)`:

- `ERR_IND` — `ntol`-by-4 cell of per-iteration error curves; columns are
  `[exact, randQB_EI, randQB_MF_Fro, randQB_AF_Fro]` (column 1 is the exact
  reference, columns 2–4 are each method's error indicator).
- `RANKS`   — `ntol`-by-4 computed ranks, columns `[t-SVD, EI, MF, AF]`.
- `REL_ERR` — `ntol`-by-4 *exact* relative errors `||A-Q*B||_F/||A||_F`,
  columns `[t-SVD, EI, MF, AF]`.
- `EI_NEG`  — `ntol`-by-1 logical; `true` at tolerances where `randQB_EI`'s
  error indicator went negative (Gramian cancellation). The drivers use it to
  hide `randQB_EI` from both plots at those tolerances.

### Drivers (reproduce the Section 5.2 figures and tables)

| File | Test matrix |
|------|-------------|
| `Drive_Fro_exp.m`    | exponential decay (`'exp'`) |
| `Drive_Fro_poly.m`   | polynomial decay (`'poly4'`) |
| `Drive_Fro_sshape.m` | S-shaped decay (`'s_shape'`) |

Each driver sets up its matrix, calls `compare_Fro_methods`, and then:

1. produces the two figure panels (computed rank vs. tolerance, and the
   per-iteration error indicators), saving them as EPS
   (e.g. `rank_Fro_exp.eps`, `err_ind_Fro_exp.eps`); and
2. prints the results table (rank `k` and relative error per method) as
   copy-paste LaTeX rows, with `randQB_EI` failures marked `fail`.

Example:

```matlab
addpath('..');        % so Generate_Test_Matrix is on the path
Drive_Fro_exp         % generates the exponential-decay figure + table
```

---

## Spectral-norm experiments (`Spec_Norm_Experiments/`)

### Reference methods (baselines)

| File | Description |
|------|-------------|
| `ref_svd_Spec.m`   | Truncated SVD (the optimal rank for a given spectral-norm tolerance). |
| `ref_randQB_HMT.m` | `randQB_HMT`, the adaptive randomized range finder of Halko, Martinsson & Tropp (Algorithm 4.2 in *SIAM Review* 2011), with `r` random probe vectors for the a posteriori error estimator. |

```matlab
[A_k, k]      = ref_svd_Spec(A, tol, tol_type)
[Q, k, M]     = ref_randQB_HMT(A, epsilon, r)   % epsilon ABSOLUTE; M = indicator history
```

`ref_svd_Spec` accepts the same `'relative'`/`'absolute'` `tol_type` option;
`ref_randQB_HMT` takes an absolute tolerance (the comparison engine passes
`tol*||A||_2`).

### Comparison engine

```matlab
[ERR_IND, RANKS, REL_ERR] = compare_Spec_methods(A, b, TOL)
```

Runs all four methods (t-SVD, `randQB_HMT`, `randQB_MF_Spec`,
`randQB_AF_Spec`) over a list of tolerances `TOL`.  With `ntol = numel(TOL)`:

- `ERR_IND` — `ntol`-by-4 cell of per-iteration error curves; columns are
  `[exact, randQB_HMT, randQB_MF_Spec, randQB_AF_Spec]` (column 1 is the exact
  reference, columns 2–4 are each method's error indicator).
- `RANKS`   — `ntol`-by-4 computed ranks, columns `[t-SVD, HMT, MF, AF]`.
- `REL_ERR` — `ntol`-by-4 *exact* relative errors `||A-Q*B||_2/||A||_2`,
  columns `[t-SVD, HMT, MF, AF]`.

### Drivers (reproduce the Section 5.3 figures and tables)

| File | Test matrix |
|------|-------------|
| `Drive_Spec_exp.m`    | exponential decay (`'exp'`, `b = 16`) |
| `Drive_Spec_poly.m`   | polynomial decay (`'poly2'`, `b = 32`) |
| `Drive_Spec_sshape.m` | S-shaped decay (`'s_shape'`, `b = 32`) |

Each driver sets up its matrix, calls `compare_Spec_methods`, and then:

1. produces the two figure panels (computed rank vs. tolerance, and the
   per-iteration error indicators), saving them as EPS
   (e.g. `rank_Spec_exp.eps`, `err_ind_Spec_exp.eps`); and
2. prints the results table (rank `k` and relative error per method) as
   copy-paste LaTeX rows. `randQB_HMT` does not fail — it simply overestimates
   the rank (sometimes up to full rank) — so the raw numbers are shown for all
   methods.

Example:

```matlab
addpath('..');        % so Generate_Test_Matrix and the methods are on the path
Drive_Spec_exp        % generates the exponential-decay figure + table
```

---

## License

This code is released under the MIT License; see [`LICENSE`](LICENSE).

## Acknowledgment

AI-based assistants were used to help develop and document this code. The
authors reviewed and tested all code and take full responsibility for it.
