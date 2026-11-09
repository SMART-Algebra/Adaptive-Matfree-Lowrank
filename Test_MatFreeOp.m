% Test_MatFreeOp  Checks MatFreeOp against the dense matrix (run from the repository root).
%
%   1. The operations the routines use (size, A*X, X*A, A') agree with the matrix.
%   2. With the same random seed, every randQB_* routine returns the same
%      factors from the wrapped operator as from the dense matrix.
%   3. A handle that accepts one column at a time can be wrapped by a loop.
%   4. An operator defined only through sparse solves (the Jacobian-style
%      use case) runs end to end into QB_to_SVD; the error is checked
%      against a dense copy formed for the check only.

rng(0);
m = 600;  n = 500;  b = 16;  tol = 1e-6;
A_dense = Generate_Test_Matrix(m, n, 'exp') + 0.3i * Generate_Test_Matrix(m, n, 'exp');
A_op    = MatFreeOp(m, n, @(X) A_dense * X, @(Y) A_dense' * Y);

% 1. The operations the routines use.
[mm, nn] = size(A_op);
assert(mm == m && nn == n && size(A_op, 1) == m && isequal(size(A_op), [m n]), 'size');
X = randn(n, 3);  Y = randn(m, 3);
assert(isequal(A_op * X, A_dense * X), 'A*X');
assert(isequal(A_op' * Y, A_dense' * Y), 'A''*Y');
assert(norm(Y' * A_op - Y' * A_dense, 'fro') <= 1e-14 * norm(Y' * A_dense, 'fro'), 'X*A');

% 2. Same factors as the dense matrix, routine by routine (same seed => same Omega).
routines = {@randQB_MF_Fro, @randQB_AF_Fro, @randQB_MF_Spec, @randQB_AF_Spec};
for j = 1:numel(routines)
    f = routines{j};
    rng(1);  [F1, E1] = f(A_dense, b, tol);     % Q and B (MF) or Q and Err_ind (AF)
    rng(1);  [F2, E2] = f(A_op,    b, tol);
    % Same rank, and the same approximation.  The basis Q itself may differ
    % at ~1e-10 near the stopping tolerance, because the MF routines form
    % Qi'*A as (A'*Qi)' for an operator and the rounding difference is
    % amplified by the conditioning of the last residual sketch.
    assert(isequal(size(F1), size(F2)) && isequal(size(E1), size(E2)), ...
           '%s: rank differs between matrix and operator', func2str(f));
    assert(norm(F1 - F2, 'fro') <= 1e-8 * norm(F1, 'fro'), ...
           '%s: Q differs between matrix and operator', func2str(f));
    assert(norm(E1(:) - E2(:)) <= 1e-12 * norm(E1(:)), ...
           '%s: second output differs between matrix and operator', func2str(f));
    if contains(func2str(f), 'MF')        % E = B: compare the approximations Q*B
        assert(norm(F1*E1 - F2*E2, 'fro') <= 1e-12 * norm(A_dense, 'fro'), ...
               '%s: Q*B differs between matrix and operator', func2str(f));
    end
end

% 3. A single-column handle, wrapped by a loop over the columns of X.
Avec   = @(x) A_dense * x;                      % accepts one column only
A_loop = MatFreeOp(m, n, @(X) apply_columns(Avec, X));
assert(isequal(A_loop * X, A_dense * X), 'column loop');
rng(1);  Q_loop = randQB_AF_Fro(A_loop, b, tol);   % adjoint-free routine: no A' needed
rng(1);  Q_ref  = randQB_AF_Fro(A_dense, b, tol);
assert(norm(Q_loop - Q_ref, 'fro') <= 1e-10 * norm(Q_ref, 'fro'), 'column loop in randQB_AF_Fro');

% 4. Solve-based operator in the Jacobian pattern (solves, then a pointwise
%    scaling; the adjoint scales first and then solves with the transposed
%    factorization):  J = K^{-2} * diag(u)  with K the sparse 2D Laplacian.
%    Two solves per product give singular values that decay like 1/j^2, so
%    the discovered rank is well below the dimension at this tolerance.
N  = 30;  K = delsq(numgrid('S', N + 2));        % N^2-by-N^2 sparse SPD
nK = size(K, 1);
u  = exp(1i * linspace(0, 2*pi, nK)');           % a complex "incident field"
dK = decomposition(K, 'chol');                   % factorized once, captured by the handles
J_op = MatFreeOp(nK, nK, @(X) dK \ (dK \ (u .* X)), @(Y) conj(u) .* (dK' \ (dK' \ Y)));
tolJ = 1e-3;
[Q, B]        = randQB_MF_Fro(J_op, b, tolJ);
[U, S, V]     = QB_to_SVD(Q, B);
Kf            = full(K);
J_dense       = Kf \ (Kf \ diag(u));            % for the check only
rel_err       = norm(J_dense - U*S*V', 'fro') / norm(J_dense, 'fro');
assert(rel_err <= 2 * tolJ, 'solve-based operator: error %.2e exceeds 2*tol', rel_err);
assert(size(S, 1) < nK / 2, 'solve-based operator: rank %d is not a reduction', size(S, 1));
P = V * (S \ U');                                % pseudoinverse of the approximation, factored
assert(norm(P * (U*S*V') * P - P, 'fro') <= 1e-8 * norm(P, 'fro'), 'pseudoinverse identity');
fprintf('Test_MatFreeOp: solve-based operator of size %d, discovered rank %d, relative error %.2e (tol %.0e)\n', ...
        nK, size(S, 1), rel_err, tolJ);
disp('Test_MatFreeOp: all checks passed');

function Y = apply_columns(fvec, X)
% Apply a single-column handle to every column of X.
    p = size(X, 2);
    y1 = fvec(X(:, 1));
    Y  = zeros(size(y1, 1), p, 'like', y1);
    Y(:, 1) = y1;
    for j = 2:p
        Y(:, j) = fvec(X(:, j));
    end
end
