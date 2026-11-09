% Test_QB_to_SVD  Self-check for QB_to_SVD (run from the repository root).
%
%   Checks, on a real and on a complex test matrix:
%     1. no truncation reproduces Q*B exactly, with orthonormal U and V;
%     2. 'fro' truncation respects its budget and the quadrature identity
%        ||A - U*S*V'||_F^2 = ||A - Q*B||_F^2 + norm(sv_drop)^2;
%     3. 'spec' and 'rank' truncation keep what they promise;
%     4. the factored pseudoinverse V*(S\U') matches pinv of the approximation.

rng(0);
m = 600;  n = 500;  b = 16;  rel_tol = 1e-6;

A_real = Generate_Test_Matrix(m, n, 'exp');
A_cplx = A_real + 0.5i * Generate_Test_Matrix(m, n, 'poly4');

for which = {'real', 'complex'}
    if strcmp(which{1}, 'real'), A = A_real; else, A = A_cplx; end

    [Q, B] = randQB_MF_Fro(A, b, rel_tol);
    k      = size(Q, 2);
    err_QB = norm(A - Q*B, 'fro');

    % 1. No truncation.
    [U, S, V, sv_drop] = QB_to_SVD(Q, B);
    assert(isempty(sv_drop) && size(S, 1) == k, 'rank changed without truncation');
    assert(norm(U'*U - eye(k)) < 1e-12 && norm(V'*V - eye(k)) < 1e-12, 'U or V not orthonormal');
    assert(norm(Q*B - U*S*V', 'fro') <= 1e-12 * norm(B, 'fro'), 'U*S*V'' ~= Q*B');

    % 2. 'fro' truncation: budget and quadrature identity.
    tau = err_QB;                                   % spend a second budget equal to the QB error
    [U, S, V, sv_drop] = QB_to_SVD(Q, B, 'fro', tau);
    r_fro   = size(S, 1);
    err_SVD = norm(A - U*S*V', 'fro');
    assert(norm(sv_drop) <= tau, '''fro'' truncation exceeded its budget');
    assert(abs(err_SVD^2 - (err_QB^2 + norm(sv_drop)^2)) <= 1e-8 * err_SVD^2, ...
           'quadrature identity violated');
    assert(err_SVD <= sqrt(2) * err_QB + 1e-15, 'truncated error larger than the combined budget');

    % 3. 'spec' and 'rank'.
    [~, S, ~, sv_drop] = QB_to_SVD(Q, B, 'spec', 1e-4);
    assert(all(diag(S) > 1e-4) && (isempty(sv_drop) || max(sv_drop) <= 1e-4), '''spec'' truncation wrong');
    [U, S, V, sv_drop] = QB_to_SVD(Q, B, 'rank', 10);
    assert(isequal(size(S), [10 10]) && numel(sv_drop) == k - 10, '''rank'' truncation wrong');

    % 4. Factored pseudoinverse against pinv of the dense approximation.
    P_fac = V * (S \ U');
    P_ref = pinv(U*S*V');
    assert(norm(P_fac - P_ref, 'fro') <= 1e-8 * norm(P_ref, 'fro'), 'pseudoinverse mismatch');

    fprintf('Test_QB_to_SVD (%s): k = %d, ||A-QB||_F/||A||_F = %.2e, rank after ''fro'' truncation = %d: passed\n', ...
            which{1}, k, err_QB / norm(A, 'fro'), r_fro);
end
disp('Test_QB_to_SVD: all checks passed');
