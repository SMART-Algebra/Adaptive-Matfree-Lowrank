function [Q, B, Err_ind, Exact_err] = randQB_MF_Fro(A, b, epsilon, tol_type, compute_exact)
% randQB_MF_Fro  Matrix-free, adaptive low-rank QB approximation in the
%                Frobenius norm (Algorithm 3.1).
%
%   [Q, B, Err_ind, Exact_err] = randQB_MF_Fro(A, b, epsilon)
%   [...] = randQB_MF_Fro(A, b, epsilon, tol_type)
%   [...] = randQB_MF_Fro(A, b, epsilon, tol_type, compute_exact)
%
%   Builds Q (orthonormal columns) and B = Q'*A block by block until the
%   randomized indicator certifies the target tolerance, then prunes the
%   final block to a near-minimal rank via Prune_Rank_Fro.  A is accessed
%   only through matrix-vector products.  The computed rank is size(Q,2).
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       tolerance (epsilon > 0)
%     tol_type      (optional, default 'relative')
%                     'relative' : stop when ||A-Q*B||_F <= epsilon*||A||_F;
%                                  ||A||_F is estimated on the fly.
%                     'absolute' : stop when ||A-Q*B||_F <= epsilon.
%     compute_exact (optional, default false) if true, also return the exact
%                   residual norm at each iteration.  This forms the dense
%                   residual A - Q*B and is NOT matrix-free; use only for
%                   diagnostics/experiments.
%
%   Outputs
%     Q, B          QB factors satisfying the requested tolerance
%     Err_ind       randomized error indicator per iteration
%     Exact_err     exact residual norm per iteration ([] unless compute_exact)

    if nargin < 4 || isempty(tol_type)
        tol_type = 'relative';
    end
    if nargin < 5 || isempty(compute_exact)
        compute_exact = false;
    end
    use_relative = strcmpi(tol_type, 'relative');

    [m, n]  = size(A);
    min_mn  = min(m, n);
    maxit   = ceil(min_mn / b);

    % Preallocate factors and per-iteration diagnostics.
    Q = zeros(m, maxit * b);
    B = zeros(maxit * b, n);
    ncol = 0;                       % number of columns currently stored in Q

    Err_ind = zeros(maxit, 1);
    if compute_exact
        Exact_err = zeros(maxit, 1);
    else
        Exact_err = [];
    end

    sum_Z = 0;                      % running sqrt(sum_k ||Z_k||_F^2) (relative only)

    for i = 1:maxit

        % Random sketch of A.  Folding 1/sqrt(b) into Omega makes ||Y||_F a
        % direct (unbiased) estimate of the residual Frobenius norm.
        Omega = (1/sqrt(b)) * randn(n, b);
        Z = A * Omega;

        % Target tolerance for this iteration.
        if use_relative
            sum_Z = hypot(sum_Z, norm(Z, 'fro'));
            tau   = epsilon * (sum_Z / sqrt(i));   % epsilon * (running ||A||_F estimate)
        else
            tau   = epsilon;
        end

        % Implicitly sketch the residual:  Y = (A - Q*B)*Omega.
        if ncol == 0
            Y = Z;
        else
            Qc = Q(:, 1:ncol);
            Bc = B(1:ncol, :);
            Y  = Z - Qc * (Bc * Omega);
        end

        if compute_exact
            if ncol == 0
                Exact_err(i) = norm(A, 'fro');
            else
                Exact_err(i) = norm(A - Qc * Bc, 'fro');
            end
        end

        % Randomized stopping criterion.
        Err_ind(i) = norm(Y, 'fro');
        if Err_ind(i) < tau
            % Prune the final block to a near-minimal rank (Algorithm 3.2).
            keep  = ncol - b;               % columns kept from the earlier blocks
            Q_end = Q(:, keep+1:ncol);      % last block Q_{i-1}
            B_end = B(keep+1:ncol, :);      % last block B_{i-1}
            s = Prune_Rank_Fro(B_end, Err_ind(i), tau);

            % s indexes WITHIN the last block (a subset of 1:b), not the global
            % column/row numbering.  Keep the earlier blocks whole and append
            % only the s-selected columns of Q_end / rows of B_end.
            Q = [Q(:, 1:keep), Q_end(:, s)];
            B = [B(1:keep, :);  B_end(s, :)];

            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end
            return
        end

        % Orthonormalize the new block, then re-orthogonalize against Q.
        [Qi, ~] = qr(Y, 'econ');
        if ncol > 0
            Qc = Q(:, 1:ncol);
            [Qi, ~] = qr(Qi - Qc * (Qc' * Qi), 'econ');
        end
        Bi = Qi' * A;

        Q(:, ncol+1:ncol+b) = Qi;
        B(ncol+1:ncol+b, :) = Bi;
        ncol = ncol + b;
    end

    % Tolerance never met: return the full factorization.
    Q = Q(:, 1:ncol);
    B = B(1:ncol, :);
end
