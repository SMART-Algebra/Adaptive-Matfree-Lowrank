function [Q, B, Err_ind, Exact_err] = randQB_MF_Spec(A, b, epsilon, tol_type, compute_exact)
% randQB_MF_Spec  Matrix-free, adaptive low-rank QB approximation in the
%                 spectral norm (Algorithm 4.1).
%
%   [Q, B, Err_ind, Exact_err] = randQB_MF_Spec(A, b, epsilon)
%   [...] = randQB_MF_Spec(A, b, epsilon, tol_type)
%   [...] = randQB_MF_Spec(A, b, epsilon, tol_type, compute_exact)
%
%   Builds Q (orthonormal columns) and B = Q'*A block by block.  Each new
%   block B_i serves as a "look-ahead" error indicator: ||B_i||_2 estimates
%   the residual norm ||A - Q_{1..i-1}B_{1..i-1}||_2 (eq. (4.x)).  Once
%   ||B_i||_2 <= tau, this look-ahead block is discarded and the previous
%   block is truncated to a near-minimal rank via Prune_Rank_Spec (Algorithm
%   4.2).  A is accessed only through matrix-vector products.  The computed
%   rank is size(Q,2).
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       tolerance (0 < epsilon < 1)
%     tol_type      (optional, default 'relative')
%                     'relative' : stop when ||B_i||_2 <= epsilon*||A||_2,
%                                  with ||A||_2 estimated by ||B_1||_2
%                                  (Remark 4.3).
%                     'absolute' : stop when ||B_i||_2 <= epsilon.
%     compute_exact (optional, default false) if true, also return the exact
%                   residual spectral norm of the running approximation.  This
%                   forms the dense residual and is NOT matrix-free; use for
%                   diagnostics only.
%
%   Outputs
%     Q, B          QB factors meeting the requested tolerance (w.h.p.)
%     Err_ind       look-ahead indicator ||B_i||_2 per iteration
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
    maxit   = ceil(min_mn / b) + 1;     % +1 for the final look-ahead block

    % Preallocate factors and per-iteration diagnostics.
    Q = zeros(m, maxit * b);
    B = zeros(maxit * b, n);
    ncol = 0;                           % columns of Q stored so far (blocks 1..i-1)

    Err_ind = zeros(maxit, 1);
    if compute_exact
        Exact_err = zeros(maxit, 1);
    else
        Exact_err = [];
    end

    norm_A_est = [];                    % set from ||B_1||_2 on the first iteration

    for i = 1:maxit

        Omega = randn(n, b);
        Z = A * Omega;

        % Implicit residual sketch  Y = (A - Q*B)*Omega.
        if ncol == 0
            Y = Z;
        else
            Qc = Q(:, 1:ncol);
            Bc = B(1:ncol, :);
            Y  = Z - Qc * (Bc * Omega);
        end

        if compute_exact
            if ncol == 0
                Exact_err(i) = norm(A);
            else
                Exact_err(i) = norm(A - Qc * Bc);
            end
        end

        % Orthonormalize the new block, then re-orthogonalize against Q.
        [Qi, ~] = qr(Y, 'econ');
        if ncol > 0
            Qc = Q(:, 1:ncol);
            [Qi, ~] = qr(Qi - Qc * (Qc' * Qi), 'econ');
        end
        Bi = Qi' * A;

        % Look-ahead indicator  ||B_i||_2 ~ ||A - Q_{1..i-1}B_{1..i-1}||_2.
        Err_ind(i) = norm(Bi);

        % Target tolerance (relative uses ||B_1||_2 as the ||A||_2 estimate).
        if use_relative
            if i == 1
                norm_A_est = Err_ind(1);
            end
            tau = epsilon * norm_A_est;
        else
            tau = epsilon;
        end

        if Err_ind(i) <= tau
            % Discard this look-ahead block; prune the previous block (i-1).
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end

            if ncol == 0
                % A already meets the tolerance: rank-0 approximation.
                Q = zeros(m, 0);
                B = zeros(0, n);
                return
            end

            keep   = ncol - b;                  % columns of blocks 1..i-2
            Q_prev = Q(:, keep+1:ncol);         % block i-1
            B_prev = B(keep+1:ncol, :);         % block i-1
            [U_r, S_r, V_r] = Prune_Rank_Spec(B_prev, tau);

            Q = [Q(:, 1:keep), Q_prev * U_r];
            B = [B(1:keep, :);  S_r * V_r'];
            return
        end

        % Append the block and continue.
        Q(:, ncol+1:ncol+b) = Qi;
        B(ncol+1:ncol+b, :) = Bi;
        ncol = ncol + b;

        if ncol >= min_mn
            % Tolerance never met: return the full factorization.
            Q = Q(:, 1:ncol);
            B = B(1:ncol, :);
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end
            return
        end
    end

    % Loop exhausted (degenerate): return what we have.
    Q = Q(:, 1:ncol);
    B = B(1:ncol, :);
    Err_ind = Err_ind(1:i);
    if compute_exact
        Exact_err = Exact_err(1:i);
    end
end
