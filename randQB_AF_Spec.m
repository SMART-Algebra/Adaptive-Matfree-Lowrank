function [Q, Err_ind, Exact_err] = randQB_AF_Spec(A, b, epsilon, tol_type, compute_exact)
% randQB_AF_Spec  Adjoint-free, matrix-free adaptive range finder in the
%                 spectral norm (Section 4.3).
%
%   [Q, Err_ind, Exact_err] = randQB_AF_Spec(A, b, epsilon)
%   [...] = randQB_AF_Spec(A, b, epsilon, tol_type)
%   [...] = randQB_AF_Spec(A, b, epsilon, tol_type, compute_exact)
%
%   Builds Q (orthonormal columns) block by block using ONLY the forward
%   action A*X (no adjoint A'*X).  Mirrors randQB_AF_Fro: at each iteration the
%   adjoint-free residual sketch Y_i = (A - Q*Q'*A)*Omega_i is formed first, and
%   its scaled spectral norm ||Y_i||_2/sqrt(b) -- a proxy for ||E_{i-1}||_2 --
%   serves as the look-ahead stopping indicator (estimated by power iteration).
%   The block is orthonormalized (QR) only if the tolerance is not yet met.
%   Once the indicator drops below tau, the most recently retained block is
%   pruned via Prune_Rank_Spec, applied to the forward sketch
%   B_tilde = Q_{i-1}'*Z_i (a proxy for B_{i-1}), exactly as in randQB_AF_Fro.
%   The computed rank is size(Q,2).
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       tolerance (epsilon > 0)
%     tol_type      (optional, default 'relative')
%                     'relative' : threshold epsilon*||Y_1||_2/sqrt(b)
%                                  (the first indicator estimates ||A||_2).
%                     'absolute' : threshold epsilon.
%     compute_exact (optional, default false) if true, also return the exact
%                   residual spectral norm per iteration (forms A - Q*Q'*A;
%                   NOT matrix-free; diagnostics only).
%
%   Outputs
%     Q             orthonormal basis meeting the requested tolerance (w.h.p.)
%     Err_ind       adjoint-free indicator ||Y_i||_2/sqrt(b) per iteration
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

    % Preallocate the basis and per-iteration diagnostics.
    Q = zeros(m, maxit * b);
    ncol = 0;                           % columns of Q currently stored

    Err_ind = zeros(maxit, 1);
    if compute_exact
        Exact_err = zeros(maxit, 1);
    else
        Exact_err = [];
    end

    % Target tolerance: fixed for an absolute tolerance; for a relative
    % tolerance it is set once, at the first iteration, from the first block
    % via tau = epsilon*||Y_1||_2/sqrt(b) ~ epsilon*||A||_2.
    if ~use_relative
        tau = epsilon;
    end

    for i = 1:maxit

        Omega = randn(n, b);
        Z = A * Omega;

        % Adjoint-free residual sketch  Y = (A - Q*Q'*A)*Omega = Z - Q*(Q'*Z).
        if ncol == 0
            Y = Z;
        else
            Qc = Q(:, 1:ncol);
            Y  = Z - Qc * (Qc' * Z);
        end

        % Look-ahead indicator ||Y_i||_2/sqrt(b) ~ ||E_{i-1}||_2, estimated to a
        % couple of digits by power iteration (cost O(mb)).  Since Y_i = Q_i R_i
        % with Q_i orthonormal, this equals ||R_i||_2/sqrt(b).
        Err_ind(i) = spec_norm_power(Y) / sqrt(b);

        % For a relative tolerance, fix tau from the first block's norm.
        if use_relative && i == 1
            tau = epsilon * Err_ind(1);     % ||A||_2 estimated by ||Y_1||_2/sqrt(b)
        end

        if compute_exact
            if ncol == 0
                Exact_err(i) = norm(A);
            else
                Exact_err(i) = norm(A - Qc * (Qc' * A));
            end
        end

        if Err_ind(i) <= tau
            % Indicator below tolerance: discard this look-ahead block and
            % prune the previous block (i-1).
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end

            if ncol == 0
                % A already meets the tolerance: rank-0 approximation.
                Q = zeros(m, 0);
                return
            end

            keep   = ncol - b;                  % columns of blocks 1..i-2
            Q_prev = Q(:, keep+1:ncol);         % block i-1
            % Adjoint-free pruning: prune the forward sketch
            % B_tilde = Q_{i-1}'*Z_i = B_{i-1}*Omega_i (proxy for B_{i-1}), as in
            % randQB_AF_Fro.  The 1/sqrt(b) makes its singular values estimate
            % those of B_{i-1}.
            B_tilde = Q_prev' * Z;
            [U_r, ~] = Prune_Rank_Spec(B_tilde / sqrt(b), tau);

            Q = [Q(:, 1:keep), Q_prev * U_r];
            return
        end

        % Indicator still above tolerance: orthonormalize the block (QR after
        % the check, as in randQB_AF_Fro), re-orthogonalize, and append.
        [Qi, ~] = qr(Y, 'econ');
        if ncol > 0
            [Qi, ~] = qr(Qi - Qc * (Qc' * Qi), 'econ');
        end

        Q(:, ncol+1:ncol+b) = Qi;
        ncol = ncol + b;

        if ncol >= min_mn
            % Tolerance never met: return the full basis.
            Q = Q(:, 1:ncol);
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end
            return
        end
    end

    Q = Q(:, 1:ncol);
    Err_ind = Err_ind(1:i);
    if compute_exact
        Exact_err = Exact_err(1:i);
    end
end
