function [Q, Err_ind, Exact_err] = randQB_AF_Spec(A, b, epsilon, tol_type, compute_exact)
% randQB_AF_Spec  Adjoint-free, matrix-free adaptive range finder in the
%                 spectral norm (Section 4.3).
%
%   [Q, Err_ind, Exact_err] = randQB_AF_Spec(A, b, epsilon)
%   [...] = randQB_AF_Spec(A, b, epsilon, tol_type)
%   [...] = randQB_AF_Spec(A, b, epsilon, tol_type, compute_exact)
%
%   Builds Q (orthonormal columns) block by block using ONLY the forward
%   action A*X (no adjoint A'*X).  The triangular factor R_i of the residual
%   sketch supplies the adjoint-free indicator ||R_i||_2/sqrt(b), a proxy for
%   ||B_i||_2.  Once the indicator drops below tau, the most recently retained
%   block is truncated via the SVD of its R factor.  The computed rank is
%   size(Q,2).
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       tolerance (epsilon > 0)
%     tol_type      (optional, default 'relative')
%                     'relative' : threshold epsilon*||R_1||_2/sqrt(b)
%                                  (the first indicator estimates ||A||_2).
%                     'absolute' : threshold epsilon.
%     compute_exact (optional, default false) if true, also return the exact
%                   residual spectral norm per iteration (forms A - Q*Q'*A;
%                   NOT matrix-free; diagnostics only).
%
%   Outputs
%     Q             orthonormal basis meeting the requested tolerance (w.h.p.)
%     Err_ind       adjoint-free indicator ||R_i||_2/sqrt(b) per iteration
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

    R_prev = [];                        % R factor of the previous (retained) block

    for i = 1:maxit

        Omega = randn(n, b);
        Z = A * Omega;

        % Adjoint-free residual sketch  Y = (A - Q*Q'*A)*Omega = Z - Q*(Q'*Z).
        if ncol == 0
            [Qbar, R] = qr(Z, 'econ');
            Qi = Qbar;
        else
            Qc = Q(:, 1:ncol);
            Y  = Z - Qc * (Qc' * Z);
            [Qbar, R] = qr(Y, 'econ');
            Qi = Qbar - Qc * (Qc' * Qbar);
            [Qi, ~] = qr(Qi, 'econ');
        end

        % Adjoint-free indicator ||R_i||_2 / sqrt(b)  (proxy for ||B_i||_2).
        Err_ind(i) = norm(R) / sqrt(b);

        if use_relative
            tau = epsilon * Err_ind(1);
        else
            tau = epsilon;
        end

        if compute_exact
            if ncol == 0
                Exact_err(i) = norm(A);
            else
                Exact_err(i) = norm(A - Qc * (Qc' * A));
            end
        end

        if Err_ind(i) < tau
            % Indicator below tolerance: prune the most recent block, finalize.
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end

            if ncol == 0
                Rk = R;  Qk = Qi;  keep = 0;        % nothing retained yet (block 1)
            else
                Rk = R_prev;  keep = ncol - b;      % blocks 1..i-2 are kept whole
                Qk = Q(:, keep+1:ncol);             % basis of the last retained block
            end

            [U, S, ~] = svd(Rk);
            svals = diag(S) / sqrt(b);
            s = find(svals >= tau);

            Q = [Q(:, 1:keep), Qk * U(:, s)];
            return
        end

        % Indicator still above tolerance: append the block.
        Q(:, ncol+1:ncol+b) = Qi;
        ncol = ncol + b;
        R_prev = R;

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
