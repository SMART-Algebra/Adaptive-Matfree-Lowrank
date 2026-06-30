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
%   NOTE: Section 4.3 is still being finalized; this routine preserves the
%   numerics of the experimental adjoint-free spectral algorithm and only
%   tidies its structure / output interface.
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       tolerance (0 < epsilon < 1)
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

    Q = zeros(m, 0);
    Err_ind = zeros(maxit, 1);
    if compute_exact
        Exact_err = zeros(maxit, 1);
    else
        Exact_err = [];
    end

    R_prev = [];                        % R factor of the previous (retained) block
    Qorth_prev = [];                    % basis of the previous (retained) block

    for i = 1:maxit

        Omega = randn(n, b);
        Z = A * Omega;

        % Adjoint-free residual sketch  Y = (A - Q*Q'*A)*Omega = Z - Q*(Q'*Z).
        if isempty(Q)
            [Qbar, R] = qr(Z, 'econ');
            Qorth = Qbar;
        else
            Y = Z - Q * (Q' * Z);
            [Qbar, R] = qr(Y, 'econ');
            Qorth = Qbar - Q * (Q' * Qbar);
            [Qorth, ~] = qr(Qorth, 'econ');
        end

        % Adjoint-free indicator ||R_i||_2 / sqrt(b)  (proxy for ||B_i||_2).
        Err_ind(i) = norm(R) / sqrt(b);

        if use_relative
            tau = epsilon * Err_ind(1);
        else
            tau = epsilon;
        end

        if compute_exact
            if isempty(Q)
                Exact_err(i) = norm(A);
            else
                Exact_err(i) = norm(A - Q * (Q' * A));
            end
        end

        if Err_ind(i) < tau
            % Indicator below tolerance: prune the most recent block, finalize.
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end

            if i == 1
                Rk = R;  Qk = Qorth;            % nothing retained yet
            else
                Rk = R_prev;  Qk = Qorth_prev;  % the last appended block
            end

            [U, S, ~] = svd(Rk);
            svals = diag(S) / sqrt(b);
            r = find(svals >= tau);

            if i == 1
                Q = Qk * U(:, r);
            else
                Q = [Q(:, 1:(i-2)*b), Qk * U(:, r)];
            end
            return
        end

        % Indicator still above tolerance: append the block.
        Q = [Q, Qorth];                 %#ok<AGROW>
        R_prev = R;  Qorth_prev = Qorth;

        if size(Q, 2) >= min_mn
            % Tolerance never met: return the full basis.
            Q = Q(:, 1:min_mn);
            Err_ind = Err_ind(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end
            return
        end
    end

    Err_ind = Err_ind(1:i);
    if compute_exact
        Exact_err = Exact_err(1:i);
    end
end
