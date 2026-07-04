function [Q, Err_ind, Exact_err] = randQB_AF_Fro(A, b, epsilon, tol_type, compute_exact)
% randQB_AF_Fro  Adjoint-free, matrix-free adaptive randomized range finder
%                in the Frobenius norm (Algorithm 3.3).
%
%   [Q, Err_ind, Exact_err] = randQB_AF_Fro(A, b, epsilon)
%   [...] = randQB_AF_Fro(A, b, epsilon, tol_type)
%   [...] = randQB_AF_Fro(A, b, epsilon, tol_type, compute_exact)
%
%   Builds Q (orthonormal columns) block by block using ONLY the forward
%   action A*X (no adjoint A'*X) until the randomized indicator certifies the
%   target tolerance, then prunes the final block via Prune_Rank_Fro.  The
%   projection B = Q'*A is never formed; the pruning row norms are estimated
%   from the forward sketch B_tilde = Q'*Z.  The computed rank is size(Q,2).
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       tolerance (epsilon > 0)
%     tol_type      (optional, default 'relative')
%                     'relative' : stop when ||A-Q*Q'*A||_F <= epsilon*||A||_F
%                                  (||A||_F estimated on the fly).
%                     'absolute' : stop when ||A-Q*Q'*A||_F <= epsilon.
%     compute_exact (optional, default false) if true, also return the exact
%                   residual norm per iteration (forms A - Q*Q'*A; NOT
%                   matrix-free; diagnostics only).
%
%   Outputs
%     Q             orthonormal basis satisfying the requested tolerance
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

    % Preallocate the basis and per-iteration diagnostics.
    Q = zeros(m, maxit * b);
    ncol = 0;                       % number of columns currently stored in Q

    Err_ind = zeros(maxit, 1);
    if compute_exact
        Exact_err = zeros(maxit, 1);
    else
        Exact_err = [];
    end

    sum_Z = 0;                      % running sqrt(sum_k ||Z_k||_F^2) (relative only)

    for i = 1:maxit

        % Random sketch of A (forward matvec only).  Folding 1/sqrt(b) into
        % Omega makes ||Y||_F a direct estimate of the residual norm.
        Omega = (1/sqrt(b)) * randn(n, b);
        Z = A * Omega;

        % Target tolerance for this iteration.
        if use_relative
            sum_Z = hypot(sum_Z, norm(Z, 'fro'));
            tau   = epsilon * (sum_Z / sqrt(i));   % epsilon * (running ||A||_F estimate)
        else
            tau   = epsilon;
        end

        % Adjoint-free implicit residual sketch:
        %   Y = (A - Q*Q'*A)*Omega = Z - Q*(Q'*Z).
        if ncol == 0
            Y = Z;
        else
            Qc = Q(:, 1:ncol);
            Y  = Z - Qc * (Qc' * Z);
        end

        if compute_exact
            if ncol == 0
                Exact_err(i) = norm(A, 'fro');
            else
                Exact_err(i) = norm(A - Qc * (Qc' * A), 'fro');
            end
        end

        % Randomized stopping criterion.
        Err_ind(i) = norm(Y, 'fro');
        if Err_ind(i) < tau
            % Adjoint-free pruning: estimate the final block's row norms from
            % the forward sketch B_tilde = Q_{i-1}'*Z (proxy for Q_{i-1}'*A).
            keep    = ncol - b;             % columns kept from the earlier blocks
            Q_end   = Q(:, keep+1:ncol);    % last block Q_{i-1}
            B_tilde = Q_end' * Z;
            s = Prune_Rank_Fro(B_tilde, Err_ind(i), tau);

            % s indexes WITHIN the last block (a subset of 1:b), not the global
            % column numbering.  Keep the earlier blocks whole and append only
            % the s-selected columns of Q_end.
            Q = [Q(:, 1:keep), Q_end(:, s)];

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

        Q(:, ncol+1:ncol+b) = Qi;
        ncol = ncol + b;
    end

    % Tolerance never met: return the full basis.
    Q = Q(:, 1:ncol);
end
