function [Q, B, k, EE, fail, neg] = ref_randQB_EI(A, tol, b, P, tol_type)
% ref_randQB_EI  Fixed-precision randomized QB with an energy-based error
%                indicator (the "randQB_EI" baseline).
%
%   [Q, B, k, EE, fail, neg] = ref_randQB_EI(A, tol, b, P, tol_type)
%
%   Independent MATLAB implementation of the randQB_EI algorithm of
%       W. Yu, Y. Gu, and Y. Li, "Efficient randomized algorithms for the
%       fixed-precision low-rank matrix approximation," SIAM J. Matrix Anal.
%       Appl. 39(3):1339-1359, 2018.
%   Written from the algorithm in that paper and provided here only as a
%   comparison baseline for the proposed methods.
%
%   The method grows a QB factorization in blocks of b columns and tracks the
%   squared residual through the energy identity
%       ||A - QB||_F^2 = ||A||_F^2 - ||B||_F^2,
%   updated incrementally as  E <- E - ||B_i||_F^2.  This "top-down"
%   subtraction is prone to catastrophic cancellation once the residual is
%   small relative to ||A||_F -- the failure mode analyzed in the paper and the
%   motivation for the sketched indicators used by the proposed methods.
%
%   Inputs
%     A         m-by-n matrix
%     tol       tolerance (see tol_type)
%     b         block size
%     P         power-iteration parameter (0 for none)
%     tol_type  (optional, default 'relative')
%                 'relative' : stop when ||A-QB||_F <= tol*||A||_F
%                 'absolute' : stop when ||A-QB||_F <= tol
%
%   Outputs
%     Q, B   QB factors (Q has orthonormal columns).  The final block is kept
%            whole, so size(Q,2) may exceed the reported rank k.
%     k      computed rank meeting the tolerance (k <= size(Q,2)); obtained by
%            peeling rows of the last block until the energy crosses the target.
%     EE     history of the squared-residual indicator E, one entry per block
%            up to termination (EE(1) = ||A||_F^2).  E may go negative once
%            cancellation dominates (see neg).
%     fail   true if the tolerance was never met within ceil(min(m,n)/b) blocks
%            (the indicator plateaued and the method ran to full rank).
%     neg    true if the indicator E became negative at some iteration
%            (Gramian cancellation).

    if nargin < 4 || isempty(P),        P = 0;                end
    if nargin < 5 || isempty(tol_type), tol_type = 'relative'; end

    [m, n]  = size(A);
    maxiter = ceil(min(m, n) / b);

    normA2 = norm(A, 'fro')^2;
    if strcmpi(tol_type, 'relative')
        thresh = tol^2 * normA2;            % (tol*||A||_F)^2
    else
        thresh = tol^2;                     % absolute
    end

    % Preallocate the factors; ncol counts the columns currently in use.
    Q = zeros(m, maxiter * b);
    B = zeros(maxiter * b, n);
    ncol = 0;

    E   = normA2;                           % running ||A - QB||_F^2
    EE  = E;                                % indicator history (starts at ||A||_F^2)
    neg       = false;
    converged = false;
    k = maxiter * b;

    for i = 1:maxiter
        Qc = Q(:, 1:ncol);
        Bc = B(1:ncol, :);

        % Sketch the current residual (A - Qc*Bc) with a Gaussian block.
        Omega = randn(n, b);
        Yi = A * Omega - Qc * (Bc * Omega);
        [Qi, ~] = qr(Yi, 0);

        % Optional power iteration on the residual operator (A - Qc*Bc).
        for j = 1:P
            [Qi, ~] = qr(A' * Qi - Bc' * (Qc' * Qi), 0);
            [Qi, ~] = qr(A * Qi - Qc * (Bc * Qi), 0);
        end

        % Re-orthogonalize the new block against the existing basis.
        [Qi, ~] = qr(Qi - Qc * (Qc' * Qi), 0);
        Bi = Qi' * A - (Qi' * Qc) * Bc;     % = Qi'*(A - Qc*Bc), numerically stabilized

        Q(:, ncol+1:ncol+b) = Qi;
        B(ncol+1:ncol+b, :) = Bi;
        ncol = ncol + b;

        % Energy-subtraction error indicator for the whole new block.
        Eblock = E - norm(Bi, 'fro')^2;
        if Eblock < 0
            neg = true;
        end

        if Eblock <= thresh
            % Tolerance crossed inside this block: peel its rows one at a time
            % to recover the smallest rank that still meets the target.
            for j = 1:b
                E = E - norm(Bi(j, :))^2;
                if E <= thresh
                    break
                end
            end
            k = (i - 1) * b + j;
            converged = true;
            break
        else
            E  = Eblock;
            EE = [EE, E];                   % record the indicator trajectory
        end
    end

    fail = ~converged;                      % never met tol -> ran to full rank

    % Return only the filled portion of the preallocated factors.
    Q = Q(:, 1:ncol);
    B = B(1:ncol, :);
end
