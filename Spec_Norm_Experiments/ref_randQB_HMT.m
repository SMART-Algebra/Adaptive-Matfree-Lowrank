function [Q, k, M] = ref_randQB_HMT(A, epsilon, r)
% ref_randQB_HMT  Adaptive randomized range finder of Halko, Martinsson &
%                 Tropp (the "randQB_HMT" spectral-norm baseline; Algorithm
%                 4.2 in [Halko, Martinsson, Tropp, SIAM Review, 2011]).
%
%   [Q, k, M] = ref_randQB_HMT(A, epsilon, r)
%
%   Grows an orthonormal basis Q one column at a time, using r random
%   Gaussian probe vectors to estimate the residual spectral norm, and stops
%   once that estimate falls below epsilon.  The probabilistic stopping rule
%   is  max_j ||y_j||_2 <= epsilon / (10*sqrt(2/pi)).
%
%   Inputs
%     A        m-by-n matrix
%     epsilon  ABSOLUTE spectral-norm tolerance (the caller scales by ||A||_2
%              for a relative criterion)
%     r        number of random probe vectors (block / oversampling parameter)
%
%   Outputs
%     Q   orthonormal basis with ||A - Q*Q'*A||_2 <~ epsilon (w.h.p.)
%     k   computed rank, size(Q,2)
%     M   history of the maximal probe-column norm (error-indicator trace)

    [m, n] = size(A);
    cap = min(m, n);                % at most min(m,n) orthonormal columns

    % Preallocate the basis, the probe buffer, and the indicator history so
    % none of them are grown inside the loop.
    Q  = zeros(m, cap);
    kq = 0;                         % basis columns accepted so far
    Y  = zeros(m, cap + r);         % probe buffer: initial r probes + one per step
    ny = r;                         % probe columns currently in use
    M  = zeros(1, cap + 1);         % indicator history
    nm = 1;

    Y(:, 1:r) = A * randn(n, r);
    M(1) = max(vecnorm(Y(:, 1:r), 2, 1));
    threshold = epsilon / (10 * sqrt(2/pi));

    j = 0;
    % Stop when the maximal probe norm drops below the threshold.
    while M(nm) > threshold && kq < cap
        j = j + 1;

        % Re-orthogonalize probe j against Q, normalize, accept it.
        if kq == 0
            y_j = Y(:, j);
        else
            Qc  = Q(:, 1:kq);
            y_j = Y(:, j) - Qc * (Qc' * Y(:, j));
        end
        q_j = y_j / norm(y_j, 2);
        kq = kq + 1;
        Q(:, kq) = q_j;

        % Draw a fresh probe, orthogonalize against Q, append to the buffer.
        y_new = A * randn(n, 1);
        Qc = Q(:, 1:kq);
        y_new = y_new - Qc * (Qc' * y_new);
        ny = ny + 1;
        Y(:, ny) = y_new;

        % Re-orthogonalize the remaining live probes against the new q_j.
        for t = (j+1):(j+r-1)
            Y(:, t) = Y(:, t) - (q_j' * Y(:, t)) * q_j;
        end

        nm = nm + 1;
        M(nm) = max(vecnorm(Y(:, (j+1):(j+r)), 2, 1));
    end

    Q = Q(:, 1:kq);
    M = M(1:nm);
    k = kq;
end
