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

    [~, n] = size(A);
    Q = zeros(size(A, 1), 0);
    j = 0;

    % Initial block of r random probes.
    Omega = randn(n, r);
    Y = A * Omega;

    M = max(vecnorm(Y(:, 1:r), 2, 1));            % running indicator history
    threshold = epsilon / (10 * sqrt(2/pi));

    % Stop when the maximal probe norm drops below the threshold.
    while M(end) > threshold && size(Q, 2) < n
        j = j + 1;

        % Re-orthogonalize probe j against Q, normalize, append.
        if isempty(Q)
            y_j = Y(:, j);
        else
            y_j = Y(:, j) - Q * (Q' * Y(:, j));
        end
        q_j = y_j / norm(y_j, 2);
        Q   = [Q, q_j];                            %#ok<AGROW>

        % Draw a fresh probe, orthogonalize against Q, append to Y.
        y_new = A * randn(n, 1);
        y_new = y_new - Q * (Q' * y_new);
        Y = [Y, y_new];                            %#ok<AGROW>

        % Re-orthogonalize the remaining live probes against the new q_j.
        for i = (j+1):(j+r-1)
            Y(:, i) = Y(:, i) - (q_j' * Y(:, i)) * q_j;
        end

        M = [M, max(vecnorm(Y(:, (j+1):(j+r)), 2, 1))];   %#ok<AGROW>
    end

    k = size(Q, 2);
end
