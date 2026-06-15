function [Q, comp_rank, M] = adapt_range_finder(A, epsilon, r)
    % Initialize variables
    Q = [];                      % Orthonormal basis
    j = 0;                       % Iteration counter
    M = [];                      % Store maximal column norms for diagnostics

    % Draw initial random Gaussian matrix Omega
    Omega = randn(size(A, 2), r);
    Y = A * Omega;               % Project onto random subspace

    % Compute the initial maximal column norm
    col_norms = vecnorm(Y(:, 1:r), 2, 1);
    m = max(col_norms);
    M = [M m];

    % Main loop: Stop when the maximal norm is below the threshold
    while m > epsilon / (10 * sqrt(2 / pi)) && size(Q,2) < size(A,2)
        j = j + 1;

        % Orthogonalization step (numerically stable using QR)
        if isempty(Q)
            y_j = Y(:, j);
        else
            % Orthogonalize y_j with respect to Q
            y_j = Y(:, j) - Q * (Q' * Y(:, j));
        end
        % Normalize and add to Q
        q_j = y_j / norm(y_j, 2);
        Q = [Q, q_j];  % Expand the basis

        % Draw a new random Gaussian vector and compute the next column
        omega_jr = randn(size(A, 2), 1);
        y_new = A * omega_jr;
        
        % Orthogonalize the new column with respect to Q (numerically stable)
        y_new = y_new - Q * (Q' * y_new);
        Y = [Y, y_new];  % Append to Y
        
        % Orthogonalize remaining columns of Y
        for i = (j + 1):(j + r - 1)
            Y(:, i) = Y(:, i) - (q_j' * Y(:, i)) * q_j;
        end
        
        % Update and store maximal column norm
        col_norms = vecnorm(Y(:, (j+1):(j+r)), 2, 1);
        m = max(col_norms);
        M = [M m];
    end

    % Compute the rank of the final basis
    comp_rank = size(Q, 2);

    if comp_rank >= size(A,2)
    disp('Adaptive Range Finder produced a full-rank basis.')
    end
    
end
