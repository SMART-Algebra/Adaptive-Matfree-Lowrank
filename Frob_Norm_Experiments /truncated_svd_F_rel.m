function [A_k, k] = truncated_svd_F_rel(A, tol)
    %==================================================================%
    % Inputs:  A   - input matrix
    %          tol - relative tolerance (0 < tol < 1)
    %
    % Outputs: A_k - truncated SVD approximation of A
    %          k   - numerical rank satisfying ||A - A_k||_F <= tol*||A||_F
    %==================================================================%
    
    % Form SVD and obtain squared singular values
    [U, S, V] = svd(A, 'econ');   
    d = diag(S).^2;                   % squared singular values
    tol= (tol * norm(A,'fro'))^2; % squared tolerance

    % Determine k
    cumulative_sum_sing = cumsum(flip(d));
    k_remove = find(cumulative_sum_sing  <= tol, 1, 'last');  % last index where cumulative sum <= tol
    if isempty(k_remove)
        k = length(d);  % keep all singular values
    else
        k = length(d) - k_remove;
    end

    % Truncated SVD
    U_k = U(:, 1:k);
    S_k = S(1:k, 1:k);
    V_k = V(:, 1:k);

    A_k = U_k * S_k * V_k';

end




