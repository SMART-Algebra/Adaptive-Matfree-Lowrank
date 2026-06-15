function [A_k, k] = truncated_svd_S_rel(A, tol)
    %==================================================================%
    % Inputs:  A   - input matrix
    %          tol - relative tolerance (0 < tol < 1)
    %
    % Outputs: A_k - truncated SVD approximation of A
    %          k   - numerical rank satisfying ||A - A_k||_2 <= tol*||A||_2
    %==================================================================%
    
    % Form SVD and obtain singular values
    [U, S, V] = svd(A, 'econ');   
    d = diag(S);         
    tol = tol * d(1) ;   % Obtain tolerance 

    % Determine k
    k = find(d > tol , 1, 'last');

    if isempty(k)
        k = 0;
        A_k = A;
        return
    end

    % Truncated SVD
    U_k = U(:, 1:k);
    S_k = S(1:k, 1:k);
    V_k = V(:, 1:k);

    A_k = U_k * S_k * V_k';

end




