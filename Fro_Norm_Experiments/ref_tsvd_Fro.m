function [A_k, k] = ref_tsvd_Fro(A, tol, tol_type)
    %==================================================================%
    % Reference truncated-SVD baseline in the Frobenius norm.
    %
    % Inputs:  A        - input matrix
    %          tol      - tolerance (0 < tol < 1)
    %          tol_type - (optional, default 'relative')
    %                       'relative' : ||A - A_k||_F <= tol*||A||_F
    %                       'absolute' : ||A - A_k||_F <= tol
    %
    % Outputs: A_k - truncated SVD approximation of A
    %          k   - numerical rank meeting the requested tolerance
    %==================================================================%

    if nargin < 3 || isempty(tol_type)
        tol_type = 'relative';
    end

    % Form SVD and obtain squared singular values
    [U, S, V] = svd(A, 'econ');
    d = diag(S).^2;                   % squared singular values

    % Squared error threshold on the discarded tail.
    if strcmpi(tol_type, 'relative')
        thr = (tol * norm(A,'fro'))^2;
    else
        thr = tol^2;
    end

    % Determine k
    cumulative_sum_sing = cumsum(flip(d));
    k_remove = find(cumulative_sum_sing  <= thr, 1, 'last');  % last index where the tail energy <= thr
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




