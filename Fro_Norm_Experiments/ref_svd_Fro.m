function [A_k, k] = ref_svd_Fro(A, tol, tol_type)
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

    % Rank-k truncated SVD.  Scale the columns of U by the singular values
    % (a length-k vector) instead of forming the k-by-k diagonal matrix S_k,
    % avoiding the dense U_k*S_k product.
    s = diag(S);
    A_k = (U(:, 1:k) .* s(1:k).') * V(:, 1:k)';

end




