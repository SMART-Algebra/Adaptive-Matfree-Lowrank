function [A_k, k] = ref_svd_Spec(A, tol, tol_type)
    %==================================================================%
    % Reference truncated-SVD baseline in the spectral norm.
    %
    % Inputs:  A        - input matrix
    %          tol      - tolerance (0 < tol < 1)
    %          tol_type - (optional, default 'relative')
    %                       'relative' : ||A - A_k||_2 <= tol*||A||_2
    %                       'absolute' : ||A - A_k||_2 <= tol
    %
    % Outputs: A_k - truncated SVD approximation of A
    %          k   - numerical rank meeting the requested tolerance
    %==================================================================%

    if nargin < 3 || isempty(tol_type)
        tol_type = 'relative';
    end

    [U, S, V] = svd(A, 'econ');
    s = diag(S);

    % The spectral-norm truncation error is sigma_{k+1}, so keep every
    % singular value strictly above the threshold.
    if strcmpi(tol_type, 'relative')
        thr = tol * s(1);             % tol * ||A||_2
    else
        thr = tol;
    end

    k = find(s > thr, 1, 'last');
    if isempty(k)
        k = 0;                        % even rank 0 meets the tolerance
        A_k = zeros(size(A));
        return
    end

    % Rank-k truncated SVD.  Scale the columns of U by the singular values
    % (a length-k vector) instead of forming the dense k-by-k diagonal S_k.
    A_k = (U(:, 1:k) .* s(1:k).') * V(:, 1:k)';

end
