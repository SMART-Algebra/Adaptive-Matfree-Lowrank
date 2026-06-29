function A = Generate_Test_Matrix(m, n, name)
% Generate_Test_Matrix  Random m-by-n test matrix A = U*S*V' with a
%                       prescribed singular-value decay profile.
%
%   A = Generate_Test_Matrix(m, n, name)
%
%   U and V have (random) orthonormal columns; the singular values are set by
%   the string 'name' (case-insensitive):
%
%     'exp'       fast/exponential decay   sigma_i = exp(-i/2)
%     'poly2'     polynomial decay         sigma_i = 1/i^2
%     'poly4'     polynomial decay         sigma_i = 1/i^4
%     's_shape'   S-shaped decay (piecewise: linear then polynomial in log10)
%
%   In the paper's notation: Matrix 1 = 'exp', Matrix 2 = 'poly4' (Frobenius)
%   or 'poly2' (spectral), and Matrix 3 = 's_shape'.

    K1 = randn(m, m);
    [U, ~] = qr(K1, 'econ');

    K2 = randn(n, n);
    [V, ~] = qr(K2, 'econ');

    p = min(m, n);
    d = zeros(1, p);

    % Generate the singular values.
    switch lower(name)
        case 'exp'
            d = exp(-(1:p)/2);                          % fast/exponential decay
        case 'poly2'
            d = (1:p).^(-2);                            % polynomial decay
        case 'poly4'
            d = (1:p).^(-4);                            % polynomial decay
        case 's_shape'
            % First part of the spectrum: linear decay (in log10) down to
            % 10^(-plin).
            plin = 2;            % linear piece ends at 10^(-plin)
            drop = 200;          % index where the linear piece ends
            piece_1 = (10^(-plin)).^(((1:drop)-1)/(drop-1));

            % Second part: polynomial decay (in log10).
            pp       = -10^(-0.8);
            sigma_N  = 5*eps;        % smallest singular value
            sigma_x2 = 10^(-plin);   % last value carried over from piece_1
            a  = (log10(sigma_N/sigma_x2)) / ((n-drop)^pp - 1);
            bb = log10(sigma_x2) - a;
            piece_2 = 10.^( a*(((drop+1):n)-drop).^pp + bb );

            d = [piece_1 piece_2];
        otherwise
            error('Generate_Test_Matrix:badName', ...
                ['Unknown test matrix ''%s''. Valid names: ', ...
                 'exp, poly2, poly4, s_shape.'], name);
    end

    S = spdiags(d', 0, m, n);
    A = U * S * V';
end
