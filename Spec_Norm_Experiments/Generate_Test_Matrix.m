function A = Generate_Test_Matrix(m, n, t)

    % Generate test matrices  in  "Insert paper name"

    %==================================================================%
    %% INPUT: Matrix dimensions m and n, value of t (1-5)
        % t = 1 (Fast Decay): singular values are  exp(-i/2)
        % t = 2 (Slow Decay 1): singular values are  1/sqrt(i)
        % t = 3 (Slow Decay 2): singular values are  1/(i^2)
        % t = 4 (Slow Decay 3): singular values are  1/(i^4) 
        % t = 5 (S-shaped Decay): singular values are  in the form of an S
    
    %% OUTPUT:  Matrix A where A  = U * Sigma * V'
    %==================================================================%


    K1 = randn(m, m);
    [U, ~] = qr(K1, 'econ');

    K2 = randn(n, n);
    [V, ~] = qr(K2, 'econ');


    p= min(m,n);
    d= zeros(1, p);

    % Generate singular values 
    switch t
        case 1 
            d= exp(-(1:p)/2);                             % fast decay
        case 2 
            d= (1:p).^(-0.5);                              % very slow decay
        case 3
             d= (1:p).^(-2);                                % moderate polynomialdecay
        case 4
             d= (1:p).^(-4);                                % fastest polynomial decay
        case 5
            d= 0.0001+1./(1+exp((1:p)-30));    %S-shaped Decay
        case 6
            % Form first part of spectrm (linear decay)
            p = 2;              % Choose where linear decay to stop (will be 10^(-p))
            drop = 200;
            piece_1 = (10^(-p)).^(((1:drop)-1)/(drop-1));

            % Form second part of spectrm (polynomial decay)
            pp = -10^(-0.8);
            xx = (drop+1):n;
            sigma_N =5*eps;  % Choose smallest singular value
            sigma_x2 = 10^(-p); % Get last singular value from linear decay 
            a = (log10(sigma_N/sigma_x2))/ (  (n-drop)^pp -1);
            bb = log10(sigma_x2) - a;
            piece_2 = 10.^( a*(((drop+1):n)-drop).^pp +bb  );

            % Form spectrum
            d = [piece_1 piece_2];
        otherwise
            disp('Not a valid t value!');
            exit;
    end

    S= spdiags(d', 0, m, n);
    A = U * S * V';