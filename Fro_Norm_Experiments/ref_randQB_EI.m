function [Q, B, k, EE, fail] = ref_randQB_EI(A, tol, b, P, tol_type)
% The fixed-precision randQB_EI algorithm (reference method).
% It produces a QB factorization of A whose approximation error fulfills
%     ||A-QB||_F <= tol*||A||_F   (tol_type = 'relative', default), or
%     ||A-QB||_F <= tol           (tol_type = 'absolute').
% b is block size, P is power parameter.
% Output k is the rank.
%
% Adapted from Wenjian Yu's reference implementation randQB_EI_auto.m:
%   https://github.com/WenjianYu/randQB_auto/blob/master/randQB_EI_auto.m
% We made only minor changes (extra outputs for the error indicator and a
% convergence-failure flag, plus an absolute-tolerance option) so it can be
% compared against our methods.
    if nargin < 5 || isempty(tol_type)
        tol_type = 'relative';
    end
    [m, n]  = size(A);
    EE  = [ ];
    Q = zeros(m, 0);
    B = zeros(0, n);
    E= norm(A, 'fro')^2; E0= E;
    EE = [EE E];
    if strcmpi(tol_type, 'relative')
        threshold = tol^2 * E;      % (tol * ||A||_F)^2
    else
        threshold = tol^2;          % absolute
    end
    maxiter= ceil(min(m,n)/b);
    flag= false;
    fail = false;
    for i=1:maxiter
        Omg = randn(n, b);
        Y = A * Omg - (Q * (B * Omg));
        [Qi, ~] = qr(Y, 0);
        
        for j = 1:P
            [Qi, ~] = qr(A'*Qi - B'*(Q'*Qi), 0);  % can skip orthonormalization for small b. 
            [Qi, ~] = qr(A*Qi - Q*(B*Qi), 0);
        end
        
        [Qi, ~] = qr(Qi - Q * (Q' * Qi), 0);
        
        Bi = Qi' * A - Qi' * Q * B;
        
        Q = [Q, Qi];
        B = [B; Bi];
        
        temp = E- norm(Bi, 'fro')^2;
        if temp<0
            warning('Error indicator became negative')
            fprintf('Tolerance: %.1e\n', tol)
        end

        if temp< threshold,   % precise rank determination
            for j=1:b,
                E= E-norm(Bi(j,:))^2;
                if E< threshold,
                    flag= true;
                    break;
                end
            end
        else
            E= temp;
            EE = [EE E];
        end
        if flag,
            k= (i-1)*b + j;
            break;
        end
        
    end
    if ~flag,
        k= i*b;
        fail = true;
        fprintf('E = %f. Fail to converge within maxiter!\n\n', sqrt(E/E0));
    end
end