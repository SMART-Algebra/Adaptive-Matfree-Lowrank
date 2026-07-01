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
    maxiter = ceil(min(m,n)/b);

    % Preallocate the factors; ncol tracks how many columns are in use.
    Q = zeros(m, maxiter*b);
    B = zeros(maxiter*b, n);
    ncol = 0;

    E = norm(A, 'fro')^2;  E0 = E;
    EE = E;                             % error-indicator history (small vector)
    if strcmpi(tol_type, 'relative')
        threshold = tol^2 * E;          % (tol * ||A||_F)^2
    else
        threshold = tol^2;              % absolute
    end
    flag = false;
    fail = false;
    for i = 1:maxiter
        Qc = Q(:, 1:ncol);
        Bc = B(1:ncol, :);

        Omg = randn(n, b);
        Y = A * Omg - (Qc * (Bc * Omg));
        [Qi, ~] = qr(Y, 0);

        for j = 1:P
            [Qi, ~] = qr(A'*Qi - Bc'*(Qc'*Qi), 0);  % can skip orthonormalization for small b.
            [Qi, ~] = qr(A*Qi - Qc*(Bc*Qi), 0);
        end

        [Qi, ~] = qr(Qi - Qc * (Qc' * Qi), 0);

        Bi = Qi' * A - Qi' * Qc * Bc;

        Q(:, ncol+1:ncol+b) = Qi;
        B(ncol+1:ncol+b, :) = Bi;
        ncol = ncol + b;

        temp = E - norm(Bi, 'fro')^2;
        if temp < 0
            warning('Error indicator became negative')
            fprintf('Tolerance: %.1e\n', tol)
        end

        if temp < threshold     % precise rank determination
            for j = 1:b
                E = E - norm(Bi(j,:))^2;
                if E < threshold
                    flag = true;
                    break;
                end
            end
        else
            E = temp;
            EE = [EE E];
        end
        if flag
            k = (i-1)*b + j;
            break;
        end
    end
    if ~flag
        k = i*b;
        fail = true;
        fprintf('E = %f. Fail to converge within maxiter!\n\n', sqrt(E/E0));
    end

    % Return only the filled portion of the preallocated factors.
    Q = Q(:, 1:ncol);
    B = B(1:ncol, :);
end
