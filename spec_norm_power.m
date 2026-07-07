function s = spec_norm_power(M)
% spec_norm_power  Estimate ||M||_2 (largest singular value) by power
% iteration on M'*M.  Only a couple of correct digits are needed for the
% stopping indicator, so a loose relative tolerance and a small iteration
% cap suffice; each iteration costs two matvecs with M.  Works for either
% orientation of M (tall or wide), since it only forms M*v and M'*w.
%
%   Shared by randQB_MF_Spec (M = B_i, short-and-wide) and randQB_AF_Spec
%   (M = Y_i, tall), so both spectral-norm methods use the same estimator.

    reltol = 1e-2;      % ~2 correct digits in the singular value
    maxit  = 100;       % safety cap (typically only a handful of iterations)

    v  = randn(size(M, 2), 1);
    nv = norm(v);
    if nv == 0
        s = 0; return
    end
    v = v / nv;

    lambda = 0;
    for k = 1:maxit
        w = M' * (M * v);               % (M'*M) * v
        lambda_new = norm(w);
        if lambda_new == 0
            s = 0; return               % M is (numerically) zero
        end
        v = w / lambda_new;
        if abs(lambda_new - lambda) <= reltol * lambda_new
            lambda = lambda_new;
            break
        end
        lambda = lambda_new;
    end

    s = sqrt(lambda);                   % lambda -> sigma_max(M)^2
end
