function [ERR_IND, RANKS, REL_ERR] = compare_Spec_methods(A, b, TOL)
% compare_Spec_methods  Compare the spectral-norm low-rank methods used in
%                       Section 5.3 across a list of relative tolerances.
%
%   [ERR_IND, RANKS, REL_ERR] = compare_Spec_methods(A, b, TOL)
%
%   For each tolerance in TOL, runs the truncated SVD baseline, randQB_HMT,
%   randQB_MF_Spec, and randQB_AF_Spec on A with block size b.  Let
%   ntol = numel(TOL).
%
%   Outputs
%     ERR_IND  ntol-by-4 cell of per-iteration error curves at each tolerance.
%              Column 1 is the exact residual spectral norm of the
%              randQB_MF_Spec run (reference, formed explicitly); columns 2-4
%              are the randomized error indicators of randQB_HMT /
%              randQB_MF_Spec / randQB_AF_Spec.
%     RANKS    ntol-by-4 computed ranks         [t-SVD, HMT, MF, AF]
%     REL_ERR  ntol-by-4 true relative errors   [t-SVD, HMT, MF, AF],
%              ||A - Q*B||_2 / ||A||_2 formed explicitly (exact, not estimated).

    norm_A = norm(A);                 % spectral norm
    ntol   = numel(TOL);

    ERR_IND = cell(ntol, 4);
    RANKS   = zeros(ntol, 4);
    REL_ERR = zeros(ntol, 4);

    for i = 1:ntol
        tol = TOL(i);

        % t-SVD baseline.
        [A_k, svd_rank] = ref_svd_Spec(A, tol, 'relative');
        svd_err = norm(A - A_k) / norm_A;

        % randQB_MF_Spec (compute_exact = true to also return the exact curve).
        [Qmf, Bmf, ErrMF, ExactMF] = randQB_MF_Spec(A, b, tol, 'relative', true);
        mf_err = norm(A - Qmf*Bmf) / norm_A;

        % randQB_AF_Spec (adjoint-free; form Q'*A only to score the error).
        [Qaf, ErrAF] = randQB_AF_Spec(A, b, tol, 'relative');
        af_err = norm(A - Qaf*(Qaf'*A)) / norm_A;

        % randQB_HMT baseline (uses an absolute tolerance tol*||A||_2).
        [Qhmt, hmt_rank, Mhmt] = ref_randQB_HMT(A, tol*norm_A, b);
        hmt_err = norm(A - Qhmt*(Qhmt'*A)) / norm_A;

        % Per-iteration curves: col 1 exact, cols 2-4 indicators.  The HMT
        % column-norm history is sampled on the same per-block grid and
        % truncated to the reference curve length for plotting.
        hmt_indicator = (10*sqrt(2/pi)) * Mhmt(1:b:end)';
        if numel(hmt_indicator) > numel(ExactMF)
            hmt_indicator = hmt_indicator(1:numel(ExactMF));
        end

        ERR_IND{i,1} = ExactMF(:)';
        ERR_IND{i,2} = hmt_indicator(:)';
        ERR_IND{i,3} = ErrMF(:)';
        ERR_IND{i,4} = ErrAF(:)';

        RANKS(i,:)   = [svd_rank, hmt_rank, size(Qmf,2), size(Qaf,2)];
        REL_ERR(i,:) = [svd_err,  hmt_err,  mf_err,      af_err];
    end

end
