function [ERR_IND, RANKS, REL_ERR, EI_NEG] = compare_Fro_methods(A, b, TOL)
% compare_Fro_methods  Compare the Frobenius-norm low-rank methods used in
%                      Section 5.2 across a list of relative tolerances.
%
%   [ERR_IND, RANKS, REL_ERR] = compare_Fro_methods(A, b, TOL)
%
%   For each tolerance in TOL, runs the truncated SVD baseline, randQB_EI,
%   randQB_MF_Fro, and randQB_AF_Fro on A with block size b.  Let
%   ntol = numel(TOL).
%
%   Outputs
%     ERR_IND  ntol-by-4 cell of per-iteration error curves at each tolerance.
%              Column 1 is the exact residual norm of the randQB_MF_Fro run
%              (reference, formed explicitly); columns 2-4 are the randomized
%              error indicators computed internally by randQB_EI /
%              randQB_MF_Fro / randQB_AF_Fro (i.e. estimates of the error).
%     RANKS    ntol-by-4 computed ranks         [t-SVD, EI, MF, AF]
%     REL_ERR  ntol-by-4 true relative errors   [t-SVD, EI, MF, AF],
%              ||A - Q*B||_F / ||A||_F formed explicitly (exact, not an estimate).
%     EI_NEG   ntol-by-1 logical; true at tolerances where randQB_EI's error
%              indicator went negative (Gramian cancellation).

    norm_A = norm(A,'fro');
    ntol   = numel(TOL);

    % ntol x 4 cell of per-iteration error curves:
    % col 1 = exact residual (from the MF run), col 2 = EI indicator,
    % col 3 = MF indicator, col 4 = AF indicator.
    ERR_IND = cell(ntol, 4);

    % ntol x 4, columns [t-SVD, randQB_EI, randQB_MF_Fro, randQB_AF_Fro].
    RANKS   = zeros(ntol, 4);
    REL_ERR = zeros(ntol, 4);
    EI_NEG  = false(ntol, 1);           % randQB_EI indicator went negative?

    for i = 1:ntol
        tol = TOL(i);

        % t-SVD baseline.
        [A_k, svd_rank] = ref_svd_Fro(A, tol, 'relative');
        svd_err = norm(A - A_k, 'fro') / norm_A;

        % randQB_MF_Fro (compute_exact = true to also return the exact curve).
        [Qmf, Bmf, ErrMF, ExactMF] = randQB_MF_Fro(A, b, tol, 'relative', true);
        mf_err = norm(A - Qmf*Bmf, 'fro') / norm_A;

        % randQB_AF_Fro (adjoint-free; form B = Q'*A only to score the error).
        [Qaf, ErrAF] = randQB_AF_Fro(A, b, tol, 'relative');
        af_err = norm(A - Qaf*(Qaf'*A), 'fro') / norm_A;

        % randQB_EI reference.
        [Qei, Bei, ei_rank, EI_err_ind, ~, ei_neg] = ref_randQB_EI(A, tol, b, 0, 'relative');
        ei_err = norm(A - Qei(:,1:ei_rank)*Bei(1:ei_rank,:), 'fro') / norm_A;
        EI_NEG(i) = ei_neg;

        % Per-iteration error curves: col 1 exact, cols 2-4 the indicators.
        ERR_IND{i,1} = ExactMF;
        ERR_IND{i,2} = sqrt(EI_err_ind(:)');
        ERR_IND{i,3} = ErrMF;
        ERR_IND{i,4} = ErrAF;

        RANKS(i,:)   = [svd_rank, ei_rank, size(Qmf,2), size(Qaf,2)];
        REL_ERR(i,:) = [svd_err,  ei_err,  mf_err,      af_err];
    end

end