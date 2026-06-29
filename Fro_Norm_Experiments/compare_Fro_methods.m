function [ERR_IND, RANKS, REL_ERR] = compare_Fro_methods(A, b, TOL)
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

    norm_A = norm(A,'fro');
    ntol = length(TOL);

    %==============================================================%
    % Storage
    %==============================================================%

    % ntol x 4 cell of per-iteration error curves:
    % col 1 = exact residual (from the MF run), col 2 = EI indicator,
    % col 3 = MF indicator, col 4 = AF indicator.
    ERR_IND = cell(ntol, 4);

    % ranks and relative errors
    RANKS = zeros(ntol, 4);   % [SVD, EI, MF, AF]
    REL_ERR = zeros(ntol, 4);

    % preallocate vectors
    SVD_rank = zeros(ntol,1);
    EI_rank  = zeros(ntol,1);
    MF_rank  = zeros(ntol,1);
    AF_rank  = zeros(ntol,1);

    SVD_rel_err = zeros(ntol,1);
    EI_rel_err  = zeros(ntol,1);
    MF_rel_err  = zeros(ntol,1);
    AF_rel_err  = zeros(ntol,1);

    EI_fail = false(ntol,1);

    %==============================================================%
    % Main loop over tolerances
    %==============================================================%
    for i = 1:ntol

        tol = TOL(i);

        %-----------------------------
        % SVD baseline
        %-----------------------------
        [A_k, SVD_rank(i)] = truncated_svd_F_rel(A, tol);
        SVD_rel_err(i) = norm(A - A_k,'fro') / norm_A;

        %-----------------------------
        % MF method (randQB_MF_Fro)
        %-----------------------------
        %rng(2024)
        [Qmf, Bmf, ErrMF, ExactMF] = randQB_MF_Fro(A, b, tol, 'relative', true);
        MF_rank(i) = size(Qmf, 2);
        MF_rel_err(i) = norm(A - Qmf*Bmf,'fro') / norm_A;

        %-----------------------------
        % AF method (randQB_AF_Fro)
        %-----------------------------
        %rng(2024)
        [Qaf, ErrAF] = randQB_AF_Fro(A, b, tol, 'relative');
        AF_rank(i) = size(Qaf, 2);
        Baf = Qaf' * A;
        AF_rel_err(i) = norm(A - Qaf*Baf,'fro') / norm_A;

        %-----------------------------
        % EI method
        %-----------------------------
        %rng(2024)
        [Qei, Bei, EI_rank(i), EI_err_ind, fail] = randQB_EI_err_rank(A, tol, b, 0);
        Qei = Qei(:, 1:EI_rank(i));
        Bei = Bei(1:EI_rank(i), :);
        EI_rel_err(i) = norm(A - Qei*Bei,'fro') / norm_A;
        EI_fail(i) = fail;

        %======================================================%
        % STORE ERROR TRAJECTORIES
        %======================================================%

        ERR_IND{i,1} = ExactMF;              % exact residual (reference)
        ERR_IND{i,2} = sqrt(EI_err_ind(:)'); % randQB_EI indicator
        ERR_IND{i,3} = ErrMF;                % randQB_MF_Fro indicator
        ERR_IND{i,4} = ErrAF;                % randQB_AF_Fro indicator

        %======================================================%
        % STORE RANKS / REL ERRORS
        %======================================================%
        RANKS(i,:) = [SVD_rank(i), EI_rank(i), MF_rank(i), AF_rank(i)];

        REL_ERR(i,:) = [SVD_rel_err(i), EI_rel_err(i), MF_rel_err(i), AF_rel_err(i)];

    end

end