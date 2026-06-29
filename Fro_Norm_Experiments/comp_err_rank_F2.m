function [E, ERRORS, RANKS, REL_ERR] = comp_err_rank_F2(A, b, TOL)

    norm_A = norm(A,'fro');
    ntol = length(TOL);

    E = TOL(:)';

    %==============================================================%
    % Storage
    %==============================================================%

    % 4 x ntol cell array:
    % row 1 = Exact error (from MF run)
    % row 2 = EI indicator
    % row 3 = MF indicator
    % row 4 = AF indicator
    ERRORS = cell(4, ntol);

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

        ERRORS{1,i} = ExactMF;        % exact residual energy curve
        ERRORS{2,i} = sqrt(EI_err_ind(:)'); % EI indicator
        ERRORS{3,i} = ErrMF;          % MF indicator
        ERRORS{4,i} = ErrAF;          % AF indicator

        %======================================================%
        % STORE RANKS / REL ERRORS
        %======================================================%
        RANKS(i,:) = [SVD_rank(i), EI_rank(i), MF_rank(i), AF_rank(i)];

        REL_ERR(i,:) = [SVD_rel_err(i), EI_rel_err(i), MF_rel_err(i), AF_rel_err(i)];

    end

end