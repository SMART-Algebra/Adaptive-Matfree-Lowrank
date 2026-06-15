function [E, ERRORS, RANKS, REL_ERR] = comp_err_rank_S2(A, b, TOL, pow_its)

    norm_A = norm(A);
    ntol = length(TOL);

    E = TOL(:)';

    %==============================================================%
    % Storage
    %==============================================================%

    % 4 x ntol cell array:
    % row 1 = Exact error (from MF run)
    % row 2 = Adaptive Range Finder indicator
    % row 3 = MF indicator
    % row 4 = AF indicator
    ERRORS = cell(4, ntol);

    % ranks and relative errors
    RANKS = zeros(ntol, 4);   % [SVD, ARF, MF, AF]
    REL_ERR = zeros(ntol, 4);

    % preallocate vectors
    SVD_rank = zeros(ntol,1);
    ARF_rank = zeros(ntol,1);
    MF_rank  = zeros(ntol,1);
    AF_rank  = zeros(ntol,1);

    SVD_rel_err = zeros(ntol,1);
    ARF_rel_err = zeros(ntol,1);
    MF_rel_err  = zeros(ntol,1);
    AF_rel_err  = zeros(ntol,1);

    %==============================================================%
    % Main loop over tolerances
    %==============================================================%
    for i = 1:ntol

        tol = TOL(i);

        %-----------------------------
        % SVD baseline
        %-----------------------------
        [A_k, SVD_rank(i)] = truncated_svd_S_rel(A, tol);
        SVD_rel_err(i) = norm(A - A_k) / norm_A;

        %-----------------------------
        % MF method
        %-----------------------------
        %rng(2024)
        [Qmf, Bmf, MF_rank(i), ErrMF, ExactMF] = ...
            randQB_tol_S_rel_err_rank(A, b, tol, pow_its);

        MF_rel_err(i) = norm(A - Qmf*Bmf) / norm_A;

        %-----------------------------
        % AF method
        %-----------------------------
        %rng(2024)
        [Qaf, AF_rank(i), ErrAF, ~] = ...
            randQB_AF_Spec_New(A, b, tol);

        AF_rel_err(i) = norm(A - Qaf*(Qaf'*A)) / norm_A;

        %-----------------------------
        % Adaptive Range Finder
        %-----------------------------
        %rng(2024)
        [Qarf, ARF_rank(i), M] = ...
            adapt_range_finder(A, tol * norm(A,2), b);

        ARF_rel_err(i) = norm(A - Qarf*(Qarf'*A)) / norm_A;
            deviation = norm(Qarf'*Qarf - eye(size(Qarf,2)));
        % If ARF gives full rank, stop future runs
        if ARF_rank(i) >= size(A,2)
            fprintf('Adaptive Range Finder produced full rank matrix.\n')
            fprintf('Tolerance: %.1e\n', tol)

            % ARF_rank(i+1:end) = size(A,2);
            % ARF_rel_err(i+1:end) = NaN;
        end

        %======================================================%
        % STORE ERROR TRAJECTORIES
        %======================================================%

        % Match lengths for plotting
        arf_indicator = (10 * sqrt(2/pi)) * M(1:b:end)';

        if length(arf_indicator) > length(ExactMF)
            arf_indicator(length(ExactMF)+1:end) = [];
        end

        ERRORS{1,i} = sqrt(ExactMF);          % exact error
        ERRORS{2,i} = arf_indicator;    % adaptive range finder indicator
        ERRORS{3,i} = sqrt(ErrMF);            % MF indicator
        ERRORS{4,i} = ErrAF;            % AF indicator

        %======================================================%
        % STORE RANKS / REL ERRORS
        %======================================================%

        RANKS(i,:) = [ ...
            SVD_rank(i), ...
            ARF_rank(i), ...
            MF_rank(i), ...
            AF_rank(i)];

        REL_ERR(i,:) = [ ...
            SVD_rel_err(i), ...
            ARF_rel_err(i), ...
            MF_rel_err(i), ...
            AF_rel_err(i)];

    end

end