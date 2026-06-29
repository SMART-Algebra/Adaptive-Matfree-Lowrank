% Drive_Fro_exp  Section 5.2 results for the fast-decaying spectrum.
%
%   Test matrix: Matrix 1, exponential decay  sigma_i = exp(-i/2).
%   Produces the two Figure-2 panels for this matrix:
%     (a) per-iteration error indicators (at the smallest tolerance), and
%     (b) computed rank vs. prescribed tolerance,
%   comparing t-SVD, randQB_EI, randQB_MF_Fro, and randQB_AF_Fro.
%
%   The heavy comparison is delegated to compare_Fro_methods; the plotting is
%   kept here so this matrix can be tuned independently of the other drivers.

clear; close all;

%% Problem setup
N   = 2^10;                     % matrix size (N x N)
b   = 16;                      % block size
p   = 14;                       % tolerances 10^-1 ... 10^-p
TOL = 10.^(-1:-1:-p);

A = Generate_Test_Matrix(N, N, 'exp');

%% Run the comparison
[ERR_IND, RANKS, REL_ERR] = compare_Fro_methods(A, b, TOL);

%% Line/marker styles (bright + distinct for the new methods; black reference)
MF  = {'Color',[0.00 0.45 0.90],'LineStyle','-.', 'Marker','o','MarkerSize',15, ...
       'LineWidth',3.5,'MarkerFaceColor',[0.00 0.45 0.90]};   % randQB_MF_Fro
AF  = {'Color',[0.95 0.45 0.05],'LineStyle','--','Marker','s','MarkerSize',9,  ...
       'LineWidth',3.5,'MarkerFaceColor',[0.95 0.45 0.05]};   % randQB_AF_Fro
EI  = {'Color',[0.35 0.35 0.35],'LineStyle','-','Marker','*','MarkerSize',20, ...
       'LineWidth',3.0};                                      % randQB_EI
REF = {'Color',[0 0 0],         'LineStyle','-', 'Marker','none','LineWidth',3.0}; % t-SVD / exact

% (a) Error indicators at the smallest tolerance
figure('Name','exp: error indicator');
hRef = semilogy(1:numel(ERR_IND{end,1}), ERR_IND{end,1}, REF{:}); hold on
hEI  = semilogy(1:numel(ERR_IND{end,2}), ERR_IND{end,2}, EI{:});
hMF  = semilogy(1:numel(ERR_IND{end,3}), ERR_IND{end,3}, MF{:});
hAF  = semilogy(1:numel(ERR_IND{end,4}), ERR_IND{end,4}, AF{:}); hold off
legend([hMF hAF hEI hRef], ...
    {'\texttt{randQB\_MF\_Fro}','\texttt{randQB\_AF\_Fro}','\texttt{randQB\_EI}','Exact error'}, ...
    'Interpreter','latex','Location','best');
xlabel('Iteration $i$','Interpreter','latex');
ylabel('Error estimator','Interpreter','latex');
grid on; set(gca,'FontSize',22);

% (b) Computed rank vs. prescribed tolerance
mEI = RANKS(:,2) ~= N;          % drop tolerances where randQB_EI returns full rank
figure('Name','exp: rank vs tolerance');
hRef = semilogx(TOL,      RANKS(:,1),   REF{:}); hold on   % t-SVD
hEI  = semilogx(TOL(mEI), RANKS(mEI,2), EI{:});            % randQB_EI
hMF  = semilogx(TOL,      RANKS(:,3),   MF{:});            % randQB_MF_Fro
hAF  = semilogx(TOL,      RANKS(:,4),   AF{:}); hold off   % randQB_AF_Fro
legend([hMF hAF hEI hRef], ...
    {'\texttt{randQB\_MF\_Fro}','\texttt{randQB\_AF\_Fro}','\texttt{randQB\_EI}','\texttt{t-SVD}'}, ...
    'Interpreter','latex','Location','best');
xlabel('Prescribed tolerance $\epsilon$','Interpreter','latex');
ylabel('Computed rank $k$','Interpreter','latex');
grid on; set(gca,'FontSize',22);

%% Results table for this matrix: rank k and relative error per method
%  (printed as copy-paste LaTeX rows).
fprintf('\n%%===== Rank and relative error -- exponential (fast) decay =====\n');
fprintf(['$\\epsilon$ & \\multicolumn{2}{c}{SVD} & \\multicolumn{2}{c}{randQB\\_EI}', ...
         ' & \\multicolumn{2}{c}{randQB\\_MF\\_Fro}', ...
         ' & \\multicolumn{2}{c}{randQB\\_AF\\_Fro} \\\\\n']);
fprintf([' & $k$ & rel.\\ err & $k$ & rel.\\ err & $k$ & rel.\\ err', ...
         ' & $k$ & rel.\\ err \\\\ \\hline\n']);
for i = 1:numel(TOL)
    fprintf('%.0e & %d & %.2e & %d & %.2e & %d & %.2e & %d & %.2e \\\\\n', ...
        TOL(i), ...
        RANKS(i,1), REL_ERR(i,1), ...
        RANKS(i,2), REL_ERR(i,2), ...
        RANKS(i,3), REL_ERR(i,3), ...
        RANKS(i,4), REL_ERR(i,4));
end
