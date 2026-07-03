% Drive_Fro_poly  Section 5.2 results for the slowly-decaying (polynomial) spectrum.
%
%   Test matrix: Matrix 2, polynomial decay  sigma_i = 1/i^4.
%   Produces the two Figure panels for this matrix:
%     (a) per-iteration error indicators (at the smallest tolerance), and
%     (b) computed rank vs. prescribed tolerance,
%   comparing t-SVD, randQB_EI, randQB_MF_Fro, and randQB_AF_Fro, plus the
%   results table.
%
%   The heavy comparison is delegated to compare_Fro_methods; the plotting is
%   kept here so this matrix can be tuned independently of the other drivers.

clear; close all;

%% Problem setup
N   = 2^10;                     % matrix size (N x N)
b   = 32;                       % block size
p   = 12;                       % tolerances 10^-1 ... 10^-p
TOL = 10.^(-1:-1:-p);

A = Generate_Test_Matrix(N, N, 'poly4');

%% Run the comparison
[ERR_IND, RANKS, REL_ERR] = compare_Fro_methods(A, b, TOL);

%% Line/marker styles (bright + distinct for the new methods; black reference)
MF  = {'Color','b','LineStyle','-.', 'Marker','o','MarkerSize',23, ...
       'LineWidth',3.5};   % randQB_MF_Fro (hollow marker so overlapped markers show through)
AF  = {'Color','r','LineStyle','--','Marker','^','MarkerSize',15,  ...
       'LineWidth',3.5,'MarkerFaceColor',[0.85 0.10 0.10]};   % randQB_AF_Fro (solid red triangle)
EI  = {'Color','g','LineStyle','-','Marker','*','MarkerSize',30, ...
       'LineWidth',3.0};                                      % randQB_EI
REF = {'Color',[0 0 0],         'LineStyle','-', 'Marker','none','LineWidth',3.0}; % t-SVD / exact

% (a) Error indicators at the smallest tolerance
f_err = figure('Name','poly: error indicator');
hMF  = semilogy(1:numel(ERR_IND{end,3}), ERR_IND{end,3}, MF{:}); hold on
hAF  = semilogy(1:numel(ERR_IND{end,4}), ERR_IND{end,4}, AF{:});
hEI  = semilogy(1:numel(ERR_IND{end,2}), ERR_IND{end,2}, EI{:});
hRef = semilogy(1:numel(ERR_IND{end,1}), ERR_IND{end,1}, REF{:});
hold off
legend([hMF hAF hEI hRef], ...
    {'\texttt{randQB\_MF\_Fro}','\texttt{randQB\_AF\_Fro}','\texttt{randQB\_EI}','Exact error'}, ...
    'Interpreter','latex','Location','best');
xlabel('Iteration $i$','Interpreter','latex');
ylabel('Error indicator','Interpreter','latex');
grid on; set(gca,'FontSize',22);

% (b) Computed rank vs. prescribed tolerance
mEI = RANKS(:,2) ~= N;          % drop tolerances where randQB_EI returns full rank
f_rank = figure('Name','poly: rank vs tolerance');
hRef = semilogx(TOL,      RANKS(:,1),   REF{:}); hold on   % t-SVD
hMF  = semilogx(TOL,      RANKS(:,3),   MF{:});            % randQB_MF_Fro
hAF  = semilogx(TOL,      RANKS(:,4),   AF{:});    % randQB_AF_Fro
hEI  = semilogx(TOL(mEI), RANKS(mEI,2), EI{:});            % randQB_EI
hold off
legend([hMF hAF hEI hRef], ...
    {'\texttt{randQB\_MF\_Fro}','\texttt{randQB\_AF\_Fro}','\texttt{randQB\_EI}','\texttt{t-SVD}'}, ...
    'Interpreter','latex','Location','best');
xlabel('Prescribed tolerance $\epsilon$','Interpreter','latex');
ylabel('Computed rank $k$','Interpreter','latex');
grid on; set(gca,'FontSize',22);

%% Save the two panels as EPS (copy these into the manuscript's FIGURES/ folder)
print(f_err,  'err_ind_Fro_poly', '-depsc');
print(f_rank, 'rank_Fro_poly',    '-depsc');

%% Results table for this matrix: rank k and relative error per method
%  (printed as copy-paste LaTeX rows).
fprintf('\n%%===== Rank and relative error -- polynomial decay =====\n');
fprintf(['$\\epsilon$ & \\multicolumn{2}{c}{SVD} & \\multicolumn{2}{c}{randQB\\_EI}', ...
         ' & \\multicolumn{2}{c}{randQB\\_MF\\_Fro}', ...
         ' & \\multicolumn{2}{c}{randQB\\_AF\\_Fro} \\\\\n']);
fprintf([' & $k$ & rel.\\ err & $k$ & rel.\\ err & $k$ & rel.\\ err', ...
         ' & $k$ & rel.\\ err \\\\ \\hline\n']);
for i = 1:numel(TOL)
    % randQB_EI fails when it cannot deliver a low-rank approximation meeting the
    % tolerance: either its error exceeds the tolerance, or it returns full rank
    % (no compression). Mark such cases 'fail' in red across its two columns
    % (needs \usepackage{xcolor} in the manuscript).
    if REL_ERR(i,2) > TOL(i) || (RANKS(i,2) == N && RANKS(i,1) < N)
        ei_cell = '\multicolumn{2}{c|}{\textcolor{red}{fail}}';
    else
        ei_cell = sprintf('%d & %.2e', RANKS(i,2), REL_ERR(i,2));
    end
    fprintf('%.0e & %d & %.2e & %s & %d & %.2e & %d & %.2e \\\\\n', ...
        TOL(i), ...
        RANKS(i,1), REL_ERR(i,1), ...
        ei_cell, ...
        RANKS(i,3), REL_ERR(i,3), ...
        RANKS(i,4), REL_ERR(i,4));
end
