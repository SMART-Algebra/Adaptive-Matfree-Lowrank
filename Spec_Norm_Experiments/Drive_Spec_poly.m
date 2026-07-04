% Drive_Spec_poly  Section 5.3 results for the slowly-decaying spectrum.
%
%   Test matrix: Matrix 2, polynomial decay  sigma_i = 1/i^2.
%   Produces the two panels for this matrix:
%     (a) computed rank vs. prescribed tolerance, and
%     (b) per-iteration error indicators (at the smallest tolerance),
%   comparing t-SVD, randQB_HMT, randQB_MF_Spec, and randQB_AF_Spec.
%
%   The heavy comparison is delegated to compare_Spec_methods; the plotting is
%   kept here so this matrix can be tuned independently of the other drivers.

clear; close all;

%% Problem setup
N   = 2^10;                     % matrix size (N x N)
b   = 2^5;                      % block size
TOL = sort([ ...                % refined tolerance grid down to 5e-7
    1e-1 5e-2 2e-2 ...
    1e-2 5e-3 2e-3 ...
    1e-3 5e-4 2e-4 ...
    1e-4 5e-5 2e-5 1e-5 5e-6 1e-6 5e-7], 'descend');

A = Generate_Test_Matrix(N, N, 'poly2');

%% Run the comparison
[ERR_IND, RANKS, REL_ERR] = compare_Spec_methods(A, b, TOL);

%% Line/marker styles (bright + distinct for the new methods; black reference)
MF  = {'Color','b','LineStyle','-.', 'Marker','o','MarkerSize',23, ...
       'LineWidth',3.5};   % randQB_MF_Spec (hollow marker so overlaps show through)
AF  = {'Color','r','LineStyle','--','Marker','^','MarkerSize',15,  ...
       'LineWidth',3.5,'MarkerFaceColor',[0.85 0.10 0.10]};   % randQB_AF_Spec (solid red triangle)
HMT = {'Color','g','LineStyle','-','Marker','*','MarkerSize',30, ...
       'LineWidth',3.0};                                      % randQB_HMT
REF = {'Color',[0 0 0],         'LineStyle','-', 'Marker','none','LineWidth',3.0}; % t-SVD / exact

% (a) Computed rank vs. prescribed tolerance
f_rank = figure('Name','poly: rank vs tolerance');
hRef = semilogx(TOL, RANKS(:,1), REF{:}); hold on   % t-SVD
hMF  = semilogx(TOL, RANKS(:,3), MF{:});            % randQB_MF_Spec
hAF  = semilogx(TOL, RANKS(:,4), AF{:});            % randQB_AF_Spec
hHMT = semilogx(TOL, RANKS(:,2), HMT{:});           % randQB_HMT
hold off
legend([hMF hAF hHMT hRef], ...
    {'\texttt{randQB\_MF\_Spec}','\texttt{randQB\_AF\_Spec}','\texttt{randQB\_HMT}','\texttt{t-SVD}'}, ...
    'Interpreter','latex','Location','best');
xlabel('Prescribed tolerance $\epsilon$','Interpreter','latex');
ylabel('Computed rank $k$','Interpreter','latex');
grid on; set(gca,'FontSize',22);

% (b) Error indicators at the smallest tolerance
f_err = figure('Name','poly: error indicator');
hMF  = semilogy(1:numel(ERR_IND{end,3}), ERR_IND{end,3}, MF{:});  hold on
hAF  = semilogy(1:numel(ERR_IND{end,4}), ERR_IND{end,4}, AF{:});
hHMT = semilogy(1:numel(ERR_IND{end,2}), ERR_IND{end,2}, HMT{:});
hRef = semilogy(1:numel(ERR_IND{end,1}), ERR_IND{end,1}, REF{:});
hold off
legend([hMF hAF hHMT hRef], ...
    {'\texttt{randQB\_MF\_Spec}','\texttt{randQB\_AF\_Spec}','\texttt{randQB\_HMT}','Exact error'}, ...
    'Interpreter','latex','Location','best');
xlabel('Iteration $i$','Interpreter','latex');
ylabel('Error indicator','Interpreter','latex');
grid on; set(gca,'FontSize',22);

%% Save the two panels as EPS (copy these into the manuscript's FIGURES/ folder)
print(f_rank, 'rank_Spec_poly',    '-depsc');
print(f_err,  'err_ind_Spec_poly', '-depsc');

%% Results table for this matrix: rank k and relative error per method
%  (printed as copy-paste LaTeX rows).
fprintf('\n%%===== Rank and relative error -- polynomial (slow) decay =====\n');
fprintf(['$\\epsilon$ & \\multicolumn{2}{c}{SVD} & \\multicolumn{2}{c}{randQB\\_HMT}', ...
         ' & \\multicolumn{2}{c}{randQB\\_MF\\_Spec}', ...
         ' & \\multicolumn{2}{c}{randQB\\_AF\\_Spec} \\\\\n']);
fprintf([' & $k$ & rel.\\ err & $k$ & rel.\\ err & $k$ & rel.\\ err', ...
         ' & $k$ & rel.\\ err \\\\ \\hline\n']);
for i = 1:numel(TOL)
    % randQB_HMT does not fail; it simply overestimates the rank (sometimes up
    % to full rank).  Print the raw rank and relative error for every method.
    fprintf('%.0e & %d & %.2e & %d & %.2e & %d & %.2e & %d & %.2e \\\\\n', ...
        TOL(i), ...
        RANKS(i,1), REL_ERR(i,1), ...
        RANKS(i,2), REL_ERR(i,2), ...
        RANKS(i,3), REL_ERR(i,3), ...
        RANKS(i,4), REL_ERR(i,4));
end
