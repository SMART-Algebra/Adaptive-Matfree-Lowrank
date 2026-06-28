clear all
close all

%% Generate A
N = 2^10;
A = Generate_Test_Matrix(N, N, 6);

%% Set block size
b = 2^5;

%% Generate array of tolerance levels
p = 15;
TOL = 10.^(-1:-1:-p);

%% Execution 
[E, ERRORS, RANKS, REL_ERR] = comp_err_rank_F2(A, b, TOL);


%% Plotting 
RANKS_2 = RANKS(:,2);
RANKS_2(9:end) = [ ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% semilogy(1:length(ERRORS{3,1}),(ERRORS{3,1}),'k', 'LineStyle', ':', 'Marker', 'o', 'MarkerSize', 12, 'LineWidth',3)
% hold on
% semilogy(1:length(ERRORS{4,1}),ERRORS{4,1}, ...
%     'Color',[0.5 0.5 0.5], ...
%     'LineStyle','-.', ...
%     'Marker','o', ...
%     'MarkerSize',12, ...
%     'LineWidth',3)
% semilogy(1:length(ERRORS{2,1}),sqrt(ERRORS{2,1}),'b', 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 12, 'LineWidth',3)
% semilogy(1:length(ERRORS{1,1}), (ERRORS{1,1}),'Color',[0.0, 0.8, 0.0], 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 12, 'LineWidth',3)
% hold off
% legend({'\texttt{randQB\_MF\_Fro}','\texttt{randQB\_AF\_Fro}','\texttt{randQB\_EI}', 'Exact error'}, 'Interpreter','latex', 'Location','best');
% xlabel('Iteration $i$', 'FontSize', 24, 'Interpreter', 'latex');
% grid on 
% set(gca, 'FontSize', 24);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% semilogx(E, RANKS(:,3),'k', 'LineStyle', ':', 'Marker', 'o', 'MarkerSize', 12, 'LineWidth',3)
% hold on
% semilogx(E, RANKS(:,4), ...
%     'Color',[0.5 0.5 0.5], ...
%     'LineStyle','-.', ...
%     'Marker','o', ...
%     'MarkerSize',12, ...
%     'LineWidth',3)
% semilogx(E(1:length(RANKS_2)),RANKS_2,'b', 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 12, 'LineWidth',3)
% semilogx(E, RANKS(:,1),'r', 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 12, 'LineWidth',3)
% hold off
% legend({'\texttt{randQB\_MF\_Fro}','\texttt{randQB\_AF\_Fro}', '\texttt{randQB\_EI}' ,'\texttt{t-SVD}'}, 'Interpreter','latex', 'Location','best');
% xlabel( 'Prescribed tolerance $\epsilon$', 'FontSize', 24, 'Interpreter', 'latex');
% ylabel('Computed rank $k$', 'FontSize', 24, 'Interpreter', 'latex');
% grid on 
% set(gca, 'FontSize', 24);


figure(1)
semilogy(1:length(ERRORS{4,end}),ERRORS{4,end}, ...
    'Color',[0.5 0.5 0.5], ...
    'LineStyle','-.', ...
    'Marker','o', ...
    'MarkerSize',12, ...
    'LineWidth',3)
hold on
semilogy(1:length(ERRORS{3,end}),(ERRORS{3,end}),'k', 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 12, 'LineWidth',3)

semilogy(1:length(ERRORS{2,end}),ERRORS{2,end},'b', 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 12, 'LineWidth',3)
semilogy(1:length(ERRORS{1,end}), (ERRORS{1,end}),'Color',[0.0, 0.8, 0.0], 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 12, 'LineWidth',3)
hold off
legend({'\texttt{randQB\_AF\_Fro}','\texttt{randQB\_MF\_Fro}','\texttt{randQB\_EI}', 'Exact error'}, 'Interpreter','latex', 'Location','best');
xlabel('Iteration $i$', 'FontSize', 24, 'Interpreter', 'latex');
grid on 
set(gca, 'FontSize', 24);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
semilogx(E, RANKS(:,4), ...
    'Color',[0.5 0.5 0.5], ...
    'LineStyle',':', ...
    'Marker','o', ...
    'MarkerSize',16, ...
    'LineWidth',3)
hold on
semilogx(E, RANKS(:,3),'k', 'LineStyle','none','Marker', 's', 'MarkerSize', 12, 'LineWidth',3)


semilogx(E(1:length(RANKS_2)),RANKS_2,'b', 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 12, 'LineWidth',3)
semilogx(E, RANKS(:,1),'r', 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 12, 'LineWidth',3)
hold off
legend({'\texttt{randQB\_AF\_Fro}','\texttt{randQB\_MF\_Fro}', '\texttt{randQB\_EI}' ,'\texttt{t-SVD}'}, 'Interpreter','latex', 'Location','best');
xlabel( 'Prescribed tolerance $\epsilon$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('Computed rank $k$', 'FontSize', 24, 'Interpreter', 'latex');
grid on 
set(gca, 'FontSize', 24);



