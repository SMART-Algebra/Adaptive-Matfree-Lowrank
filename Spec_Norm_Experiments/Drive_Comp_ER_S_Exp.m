clear all
close all

%% Generate A
N = 2^10;
A = Generate_Test_Matrix(N, N, 'exp');

%% Set block size
b = 2^3;

%% Generate array of tolerance levels
p = 14;
TOL = 10.^(-1:-1:-p);

% Choose number of power iterations 
pow_its = 3;


[E, ERRORS, RANKS, REL_ERR] = comp_err_rank_S2(A, b, TOL, pow_its);

RANKS_2 =  RANKS(:,2);
%REL_ERR_2 = REL_ERR(:,2);
%REL_ERR_2(14) = [ ]
RANKS_2(14) = [ ];

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
legend({'\texttt{randQB\_AF\_Fro}','\texttt{randQB\_MF\_Fro}','\texttt{A.R.F}', 'Exact error'}, 'Interpreter','latex', 'Location','best');
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
legend({'\texttt{randQB\_AF\_Fro}','\texttt{randQB\_MF\_Fro}', '\texttt{A.R.F}' ,'\texttt{t-SVD}'}, 'Interpreter','latex', 'Location','best');
xlabel( 'Prescribed tolerance $\epsilon$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('Computed rank $k$', 'FontSize', 24, 'Interpreter', 'latex');
grid on 
set(gca, 'FontSize', 24);




%% Plotting 
% R2 =  RANKS(:,2);
% REL_ERR_2 = REL_ERR(:,2);
% REL_ERR_2(14) = [ ]
% R2(14) = [ ];
% figure(1)
% % Error indicators
% semilogy(1:length(ERRORS{3,1}),sqrt(ERRORS{3,1}),'k', 'LineStyle', ':', 'Marker', 's', 'MarkerSize', 12, 'LineWidth',3)
% hold on
% semilogy(1:length(ERRORS{2,1}),(ERRORS{2,1}),'b', 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 12, 'LineWidth',3)
% semilogy(1:length(ERRORS{1,1}), sqrt(ERRORS{1,1}),'Color',[0.0, 0.8, 0.0], 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 12, 'LineWidth',3)
% hold off
% legend({'\texttt{randQB\_MF\_Spec}','\texttt{Adaptive Range Finder}', 'Exact error'}, 'Interpreter','latex', 'Location','best', 'FontSize', 18);
% xlabel('Iteration $i$', 'FontSize', 24, 'Interpreter', 'latex');
% grid on 
% set(gca, 'FontSize', 24);





% figure(1)
% 
% semilogy(1:length(ERRORS{3,1}),sqrt(ERRORS{3,1}), ...
%     'k','LineStyle',':','Marker','s', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% hold on
% 
% semilogy(1:length(ERRORS{4,1}),ERRORS{4,1}, ...
%     'Color',[0.5 0.5 0.5], ...
%     'LineStyle','-.', ...
%     'Marker','o', ...
%     'MarkerSize',12, ...
%     'LineWidth',3)
% 
% semilogy(1:length(ERRORS{2,1}),ERRORS{2,1}, ...
%     'b','LineStyle',':','Marker','^', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% semilogy(1:length(ERRORS{1,1}),sqrt(ERRORS{1,1}), ...
%     'Color',[0.0,0.8,0.0], ...
%     'LineStyle','--','Marker','x', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% hold off
% 
% legend({'\texttt{randQB\_MF\_Spec}', ...
%          '\texttt{randQB\_AF\_Spec}', ...
%          '\texttt{Adaptive Range Finder}', ...
%          'Exact error'}, ...
%          'Interpreter','latex', ...
%          'Location','best', ...
%          'FontSize',18);
% xlabel('Iteration $i$', 'FontSize', 24, 'Interpreter', 'latex');
% grid on 
% set(gca, 'FontSize', 24);
% % figure(2)
% % % Ranks
% % semilogx(E, RANKS(:,3),'k', 'LineStyle', ':', 'Marker', 's', 'MarkerSize', 12, 'LineWidth',3)
% % hold on
% % semilogx(E(1:length(R2)), R2,'b', 'LineStyle', ':', 'Marker', '^', 'MarkerSize', 12, 'LineWidth',3)
% % semilogx(E, RANKS(:,1),'r', 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 12, 'LineWidth',3)
% % hold off
% % legend({'\texttt{randQB\_MF\_Spec}','\texttt{Adaptive Range Finder}' , '\texttt{t-SVD}'}, 'Interpreter','latex', 'Location','best', 'FontSize', 18);
% % xlabel( 'Prescribed tolerance $\epsilon$', 'FontSize', 24, 'Interpreter', 'latex');
% % ylabel('Computed rank $k$', 'FontSize', 24, 'Interpreter', 'latex');
% % grid on 
% % set(gca, 'FontSize', 24);
% 
% figure(2)
% 
% semilogx(E, RANKS(:,3), ...
%     'k','LineStyle',':','Marker','s', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% hold on
% 
% semilogx(E, RANKS(:,4), ...
%     'Color',[0.5 0.5 0.5], ...
%     'LineStyle','-.', ...
%     'Marker','o', ...
%     'MarkerSize',12, ...
%     'LineWidth',3)
% 
% semilogx(E(1:length(R2)),R2, ...
%     'b','LineStyle',':','Marker','^', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% semilogx(E,RANKS(:,1), ...
%     'r','LineStyle','--','Marker','x', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% hold off
% legend({'\texttt{randQB\_MF\_Spec}', ...
%          '\texttt{randQB\_AF\_Spec}', ...
%          '\texttt{Adaptive Range Finder}', ...
%          '\texttt{t-SVD}'}, ...
%          'Interpreter','latex', ...
%          'Location','best', ...
%          'FontSize',18);
% xlabel( 'Prescribed tolerance $\epsilon$', 'FontSize', 24, 'Interpreter', 'latex');
% ylabel('Computed rank $k$', 'FontSize', 24, 'Interpreter', 'latex');
% grid on 
% set(gca, 'FontSize', 24);
% 
% 
% 
% 
% 
% figure(3)
% 
% semilogx(E, REL_ERR(:,3), ...
%     'k','LineStyle',':','Marker','s', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% hold on
% 
% semilogx(E, REL_ERR(:,4), ...
%     'Color',[0.5 0.5 0.5], ...
%     'LineStyle','-.', ...
%     'Marker','o', ...
%     'MarkerSize',12, ...
%     'LineWidth',3)
% 
% semilogx(REL_ERR_2,...
%     'b','LineStyle',':','Marker','^', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% semilogx(E,REL_ERR(:,1), ...
%     'r','LineStyle','--','Marker','x', ...
%     'MarkerSize',12,'LineWidth',3)
% 
% hold off
% legend({'\texttt{randQB\_MF\_Spec}', ...
%          '\texttt{randQB\_AF\_Spec}', ...
%          '\texttt{Adaptive Range Finder}', ...
%          '\texttt{t-SVD}'}, ...
%          'Interpreter','latex', ...
%          'Location','best', ...
%          'FontSize',18);
% xlabel( 'Prescribed tolerance $\epsilon$', 'FontSize', 24, 'Interpreter', 'latex');
% ylabel('Computed rank $k$', 'FontSize', 24, 'Interpreter', 'latex');
% grid on 
% set(gca, 'FontSize', 24);




