% Generate_Sec5_2_Figures  Reproduce Figure 2 of Section 5.2.
%
%   Produces, for each Frobenius-norm test matrix, two panels:
%     (i)  the per-iteration error indicators (t-SVD/EI/MF/AF) at the
%          smallest tolerance, and
%     (ii) the computed rank vs. prescribed tolerance.
%
%   This single script replaces the per-matrix Drive_Comp_ER_* drivers.
%
%   Visual scheme: the two proposed methods get bright, well-separated
%   colors/markers (MF = blue circles, AF = orange squares); the baselines
%   stay subordinate (randQB_EI = gray triangles, reference = black line).

clear; close all;

N = 2^10;                           % matrix size (N x N)
b = 2^5;                            % block size

% Each row: {decay-profile name, number of tolerance decades 10^-1 ... 10^-p}.
specs = { 'exp',     14;            % Matrix 1: exponential decay
          'poly4',   12;            % Matrix 2: polynomial decay (p = 4)
          's_shape', 15 };          % Matrix 3: S-shaped decay

for t = 1:size(specs, 1)
    name = specs{t, 1};
    p    = specs{t, 2};

    A   = Generate_Test_Matrix(N, N, name);
    TOL = 10.^(-1:-1:-p);

    [E, ERRORS, RANKS, ~] = compare_Fro_methods(A, b, TOL);

    figure('Name', sprintf('%s: error indicator', name));
    plot_error_indicator(ERRORS);

    figure('Name', sprintf('%s: rank vs tolerance', name));
    plot_rank_vs_tol(E, RANKS, N);
end


% ===================== local plotting functions =====================

function plot_error_indicator(ERRORS)
% Per-iteration error indicators at the smallest tolerance (last column).
    st = sec52_styles();
    hRef = plot_series(@semilogy, 1:numel(ERRORS{1,end}), ERRORS{1,end}, st.REF);
    hold on
    hEI  = plot_series(@semilogy, 1:numel(ERRORS{2,end}), ERRORS{2,end}, st.EI);
    hMF  = plot_series(@semilogy, 1:numel(ERRORS{3,end}), ERRORS{3,end}, st.MF);
    hAF  = plot_series(@semilogy, 1:numel(ERRORS{4,end}), ERRORS{4,end}, st.AF);
    hold off
    legend([hMF hAF hEI hRef], ...
        {'\texttt{randQB\_MF\_Fro}', '\texttt{randQB\_AF\_Fro}', ...
         '\texttt{randQB\_EI}', 'Exact error'}, ...
        'Interpreter', 'latex', 'Location', 'best');
    xlabel('Iteration $i$', 'Interpreter', 'latex');
    ylabel('Error indicator/estimator', 'Interpreter', 'latex');
    grid on
    set(gca, 'FontSize', 22);
end

function plot_rank_vs_tol(E, RANKS, N)
% Computed rank vs. prescribed tolerance.  randQB_EI returns full rank (N)
% on the tolerances where it fails to converge; those points are omitted.
    st  = sec52_styles();
    mEI = RANKS(:,2) ~= N;
    hRef = plot_series(@semilogx, E,        RANKS(:,1),   st.REF);   % t-SVD
    hold on
    hEI  = plot_series(@semilogx, E(mEI),   RANKS(mEI,2), st.EI);
    hMF  = plot_series(@semilogx, E,        RANKS(:,3),   st.MF);
    hAF  = plot_series(@semilogx, E,        RANKS(:,4),   st.AF);
    hold off
    legend([hMF hAF hEI hRef], ...
        {'\texttt{randQB\_MF\_Fro}', '\texttt{randQB\_AF\_Fro}', ...
         '\texttt{randQB\_EI}', '\texttt{t-SVD}'}, ...
        'Interpreter', 'latex', 'Location', 'best');
    xlabel('Prescribed tolerance $\epsilon$', 'Interpreter', 'latex');
    ylabel('Computed rank $k$', 'Interpreter', 'latex');
    grid on
    set(gca, 'FontSize', 22);
end


% ===================== shared style definitions =====================

function st = sec52_styles()
% Per-series line/marker styles.  Bright, well-separated colors/markers for
% the two proposed methods; subdued styling for the baselines.  The two new
% methods nearly coincide, so they differ in color, marker shape, marker
% size, and line style to stay legible where they overlap.
    st.MF  = make_style([0.00 0.45 0.90], '-',  'o',    15, 3.5, true);   % new: bright blue, large circle
    st.AF  = make_style([0.95 0.45 0.05], '--', 's',     9, 3.5, true);   % new: bright orange, small square
    st.EI  = make_style([0.35 0.35 0.35], '-.', 'd',    12, 3.0, false);  % baseline randQB_EI: dark-gray diamond
    st.REF = make_style([0.00 0.00 0.00], '-',  'none',  6, 3.0, false);  % reference (t-SVD / exact): black line
end

function s = make_style(color, ls, mk, ms, lw, fill)
    s = struct('color', color, 'ls', ls, 'mk', mk, 'ms', ms, 'lw', lw, 'fill', fill);
end

function h = plot_series(plotfun, x, y, s)
% Plot one series with the given style (plotfun is @semilogy or @semilogx).
% MF is drawn with a larger marker and AF with a smaller one, so the two
% nearly-coincident curves remain distinguishable where they overlap.
    h = plotfun(x, y, 'Color', s.color, 'LineStyle', s.ls, 'Marker', s.mk, ...
                'MarkerSize', s.ms, 'LineWidth', s.lw);
    if s.fill && ~strcmp(s.mk, 'none')
        set(h, 'MarkerFaceColor', s.color);
    end
end
