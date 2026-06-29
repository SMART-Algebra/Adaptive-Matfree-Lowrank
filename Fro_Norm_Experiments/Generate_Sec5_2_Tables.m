% Generate_Sec5_2_Tables  Reproduce the data in Tables 2-4 of Section 5.2.
%
%   For each Frobenius-norm test matrix, prints the relative approximation
%   errors of t-SVD, randQB_EI, randQB_MF_Fro, and randQB_AF_Fro across the
%   prescribed tolerances, formatted as copy-paste LaTeX table rows.

clear; close all;

N = 2^10;                           % matrix size (N x N)
b = 2^5;                            % block size

% Each row: {decay-profile name, number of tolerance decades, table caption}.
specs = { 'exp',     14, 'Exponential Decay';     % Table 2
          'poly4',   12, 'Polynomial Decay';      % Table 3
          's_shape', 15, 'S-shaped Decay' };      % Table 4

for t = 1:size(specs, 1)
    name  = specs{t, 1};
    p     = specs{t, 2};
    label = specs{t, 3};

    A   = Generate_Test_Matrix(N, N, name);
    TOL = 10.^(-1:-1:-p);

    [~, ~, REL_ERR] = compare_Fro_methods(A, b, TOL);

    fprintf('\n%%===== Relative errors -- %s =====\n', label);
    fprintf(['$\\epsilon$ & SVD & randQB\\_EI & randQB\\_MF\\_Fro ', ...
             '& randQB\\_AF\\_Fro \\\\ \\hline\n']);
    for i = 1:numel(TOL)
        fprintf('%.0e & %.4e & %.4e & %.4e & %.4e \\\\\n', ...
            TOL(i), REL_ERR(i,1), REL_ERR(i,2), REL_ERR(i,3), REL_ERR(i,4));
    end
end
