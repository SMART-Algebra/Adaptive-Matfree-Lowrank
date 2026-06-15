function [Q, B, k, Err_ind, Exact_Err]  = randQB_tol_S_rel_err_rank(A, b, epsilon, pow_its)
    
    %% INPUT: Matrix A, 
    %                 block size b, 
    %                 tolerance epsilon such that 0 < epsilon < 1 
    %                 number of power iterations to apply to A pow_its

    %% OUTPUT: Orthonormal matrix Q 
    %                     Full matrix B = Q^* * A
    %                     rank k
    %                     such that ||A - QB||_2 < epsilon ||A||_2 with
    %                     high probability

    %% Execution
    % Initialize empty arrays
    Q = [];             B = [];   

    n = size(A,2);
    %safety_factor = 10*sqrt(2/pi); 
    safety_factor = 1; 
    for i = 1:(ceil(n/b) + 1)
        % Draw random matrix
        Omega_i = randn(n,b);
        
        % Approximate ||A||_2
        if i ==1
            x = Omega_i(:,1);
            x = x / norm(x);

            for j = 1:pow_its
                z = A' * (A * x);
                x = z / norm(z);
            end % end for 
            norm_A_est = norm(A*x);

        end % end if

        
        
        % Get random projection
        Z_i = A * Omega_i;

        %% Compute orthonormal basis 
        if isempty(Q)
            Y_i = Z_i;
            [Q_i, ~] = qr(Y_i, 'econ');
            M = A;
            Exact_Err(i) = norm(M)^2;
        else
            % Project residual onto random subspace
            Y_i = Z_i - Q * (B * Omega_i);
            [Q_Y_i, ~] = qr(Y_i, 'econ');

            % Reorthogonalize 
            [Q_i, ~] = qr(Q_Y_i - Q * (Q' * Q_Y_i), 'econ');
            M = A - Q*B;
            Exact_Err(i) = norm(M)^2;
        end  % end if   

        %% Compute block B_i of B
        B_i = Q_i' * A;

        % Track error in low-rank approximation 
        Err_ind(i) = norm(B_i)^2;
        
        if  Err_ind(i) < (epsilon*norm_A_est / safety_factor)^2
            %% Rank - refinement
            % Extract last appended block from Q and B 
            Q_temp = Q(:,(i-2)*b+1:(i-1)*b);
            B_temp = B((i-2)*b+1:(i-1)*b, :);

            % Compute SVD of last block of B
             [U_B_i, S_B_i, V_B_i] = svd(B_temp, 'econ');

            % Determine which singular values to keep
             sing_vals_B_i = diag(S_B_i);
             keep = find(sing_vals_B_i >= epsilon);

            % Adjust last appended blocks
             Q_temp = Q_temp * U_B_i(:,keep);
             B_temp = S_B_i(keep, keep) * V_B_i(:, keep)';
             
             % Complete rank refinement 
             Q = [Q(:, 1:(i-2)*b), Q_temp(:, keep)];
             B = [B(1:(i-2)*b, :); B_temp(keep, :)];
             k = size(Q,2);
            return 
        else
            % Append block to Q and B
            Q = [Q, Q_i];
            B = [B; B_i];

            % Return matrix if no compression is needed
            if size(Q,2)>= size(A,2)
                % Q = eye(size(A,1));
                % B = A;
                k = size(Q,2);
                disp('Matrix is full rank OR choose larger tolerance.')
                fprintf('Tolerance: %.1e\n', epsilon)
                fprintf('Smallest error: %.4e\n', min(sqrt(Err_ind)))
                fprintf('Happened on iteration %.0f\n', find(Err_ind==min(Err_ind)))
                fprintf('rank: %.0f\n',k)
                disp('-----------------------------------------------------------------------')
                return
            end % end if
        end % end if 
    end % end for 
   
end % end function 