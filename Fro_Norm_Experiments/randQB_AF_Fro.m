function [Q, k, Err_ind, Exact_err, S] = randQB_AF_Fro(A, b, epsilon)
% randQB_AF_Fro  Adjoint-free, matrix-free adaptive low-rank range finder in
%                the Frobenius norm (Algorithm 3.3).
    Q = [ ];
    [m ,n] = size(A);
    min_mn = min(m,n);
    Err_ind = zeros(ceil(min_mn/b),1);
    Exact_err = zeros(ceil(min_mn/b),1);
    S = zeros(ceil(min_mn/b),1);
    sum_Z_i = 0;

    for i = 1:ceil(min_mn/b)
        
        % Draw random matrix
        Omega_i = (1/sqrt(b)) * randn(n,b);

        % Get random projection
        Z_i = A * Omega_i;
        sum_Z_i = hypot(sum_Z_i, norm(Z_i, 'fro'));
        S(i) = sum_Z_i / sqrt(i);
        % Implicit residual 
        if isempty(Q)
            M = A;
            Exact_err(i) =  norm(M, 'fro');
            Y_i = Z_i;
        else
            M = A - Q * (Q' *A);
            Exact_err(i) =  norm(M, 'fro');
            Y_i = Z_i  - Q * (Q'*Z_i);

        end % end if 
        
        %Check error 
        Err_ind(i) = norm(Y_i, 'fro');
        
        if Err_ind(i) < epsilon * S(i)
            rel_err_before = norm(A - Q*(Q'*A),'fro') / norm(A,'fro');
            alpha = Err_ind(i);
            tau = epsilon * S(i);
            size_Q_b4_rr = size(Q);
            

            Q_end = Q(:, (i-2)*b+1:(i-1)*b);
            B_i_1_tilde = Q_end' * Z_i;

            % Compute  norms of each row in B_i_1_tilde
            w = vecnorm(B_i_1_tilde, 2, 2);
            [w_sorted, pi] = sort(w,'ascend');

            % Initialize empty set 
            r = [ ];

            alpha_hist = zeros(b,1);

            for j = 1:b
                % Update estimator
                alpha_new = sqrt(alpha^2 + w_sorted(j)^2);

                if alpha_new >= tau
                    break
                end % end if 


                alpha = alpha_new;
                alpha_hist(j) = alpha;
                r = [r, pi(j)];
            end % end for
            
            % Keep non-discarded indices
            s = setdiff(1:b, r, 'stable');
            alpha_hist;
            Q_i_s = Q_end(:,s);

             
            Q = [Q(:, 1:(i-2)*b) Q_i_s];
            size_Q_after = size(Q,2);


            % Remove unused entries
            Err_ind = Err_ind(1:i);
            Exact_err = Exact_err(1:i);
            S = S(1:i);
            
            %% Rank-refinement 



            k = size(Q,2);
            break

        end % end if 

        [Q_i, ~] = qr(Y_i, 'econ');

        % Reorthogonalization : Q_i = Q_i - Q (Q' * Q_i)
        if isempty(Q)
            Q_i = Q_i;
        else
            Q_i = Q_i - Q * (Q'*Q_i);
            [Q_i, ~] = qr(Q_i, 'econ');
        end % end if 
    
      
        
        Q = [Q Q_i];

    end % end for 

    if i==ceil(min_mn/b) && size(Q,2) >= min_mn
        k = min_mn;
    end

end % end function 