function [Q, k, Err_ind, Exact_err] = randQB_AF_Spec_New(A, b, epsilon)
    Q = [];
    [m,n] = size(A);

    max_iter = ceil(min(m,n)/b);

    R = cell(max_iter,1);
    R_reorth = R;
    Q_bar = R;
    Q_orth = R;
    
    Err_ind = zeros(max_iter,1);

    for i = 1:max_iter
        
        %% Generate Gaussian sample
        Omega_i = randn(n,b);

        %% Apply matrix
        Z_i = A * Omega_i;

        if isempty(Q)
            [Q_bar{i,1}, R{i,1}] = qr(Z_i, 'econ');
            Q_orth{i,1} = Q_bar{i,1};
        else
            %% Implicit residual
            Y_i = Z_i - Q*(Q'*Z_i);

            %% Local QB
            [Q_bar{i,1}, R{i,1}] = qr(Y_i, 'econ');

            %% Reorthogonalization
            Q_i = Q_bar{i,1} - Q*(Q'*Q_bar{i,1});
            [Q_orth{i,1},~] = qr(Q_i, 'econ');
        end % end if 

        Err_ind(i) = norm(R{i,1})/sqrt(b);
        % Err_ind(i) = norm(Q_orth{i,1}'*Z_i)/sqrt(b);
        
        if Err_ind(i) >= epsilon * Err_ind(1)
            %% Append basis
            Q = [Q Q_orth{i,1}];
            M = A - Q *(Q'*A);
            Exact_err(i) = norm(M);
        end
        if Err_ind(i) < epsilon * Err_ind(1)
            Err_ind = Err_ind(1:i);
            % disp("BELOW TOLERANCE")
            err_b4_refinement = norm(A - Q*(Q'*A))/norm(A);
            Exact_err = [norm(A) Exact_err];
            if i == 1

                Q_test = Q_orth{i,1};
                R_test = R{i,1};

                [U,S,~] = svd(R_test);

                S = S/sqrt(b);
                svals = diag(S);

                tol_rank = epsilon * Err_ind(1);

                r = find(svals >= tol_rank);

                U_keep = U(:,r);

                Q = Q_test * U_keep;
                k = size(Q,2);
                return
            end % end if 
            i;
            sz_Q_b4_rr = size(Q);
            Q_test = Q_orth{i-1,1};

            R_test = R{i-1,1};
            [U,S,~] = svd(R_test);
            % [U,S,~] = svd(Q_test'*Z_i);
            k = size(Q,2);
            S = S/sqrt(b);
            svals = diag(S);
            
            % tol_rank = epsilon*Err_ind(1);
            tol = epsilon * Err_ind(1);

            r = find(svals >= tol);
            k = size(Q,2);
            U_keep = U(:,r);
            size(U_keep)
            Q = [Q(:, 1:(i-2)*b) Q_test*U_keep];
            k = size(Q,2);


            return
        end % end if 

        if size(Q,2) >= min(m,n) && Err_ind(i)>= epsilon * Err_ind(1)
            k = min(m,n);
            return
        end
    end % end for






    end % end function 