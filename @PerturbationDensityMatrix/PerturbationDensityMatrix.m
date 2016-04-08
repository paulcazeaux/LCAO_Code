classdef PerturbationDensityMatrix
    %DENSITYMATRIXPERTURBATION class enables to store relevant information and
    % compute coefficients of the first-order perturbation of a periodic density
    % matrix by a potential of a possibly different periodicity
    
    properties
        
        basis
        % Unperturbed density matrix, storing the periodic system information
        D0
        
        % Perturbation data
        patch
        K_list
        W
        
        % Cached factors
        Qgrid
        Factors
        
        % Numerical parameters
        Nq
    end
    
    methods
        function obj = PerturbationDensityMatrix(D0, patch, V, Nq)
            obj.basis = D0.basis;
            obj.D0 = D0;
            obj.patch = patch;
            obj.K_list = V.K_list;
            obj.W = obj.basis.W(V, obj.patch);
            obj.Nq = Nq;
            
            m = obj.basis.size;
            gap_position = obj.D0.gap_position;
            obj.Factors = zeros(m,m,Nq,Nq,size(obj.K_list,2));
            [obj.Qgrid, eigenvalues_q, eigenvectors_q] = obj.D0.Eigenpairs(Nq); 

            for k = 1:size(obj.K_list,2)
                K = obj.K_list(:,k);
                [~, eigenvalues_q_p, eigenvectors_q_p] = obj.D0.Eigenpairs(Nq, -K); 

                for q1 = 1:obj.Nq
                    for q2 = 1:obj.Nq                
                        Wq = zeros(m,m);
                        for r=1:obj.patch.size
                            Wq = Wq + obj.W(:,:,k,r)*exp(-1i*sum(obj.Qgrid(:,q1,q2).*obj.patch.positions(:,r)));
                        end
                        for n = 1:gap_position
                            for n_p = (gap_position+1):m
                                obj.Factors(:,:,q1,q2,k) = obj.Factors(:,:,q1,q2,k) - eigenvectors_q(:,n,q1,q2)'*Wq*eigenvectors_q_p(:,n_p,q1,q2) ...
                                /(eigenvalues_q_p(n_p,q1,q2)-eigenvalues_q(n,q1,q2)) * eigenvectors_q(:,n,q1,q2)*eigenvectors_q_p(:,n_p,q1,q2)';
                            end
                        end
                    end
                end
            end
        end
        
        function result = Block(obj, R1, R2)
            m = obj.basis.size;
            
            result = zeros(m);
            for q1 = 1:obj.Nq
                for q2 = 1:obj.Nq
                    C = exp(1i*sum(obj.Qgrid(:,q1,q2).*(R1-R2)));
                    for k = 1:size(obj.K_list,2)
                        result = result +  (C*exp(1i*sum(obj.K_list(:,k).* R2)))*obj.Factors(:,:,q1,q2,k)...
                            + (C*exp(-1i*sum(obj.K_list(:,k).*R1)))*obj.Factors(:,:,q1,q2,k)';
                    end
                end
            end
            dq = 1/obj.Nq^2;
            result = dq*result;
        end
    end
end

