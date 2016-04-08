classdef PeriodicDensityMatrix
    %PERIODICDENSITYMATRIX class enables to store relevant information and
    % compute coefficients of the density matrix of a given tight-binding
    % periodic system.

    
    properties
        % Geometric properties
        lattice
        basis
        patch
        
        % System parameterization
        h
        s
        gap_position
        
        % Numerical data
        Nq
        Cached
        Qgrid
        eigenvalues_q
        eigenvectors_q
    end
    
    methods
        function obj = PeriodicDensityMatrix(lattice, basis, patch, h, s, gap_position, Nq)
            if nargin == 7
                obj.lattice = lattice;
                obj.basis = basis;
                obj.patch = patch;
                obj.h = h;
                obj.s = s;
                obj.gap_position = gap_position;
                obj.Nq = Nq;
                obj.Cached = 0;
                [obj.Qgrid, obj.eigenvalues_q, obj.eigenvectors_q] = Eigenpairs(obj, Nq);
                obj.Cached = 1;
            end
        end
        
        [Qgrid, eigenvalues_q, eigenvectors_q] = Eigenpairs(obj, Nq, k);
        
        function result = Block(obj, R1, R2)
            result = zeros(obj.basis.size);
            for q1 = 1:obj.Nq
                for q2 = 1:obj.Nq
                    factor = exp(1i*sum(obj.Qgrid(:,q1,q2).*(R1-R2)));
                    for n = 1:obj.gap_position
                        result = result + factor*obj.eigenvectors_q(:,n,q1,q2)*obj.eigenvectors_q(:,n,q1,q2)';
                    end
                end
            end
            dq = 1/obj.Nq^2;
            result = dq*result;
        end
    end
end

