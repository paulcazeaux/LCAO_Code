function [Qgrid, eigenvalues_q, eigenvectors_q] = Eigenpairs(obj, Nq, k)

    if (nargin == 2 && obj.Cached && Nq == obj.Nq)
        Qgrid = obj.Qgrid;
        eigenvalues_q = obj.eigenvalues_q;
        eigenvectors_q = obj.eigenvectors_q;
        return;
    elseif nargin == 2
        Qgrid = obj.lattice.Qgrid(Nq);
    elseif nargin == 3
        if all(k == 0)
            [Qgrid, eigenvalues_q, eigenvectors_q] = obj.Eigenpairs(Nq);
            return;
        else
            Qgrid = obj.lattice.Qgrid(Nq, k);
        end
    end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Band structure calculation
    eigenvectors_q = zeros(obj.basis.size,obj.basis.size, Nq, Nq);
    eigenvalues_q = zeros(obj.basis.size, Nq, Nq);
    
    
    for q1 = 1:Nq
        for q2 = 1:Nq
            s_q = zeros(obj.basis.size);
            h_q = zeros(obj.basis.size);
            
            % Bloch transform
            Bloch = exp(-1i*(Qgrid(:,q1,q2)'*obj.patch.positions));
            
            for m=1:obj.patch.size
                s_q = s_q+Bloch(m)*obj.s(:,:,m);
                h_q = h_q+Bloch(m)*obj.h(:,:,m);
            end
            s_q = (s_q + s_q')/2;       % Ensure Hermitian matrices
            h_q = (h_q + h_q')/2;       % Ditto
            
            % Diagonalization
            [eigenvectors,eigenvalues] = eig(h_q, s_q, 'chol', 'vector');
            
            % Eigenvalue sort and save
            eigenvalues = real(eigenvalues);
            [eigenvalues_q(:,q1,q2), I] = sort(eigenvalues);
            eigenvectors_q(:,:,q1,q2) = eigenvectors(:,I);
        end
    end
end