function result = W(obj, V, patch)
% calculates (40) from Perturbation Notes
% N is the size of the Gauss-Hermite discretization
% K is a vector of the reciprocal perturbing lattice, expressed in the
%   basis given by the columns of Tr
% V is a vector giving the values of the function V(K,z) for the partial Fourier 
% transform in 2D K coefficient as a function of z

% Uses "Fourier transforms of Gaussian orbital products," by 
% G. S. Chanlder and M. A. Spackman
    result = zeros(obj.size, obj.size, size(V.K_list,2), patch.size);
    Tr = V.lattice.Lr;
    
    for ind1 = 1:obj.size
        for ind2 = 1:obj.size
            a1 = obj.degrees(1,ind1);
            b1 = obj.degrees(2,ind1);
            c1 = obj.degrees(3,ind1);
            alpha1 = obj.alpha(ind1);
            z1 = obj.centers(3,ind1);
            
            a2 = obj.degrees(1,ind2);
            b2 = obj.degrees(2,ind2);
            c2 = obj.degrees(3,ind2);
            alpha2 = obj.alpha(ind2);
            z2 = obj.centers(3,ind2);
            gamma = alpha1+alpha2;
            
            N = obj.normalization(ind1)*obj.normalization(ind2)/sqrt(pi*abs(det(Tr)));

            % Change of variables
            z = V.Quad.z/sqrt(gamma) + (alpha1*z1 + alpha2*z2)/gamma;

            % w = w / sqrt(gamma)*exp(-alpha1*alpha2/gamma*(z1-z2)^2);
            % Note: the change of variables cancels out the normalization constant
            % arising from the use of LCAO_fourier_aux (3-D integral), see LCAO_Fourier_2D
            % below calculates the z-component, amputed of the exponential weight
            % term
            
            R2 = obj.centers(:,ind2);
            for k = 1:size(V.K_list,2)
                K = [V.K_list(:,k);0];
                Iz = dot((z - z1).^c1 .* (z - z2).^c2 .* V.Evaluate(z,k), V.Quad.w);
                for r = 1:patch.size
                    R1 = [patch.positions(:,r);0] + obj.centers(:,ind1);
                    Ixy = LCAO.fourier_aux(a1, b1, 0, alpha1, R1, a2, b2, 0, alpha2, R2, K);
                    result(ind1, ind2, k,r) = N*Ixy*Iz ;
                end
            end
        end
    end
end