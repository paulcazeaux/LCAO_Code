function result = Fourier_2D(obj, R1, R2, K)
    % Normalizes the (x,y) plane overlap of the (x,y) plane portions of
    % two orbitals with a Fourier phase term exp^(+i K.x)
    
    result = zeros(obj.size);
    for ind1 = 1:obj.size
        for ind2 = 1:obj.size
            
            r1 = R1 + obj.centers(:,ind1);
            z1 = obj.centers(3,ind1);
            a1 = obj.degrees(1,ind1);
            b1 = obj.degrees(2,ind1);
            alpha1 = obj.alpha(ind1);
            
            r2 = R2 + obj.centers(:,ind2);
            z2 = obj.centers(3,ind2);
            a2 = obj.degrees(1,ind2);
            b2 = obj.degrees(2,ind2);
            alpha2 = obj.alpha(ind2);
            gamma = alpha1+alpha2;
            
            N = sqrt(gamma/pi)*exp(alpha1*alpha2/gamma*(z1-z2)^2)*obj.normalization(ind1)*obj.normalization(ind2);
            
            result(ind1,ind2) = N*LCAO.fourier_aux(a1, b1, 0, alpha1, r1, a2, b2, 0, alpha2, r2, [K(1:2);0]);
        end
    end
end