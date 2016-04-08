function result = Kinetic(obj, R1, R2)
    % Computes the Kinetic energy matrix element T(R1,R2)
    % according to formulae (3.1) and (3.13) in Mark A. Wicks's writeup
    
    result = zeros(obj.size);
    for ind1 = 1:obj.size
        for ind2 = 1:obj.size
            r1 = R1 + obj.centers(:,ind1);
            a1 = obj.degrees(1,ind1);
            b1 = obj.degrees(2,ind1);
            c1 = obj.degrees(3,ind1);
            alpha1 = obj.alpha(ind1);

            r2 = R2 + obj.centers(:,ind2);
            a2 = obj.degrees(1,ind2);
            b2 = obj.degrees(2,ind2);
            c2 = obj.degrees(3,ind2);
            alpha2 = obj.alpha(ind2);

            N = obj.normalization(ind1)*obj.normalization(ind2); 

            Ix = .5*a1*a2              *LCAO.overlap_aux(a1-1, b1, c1, alpha1, r1, a2-1, b2, c2, alpha2, r2)...  (3.13)
                    + 2*alpha1*alpha2  *LCAO.overlap_aux(a1+1, b1, c1, alpha1, r1, a2+1, b2, c2, alpha2, r2) ...
                    - alpha1*a2        *LCAO.overlap_aux(a1+1, b1, c1, alpha1, r1, a2-1, b2, c2, alpha2, r2) ...  
                    - alpha2*a1        *LCAO.overlap_aux(a1-1, b1, c1, alpha1, r1, a2+1, b2, c2, alpha2, r2);

            Iy = .5*b1*b2              *LCAO.overlap_aux(a1, b1-1, c1, alpha1, r1, a2, b2-1, c2, alpha2, r2)...  (3.13)
                    + 2*alpha1*alpha2  *LCAO.overlap_aux(a1, b1+1, c1, alpha1, r1, a2, b2+1, c2, alpha2, r2) ... 
                    - alpha1*b2        *LCAO.overlap_aux(a1, b1+1, c1, alpha1, r1, a2, b2-1, c2, alpha2, r2) ...
                    - alpha2*b1        *LCAO.overlap_aux(a1, b1-1, c1, alpha1, r1, a2, b2+1, c2, alpha2, r2);

            Iz = .5*c1*c2              *LCAO.overlap_aux(a1, b1, c1-1, alpha1, r1, a2, b2, c2-1, alpha2, r2)...  (3.13)
                    + 2*alpha1*alpha2  *LCAO.overlap_aux(a1, b1, c1+1, alpha1, r1, a2, b2, c2+1, alpha2, r2) ...
                    - alpha1*c2        *LCAO.overlap_aux(a1, b1, c1+1, alpha1, r1, a2, b2, c2-1, alpha2, r2) ...
                    - alpha2*c1        *LCAO.overlap_aux(a1, b1, c1-1, alpha1, r1, a2, b2, c2+1, alpha2, r2);

            result(ind1,ind2) = N*(Ix + Iy + Iz);   % (3.1)
        end
    end
end