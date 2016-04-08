function result = Overlap(obj, R1, R2)
    % Normalizes the overlap matrix element S(R1, R2)
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

            result(ind1, ind2) = N*LCAO.overlap_aux(a1, b1, c1, alpha1, r1, a2, b2, c2, alpha2, r2);
        end
    end
end