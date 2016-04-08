function result = Supercell_Block(obj, D, R1, R2)
    [new_basis, X] = obj.Basis_Setup(D.basis);
    result = zeros(new_basis.size, new_basis.size);
    for r1 = 1:size(X,2)
        for r2 = 1:size(X,2)
            r1_range = D.basis.size*(r1-1) + (1:D.basis.size);
            r2_range = D.basis.size*(r2-1) + (1:D.basis.size);
            result(r1_range,r2_range) = D.Block(R1+obj.lattice.L*X(:,r1), R2+obj.lattice.L*X(:,r2));
        end
    end
end