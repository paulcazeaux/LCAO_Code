function [new_basis, X] = Basis_Setup(obj, basis)
    % This function allows to change a basis to go from a description of a
    % given lattice periodic cell to a description of the supercell of a
    % superlattice.

    % new_basis is the new basis (!)
    % X is a list of coordinates of original cells 'composing' the supercell

    % Find points of the original lattice in the supercell
    j=1;
    X = zeros(2,obj.N);
    hull = [zeros(2,1) obj.S sum(obj.S,2)];
    range1 = min(hull(1,:)):max(hull(1,:));
    range2 = min(hull(2,:)):max(hull(2,:));
    for m = range1
        for n = range2
            x = obj.adjS*[m;n];
            if all(x>=0 & x<obj.N)
                X(:,j) = [m;n];
                j = j+1;
            end
        end
    end

    assert(j == obj.N+1);
    % Now we set up the new basis
    new_basis = LCAO();
    new_basis.max_degree = basis.max_degree;
    new_basis.size = obj.N*basis.size;
    new_basis.alpha = repmat(basis.alpha, obj.N, 1);
    new_basis.degrees = repmat(basis.degrees, 1, obj.N);
    new_basis.centers = zeros(3,new_basis.size);


    for j = 1:obj.N
        % First we write the new centers in the basis formed by the columns
        % of L
        new_basis.centers(1:2, basis.size*(j-1) + (1:basis.size)) = obj.lattice.L\basis.centers(1:2,:) + repmat(X(:,j), 1, basis.size);
        % In the vertical direction we keep the same values
        new_basis.centers(3, basis.size*(j-1) + (1:basis.size)) = basis.centers(3,:);
    end
    % Now we transform the basis formed by the columns of SL and we wrap
    % around in the supercell
    new_basis.centers(1:2,:) = obj.L*mod(obj.adjS*new_basis.centers(1:2,:), obj.N)/obj.N;
    new_basis.normalization = new_basis.Normalization;
end
