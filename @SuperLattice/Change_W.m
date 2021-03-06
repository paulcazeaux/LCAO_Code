function [new_W, new_patch] = Change_W(obj, basis, patch, K_list, W)
    % The change of basis for W_K(R) requires some tweaks, namely
    % attention to the phase changes
    [new_basis, X] = obj.Basis_Setup(basis);

    % First, we identify the new patch of supercells
    new_patch = Patch(obj.L);
    positions_on_lattice = [];
    for r = 1:size(X,2)
        x = obj.adjS*(patch.positions_on_lattice + repmat(X(:,r), 1, patch.size));
        x = round((x-mod(x, obj.N))/obj.N);
        positions_on_lattice = unique([positions_on_lattice x]', 'rows')';
    end
    new_patch.size = size(positions_on_lattice, 2);
    new_patch.positions_on_lattice = positions_on_lattice;
    new_patch.positions = obj.L*positions_on_lattice;

    % Now, we set up the new operator
    new_W = zeros(new_basis.size, new_basis.size, size(K_list,2), new_patch.size);
    for r2 = 1:size(X,2)
        for i=1:patch.size
            x = patch.positions_on_lattice(:,i)+X(:,r2);
            xR = round((obj.adjS*x-mod(obj.adjS*x, obj.N))/obj.N);
            xr = round(x - obj.S*xR);
            [~,R] = ismember(xR', new_patch.positions_on_lattice', 'rows');
            [~,r1] = ismember(xr', X', 'rows');
            r1_range = basis.size*(r1-1) + (1:basis.size);
            r2_range = basis.size*(r2-1) + (1:basis.size);
            for k=1:size(K_list,2)
                new_W(r1_range,r2_range,k,R) = W(:,:,k,i)*exp(1i*K_list(:,k)'*obj.lattice.L*X(:,r2));
            end
        end
    end

end