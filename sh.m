function [s,h] = sh(basis,patch,Atomic_Cutoff,Z)
    %%%%%%%%%%%%%
    % build h & s
    h = zeros(basis.size,basis.size,patch.size); % (26) in DensityPerturbation
    s = zeros(basis.size,basis.size,patch.size);
    
    O = [0;0;0];
    patch_O = Patch(patch.lattice, Atomic_Cutoff, O);
    
    for m = 1:patch.size
        R = [patch.positions(:,m) ;0];
        s(:,:,m) = basis.Overlap(O, R);
        h(:,:,m) = basis.Kinetic(O, R);
        
        patch_atoms = patch_O + Patch(patch.lattice,  Atomic_Cutoff, R);
        for c = 1:patch_atoms.size
            C = [patch_atoms.positions(:,c) ;0];
            h(:,:, m) = h(:,:, m) + basis.Nuclear(O, R, Z, C);
        end
    end
end