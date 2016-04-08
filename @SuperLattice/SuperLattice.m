classdef SuperLattice < Lattice
    %SUPERLATTICE Summary of this class goes here
    %   Detailed explanation goes here
    properties
        lattice % Underlying original lattice
        S
        N
        adjS
    end
    
    methods
        function obj = SuperLattice(lattice, S)
            obj@Lattice(lattice.L*S);
            obj.lattice = lattice;
            obj.S = S;
            obj.N = round(abs(det(obj.S)));
            obj.adjS = round(obj.N*inv(obj.S));
        end
        
        [new_basis, X] = Basis_Setup(obj, basis)
        [new_op, new_patch] = Change(obj, basis, patch, op)
        [new_W, new_patch] = Change_W(obj, basis, patch, K_list, W)
        result = Supercell_Block(obj, D, R1, R2)
        
    end
    
end

