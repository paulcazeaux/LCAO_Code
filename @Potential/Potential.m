classdef Potential
    %PERTURBINGPOTENTIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lattice
        K_list
        Quad
        FourierCoefficients
    end
    
    methods (Static)
        [X, W] = hermquad(N);
    end
    
    methods
        function obj = Potential(lattice, K_list, FourierCoefficients, Nz)
            obj.lattice = lattice;
            obj.K_list = K_list;
            [obj.Quad.z, obj.Quad.w] = obj.hermquad(Nz);
            obj.FourierCoefficients = FourierCoefficients;
        end
        
        function result = Evaluate(obj,z,k)
            result = obj.FourierCoefficients{k}(z);
        end
    end
end

