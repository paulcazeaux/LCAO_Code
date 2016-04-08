classdef Lattice
    %LATTICE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        L
        Lr
    end
    
    methods
        function obj = Lattice(L)
            if nargin >0
                obj.L = L;
                obj.Lr = 2*pi*inv(L)';
            end
        end
        
        function bool = eq(obj1, obj2)
            bool = all(obj1.L(:) == obj2.L(:));
        end
        
        function [q, K] = reciprocal_decomposition(obj, pos)
            %   for given position (pos), calculate the q such that Lr*K + q = pos,
            %   with Lr describing the reciprocal lattice of L,
            %   for some integer-valued vector K and q in the Brillouin zone.
                        
            K = floor(obj.L*pos/(2*pi)+.5);  % Brillouin Zone is centered at O
            q = pos - obj.Lr*K;
        end
        
        function result = Qgrid(obj, Nq, k)
            
            result = zeros(2, Nq, Nq);
            if nargin == 2
                for i=1:Nq
                    for j=1:Nq
                        result(:,i,j) = obj.Lr*(([i;j]-1)/Nq-.5);
                    end
                end
            elseif nargin == 3
                for i = 1:Nq
                    for j = 1:Nq
                        pos = k + obj.Lr*(([i;j]-1)/Nq-.5);
                        result(:,i,j) = obj.reciprocal_decomposition(pos);
                    end
                end
            else
                error('Wrong number of inputs!')
            end
            
        end
    end
    
end

