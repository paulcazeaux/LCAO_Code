classdef Patch
    % Gives the number and positions of some 2D lattice points on which some
    % quantity is defined.
    
    properties
        lattice
        size
        positions
        positions_on_lattice
    end
    
    methods
        function obj = plus(obj1, obj2)
            assert(obj1.lattice == obj2.lattice);
            
            obj.lattice = obj1.lattice;
            obj.positions_on_lattice = unique([obj1.positions_on_lattice';obj2.positions_on_lattice'], 'rows')';
            obj.positions = obj.lattice.L*obj.positions_on_lattice;
            obj.size = size(obj.positions,2);
        end
        
        function obj = Patch(lattice, Radius, C)
            % Computes the number and positions of 2D lattice points situated at less than
            % Radius distance of C.
            % If C is not specified, it is supposed to be located at [0;0].
            if nargin == 1
                obj.lattice = lattice;
                obj.size = 0;
                
            elseif nargin == 2
                obj.lattice = lattice;
                Nx = floor(Radius/norm(lattice.L(:,1)));
                Ny = floor(Radius/norm(lattice.L(:,2)));
                
                [M,N] = meshgrid(-Nx:Nx, -Ny:Ny);
                X = lattice.L(1,1)*M + lattice.L(1,2)*N;
                Y = lattice.L(2,1)*M + lattice.L(2,2)*N;
                R2 = X.^2 + Y.^2;
                
                I = R2 <= Radius^2;
                obj.size = sum(I(:));
                obj.positions = [X(I) Y(I)]';
                obj.positions_on_lattice = [M(I) N(I)]';
            elseif nargin >= 3
                obj.lattice = lattice;
                
                Nx2 = floor((Radius+C(1))/norm(lattice.L(:,1)));
                Nx1 = ceil((-Radius+C(1))/norm(lattice.L(:,1)));
                Ny1 = ceil((-Radius+C(2))/norm(lattice.L(:,2)));
                Ny2 = floor((Radius+C(2))/norm(lattice.L(:,2)));
                
                [M,N] = meshgrid(Nx1:Nx2, Ny1:Ny2);
                X = lattice.L(1,1)*M + lattice.L(1,2)*N;
                Y = lattice.L(2,1)*M + lattice.L(2,2)*N;
                
                R2 = (X-C(1)).^2 + (Y-C(2)).^2;
                I = (R2 <= Radius^2);
                obj.size = sum(I(:));
                obj.positions = [X(I) Y(I)]';
                obj.positions_on_lattice = [M(I) N(I)]';
            end
        end
    end
end

