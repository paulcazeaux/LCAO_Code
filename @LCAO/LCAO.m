classdef LCAO < Basis
    %LCAO Basis class - implementation of LCAO data structure and routines
    
    properties
        size;
        centers;
        normalization;
        max_degree;
        degrees;
        alpha;
    end
    
    methods (Static)
        result = overlap_aux(a1, b1, c1, alpha1, R1, a2, b2, c2, alpha2, R2);
        result = fourier_aux(a1, b1, c1, alpha1, R1, a2, b2, c2, alpha2, R2, K);
    end
    
    methods
        result = Normalization(obj);
        result = Evaluation(obj, Grid, R, ind);
        result = Overlap(obj, R1, R2);
        result = Kinetic(obj, R1, R2);
        result = Nuclear(obj, R1, R2, Z, C);
        result = Fourier(obj, R1, R2, K)
        result = Fourier_2D(obj, R1, R2, K);
        
        value = W(obj, V, patch)
    end
    
    methods % Constructor
        function obj = LCAO(max_degree, alpha, centers)
            if nargin >0
                obj.max_degree = max_degree;
                local_size = nchoosek(max_degree+3, max_degree);
            end
            if nargin == 1
                % Default value for alpha : 1
                obj.size = local_size;
                obj.alpha = ones(obj.size, 1);
                
                obj.degrees = zeros(3,obj.size);
                ind = 1;
                for l = 1:max_degree
                    for c=0:l
                        for b=0:(l-c)
                            ind = ind+1;
                            a = l-b-c;
                            obj.degrees(:,ind) = [a;b;c];
                        end
                    end
                end
                
                % Default orbital center at the origin
                obj.centers = zeros(3, obj.size);
                obj.normalization = obj.Normalization;
                
            elseif nargin == 2
                obj.size = local_size;
                if length(alpha) == obj.size
                    obj.alpha = alpha;
                elseif isscalar(alpha)
                    obj.alpha = alpha*ones(obj.size,1);
                else
                    error('alpha input has the wrong size or format');
                end
                
                obj.degrees = zeros(3, local_size);
                ind = 1;
                for l = 1:obj.max_degree
                    for c=0:l
                        for b=0:(l-c)
                            ind = ind+1;
                            a = l-b-c;
                            obj.degrees(:,ind) = [a;b;c];
                        end
                    end
                end
                
                % Default orbital center at the origin
                obj.centers = zeros(3, obj.size);
                obj.normalization = obj.Normalization;
                
            elseif nargin >= 3
                % Set up a basis replicating a LCAO basis with multiple centers
                assert(size(centers, 1) == 3); %Check that centers has the right format
                
                obj.size = local_size * size(centers, 2);
                if size(alpha,1) == local_size
                    obj.alpha = repmat(alpha, size(centers,2),1);
                elseif isscalar(alpha)
                    obj.alpha = alpha*ones(obj.size,1);
                else
                    error('alpha input has the wrong size or format');
                end
                
                centers_number = size(centers,2);
                
                obj.degrees = zeros(3, local_size);
                ind = 1;
                for l = 1:obj.max_degree
                    for c=0:l
                        for b=0:(l-c)
                            ind = ind+1;
                            a = l-b-c;
                            obj.degrees(:,ind) = [a;b;c];
                        end
                    end
                end
                obj.degrees = repmat(obj.degrees, 1, centers_number);
                
                obj.centers = repmat(centers, 1, local_size);
                I = reshape(reshape(1:obj.size, centers_number, local_size)', 1, obj.size);
                obj.centers = obj.centers(:,I);
                obj.normalization = obj.Normalization;
            end
            
        end
    end
end

