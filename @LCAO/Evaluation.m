function result = Evaluation(obj, Grid, R, ind)
        % computes the values of the LCAO basis based on Gaussian
        % and monomials up to degree L on
        % each point listed in Grid.
        N = obj.size;
        
        if nargin == 3
            C = obj.centers(:,1);
            GridR =  Grid - repmat(R + C, 1, size(Grid,2));
            I = obj.centers(1,:) == C(1) & obj.centers(2,:) == C(2) & obj.centers(3,:) == C(3);

            result = ones(N, size(Grid,2));
            % Unnormalized values
            
            for ind = 1:N
                d = obj.degrees(:,ind);
                if ~all(C == obj.centers(:,ind))
                    C = obj.centers(:,ind);
                    GridR =  Grid - repmat(R + C, 1, size(Grid,2));
                    I = obj.centers(1,:) == C(1) & obj.centers(2,:) == C(2) & obj.centers(3,:) == C(3);
                end
                if d(1) > 0
                    prev = find(I & obj.degrees(1,:) == d(1) - 1 & obj.degrees(2,:) == d(2) & obj.degrees(3,:) == d(3), 1);
                    result(ind,:) = GridR(1,:) .* result(prev,:);
                elseif d(2) > 0
                    prev = find(I & obj.degrees(1,:) == d(1) & obj.degrees(2,:) == d(2) - 1 & obj.degrees(3,:) == d(3), 1);
                    result(ind,:) = GridR(2,:) .* result(prev,:);
                elseif d(3) > 0
                    prev = find(I & obj.degrees(1,:) == d(1) & obj.degrees(2,:) == d(2) & obj.degrees(3,:) == d(3) - 1, 1);
                    result(ind,:) = GridR(3,:) .* result(prev,:);
                end
            end
            
            for ind = 1:N
                if ~all(C == obj.centers(:,ind))
                    C = obj.centers(:,ind);
                    GridR =  Grid - repmat(R + C, 1, size(Grid,2));
                end
                result(ind,:) = obj.normalization(ind) * exp(- obj.alpha(ind) * dot(GridR,GridR,1)) .* result(ind,:);
            end
        elseif nargin == 4
                a = obj.degrees(1,ind);
                b = obj.degrees(2,ind);
                c = obj.degrees(3,ind);
                C = obj.centers(:,ind);
                GridR =  Grid - repmat(R + C, 1, size(Grid,2));
                result = obj.normalization(ind) * exp(- obj.alpha(ind) * dot(GridR,GridR,1)) .* GridR(1,:).^a .* GridR(2,:).^b .* GridR(3,:).^c;
        else
            error('Wrong number of arguments')
        end
end