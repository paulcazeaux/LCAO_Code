classdef (Abstract) Basis
    %BASIS Virtual class for Basis construction
    
    properties (Abstract)
        size;
        centers;
    end
    
    methods (Abstract)
        result = Evaluation(obj, Grid, R, ind);
        result = Overlap(obj, R1, R2);
        result = W(obj, V, patch);
    end
end

