function result = Normalization(obj)
    % Computes the normalization constants of the LCAO basis
    % according to formula (1.9) in Mark A. Wicks's writeup
    
    result = zeros(obj.size,1);
    for ind = 1:obj.size
        a = obj.degrees(1, ind);
        b = obj.degrees(2, ind);
        c = obj.degrees(3, ind);
        l = a+b+c;
        result(ind) = (2/pi)^.75*2^l*obj.alpha(ind)^(l/2 + .75) / sqrt(doublefact(2*a-1) * doublefact(2*b-1) * doublefact(2*c-1));
    end
end

function result = doublefact(n)
    result = prod(n:-2:1);
end