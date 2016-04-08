function result = overlap_aux(a1, b1, c1, alpha1, R1, a2, b2, c2, alpha2, R2)
    % Computes the overlap matrix element S_{ij}
    % according to formulae (2.12) in Mark A. Wicks's writeup
    
    if any([a1 b1 c1 a2 b2 c2] < 0)
        result = 0;
        return;
    end

    gamma = alpha1 + alpha2;
    P = (alpha1*R1 + alpha2*R2)/gamma;

    Ix = integral_aux(a1, a2, P(1)-R1(1), P(1)-R2(1), gamma);
    Iy = integral_aux(b1, b2, P(2)-R1(2), P(2)-R2(2), gamma);
    Iz = integral_aux(c1, c2, P(3)-R1(3), P(3)-R2(3), gamma);

    result = Ix*Iy*Iz*exp(-alpha1*alpha2/gamma*dot(R2-R1, R2-R1)); % (2.12)
end

function result = coeff_aux(k, l1, l2, pa, pb)
    % Computes the integral f_k
    % according to formulae (1.30) in Mark A. Wicks's writeup
    result = 0;
    for q = max(-k, k-2*l2):2:min(k,2*l1-k)
        i = (k+q)/2;
        j = (k-q)/2;
        result = result + nchoosek(l1,i)*nchoosek(l2,j)*pa^(l1-i)*pb^(l2-j); % (1.30)
    end
end

function result = integral_aux(l1, l2, pa, pb, gamma)
    % Computes the integral Ix
    % according to formulae (2.15) in Mark A. Wicks's writeup
    result = 0;
    for i=0:2:(l1+l2)
        result = result + coeff_aux(i,l1,l2,pa,pb)*doublefact(i-1)/(2*gamma)^(i/2)*sqrt(pi/gamma); % (2.15)
    end

end

function result = nchoosek(n,k)
    if k > n/2
        result = prod((n+1-(1:n-k))./(1:n-k));
    else
        result = prod((n+1-(1:k))./(1:k));
    end
end

function result = doublefact(n)
    result = prod(n:-2:1);
end

function result = dot(a,b) % No complex numbers here
    result =  sum(a.*b);
end