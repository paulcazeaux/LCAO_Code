function result = fourier_aux(a1, b1, c1, alpha1, R1, a2, b2, c2, alpha2, R2, K)
    % Computes the Fourier matrix element
    % according to formulae (7) in Fourier transforms of Gaussian orbital
    % products, by G. S. CHANDLER and M. A. SPACKMAN, 1977
    
    if any([a1 b1 c1 a2 b2 c2] < 0)
        result = 0;
        return;
    end

    gamma = alpha1 + alpha2;
    P = (alpha1*R1 + alpha2*R2)/gamma;

    Ix = integral_aux(a1, a2, P(1)-R1(1), P(1)-R2(1), gamma, K(1));
    Iy = integral_aux(b1, b2, P(2)-R1(2), P(2)-R2(2), gamma, K(2));
    Iz = integral_aux(c1, c2, P(3)-R1(3), P(3)-R2(3), gamma, K(3));

    result = Ix*Iy*Iz*(pi/gamma)^1.5*exp(-alpha1*alpha2/gamma*dot(R2-R1, R2-R1) + 1i*dot(K,P) - dot(K,K)/(4*gamma)); % (2.12)
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

function result = integral_aux(l1, l2, pa, pb, gamma, k)
    % Computes the integral Ix
    % according to formulae (2.15) in Mark A. Wicks's writeup
    result = 0;
    for i=0:(l1+l2)
        result = result + coeff_aux(i,l1,l2,pa,pb)*HermiteH(i, k/(2*sqrt(gamma)))/(-2*1i*sqrt(gamma))^i; % (2.15)
    end

end

function result = nchoosek(n,k)
    if k > n/2
        result = prod((n+1-(1:n-k))./(1:n-k));
    else
        result = prod((n+1-(1:k))./(1:k));
    end
end


function result = dot(a,b) % No complex numbers here
    result =  sum(a.*b);
end

function result = HermiteH(i,x)
    assert(i>=0)
    if i==0
        result = ones(size(x)); return;
    elseif i==1
        result = 2*x; return;
    else
        n = length(x);
        X = reshape(2*x, n, 1);
        H = zeros(n,i);
        H(:,1) = X;
        H(:,2) = X.^2 - 2;
        for j=3:i
            H(:,j) = X.*H(:,j-1) - 2*(j-1)*H(:,j-2);
        end
    end
    result = reshape(H(:,i), size(x));
end