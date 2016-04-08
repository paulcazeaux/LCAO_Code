function result = Nuclear(obj, R1, R2, Z, C)
    % Computes the Nuclear attraction energy matrix element 
    % V(R1, R2)(C)
    % according to formula (56) in T. Petersson, B. Hessing, Detailed
    % derivation of Gaussian orbital based matrix elements in electron
    % structure calculations.
    
    result = zeros(obj.size);
    for ind1 = 1:obj.size
        for ind2 = 1:obj.size
            r1 = R1 + obj.centers(:,ind1);
            a1 = obj.degrees(1,ind1);
            b1 = obj.degrees(2,ind1);
            c1 = obj.degrees(3,ind1);
            l1 = a1 + b1 + c1;
            alpha1 = obj.alpha(ind1);

            r2 = R2 + obj.centers(:,ind2);
            a2 = obj.degrees(1,ind2);
            b2 = obj.degrees(2,ind2);
            c2 = obj.degrees(3,ind2);
            l2 = a2 + b2 + c2;
            alpha2 = obj.alpha(ind2);

            gamma = alpha1 + alpha2;
            P = (alpha1*r1 + alpha2*r2)/gamma;

            % Precomputation
            fact = factorial(0:max([a1+a2 b1+b2 c1+c2]));


            K = -2*Z*pi*obj.normalization(ind1)*obj.normalization(ind2)/gamma...
                    *(-1)^(l1+l2)*fact(1+a1)*fact(1+a2)*fact(1+b1)*fact(1+b2)*fact(1+c1)*fact(1+c2)...
                    *exp(-alpha1*alpha2/gamma*dot(r2-r1, r2-r1));

            F = Fnu(l1+l2, gamma*dot(P-C, P-C));
            Ix = integrals_aux(a1, a2, r1(1)-r2(1), P(1)-C(1), alpha1, alpha2, fact);
            Iy = integrals_aux(b1, b2, r1(2)-r2(2), P(2)-C(2), alpha1, alpha2, fact);
            Iz = integrals_aux(c1, c2, r1(3)-r2(3), P(3)-C(3), alpha1, alpha2, fact);

            for nx = 0:a1+a2
                for ny=0:b1+b2
                    for nz = 0:c1+c2
                        result(ind1,ind2)  = result(ind1,ind2) + F(nx+ny+nz+1)*Ix(nx+1)*Iy(ny+1)*Iz(nz+1);
                    end
                end
            end
            result(ind1,ind2) = K*result(ind1,ind2);   % (3.1)
        end
    end
end


function result = Fnu(n, u)
    % Returns the integral F_nu(u) (29) in Petersson & Hessing
    % for all 0 <= j <= n
    if u == 0
        result = 1./(2*(0:n)+1)';
        return;
    end
    
    result = zeros(n+1,1);
    result(end) = prod(2*n-1:-2:1) ...
            * ( ...
                sqrt(pi/2)/(2*u)^(n+.5)*erf(sqrt(u)) ...
                - exp(-u)*sum(1./(flip(cumprod(flip(2*n-1:-2:1))).*(2*u).^(1:n))) ...
            );
    for j = n:-1:1
        result(j) = (2*u*result(j+1) + exp(-u))/(2*j-1);
    end
    
end

function result = integrals_aux(l1, l2, ab, pc, alpha1, alpha2, fact)
    % We compute some factors in formula (56) in T. Petersson, B. Hessing
    % representing one-dimensional integrals, for nu_x = mu-u fixed
    
    if nargin == 6
        fact = factorial(0:l1+l2);
    end

    % First we compute some useful factors
    % This takes care of the sum over r
    gamma = alpha1+alpha2;
    c = zeros(l1+1, l2+1);
    if abs(ab)<1e-10
        % For ab = 0, the formula simplifies since o1+o2 = 2*r is the only
        % contribution
        for r = 0:floor((l1+l2)/2)
            for q = max(-r,r-l2):min(r,l1-r)
                o1 = r+q;
                o2 = r-q;
                c(o1+1, o2+1) =  (-alpha1)^o2 * alpha2^o1  ...
                        * fact(1+o1+o2)/(fact(1+o1)*fact(1+o2)*fact(1+r)*(-4*alpha1*alpha2*gamma)^r);
            end
        end
    else
        for o1 = 0:l1
            for o2 = 0:l2
                for r=0:floor((o1+o2)/2)
                    c(o1+1, o2+1) = c(o1+1, o2+1)+1/((-4*alpha1*alpha2*ab^2/gamma)^r*fact(1+r)*fact(1+o1+o2-2*r));
                end
                c(o1+1, o2+1) =  (-alpha1)^o2 * alpha2^o1 * (ab/gamma)^(o1+o2) * fact(1+o1+o2)/(fact(1+o1)*fact(1+o2)) * c(o1+1, o2+1);
            end
        end
    end
    
    % Now some factors depending only on i1 and i2
    fi1 = zeros(floor(l1/2)+1,1);
    fi2 = zeros(floor(l2/2)+1,1);
    for i1=0:l1/2
        fi1(i1+1) = fact(1+i1)*(4*alpha1)^i1;
    end
    for i2=0:l2/2
        fi2(i2+1) = fact(1+i2)*(4*alpha2)^i2;
    end

    % Now we compute fixed factors for fixed mu or u
    a = zeros(l1+l2+1,1);
    b = zeros(l1+l2+1,1);
    for mu = 0:(l1+l2)
        for k=0:min(l1,mu)
            A = 0;
            for i1=0:floor((l1-k)/2)
                for i2=0:floor((l2+k-mu)/2)
                    o1 = l1-2*i1-k;
                    o2 = l2-2*i2+k-mu;
                    A = A + c(o1+1,o2+1)/(fi1(i1+1)*fi2(i2+1));
                end
            end
            a(mu+1) = a(mu+1)+A/(fact(1+k)*fact(1+mu-k));
        end
        a(mu+1) = a(mu+1)*fact(1+mu);
    end
    for mu=0:(l1+l2)
        b(mu+1) = pc^mu/fact(1+mu);
    end
    
    g = zeros(floor((l1+l2)/2)+1,1);
    for u = 0:floor((l1+l2)/2)
        g(u+1) = fact(1+u)*(-4*gamma)^u;
    end
    g = 1./g;
    
    % Finally we compute a factor for given nu = mu-u
    result = zeros(l1+l2+1, 1);
    for nu = 0:(l1+l2)
        for mu = nu:min(2*nu,l1+l2)
            u = mu - nu;
            result(nu+1) = result(nu+1) + a(mu+1)*b(mu-2*u+1)*g(u+1);
        end
    end
end

function result = dot(a,b) % No complex numbers here
    result =  sum(a.*b);
end