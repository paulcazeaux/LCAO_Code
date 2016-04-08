
%%%%%%%%%%%%%%%%%%%
% Overlap integrals

% set up the LCAO basis
degree = 3;
alpha = 4;
basis = LCAO.Setup(degree, alpha);
ind1 = 2; R1 = [0;0;0];
ind2 = 2; R2 = [1;0;0];

f = @(x,y,z) LCAO.Evaluation([x;y;z], basis, R1, ind1) .*LCAO.Evaluation([x;y;z], basis, R2, ind2);

disp('Test Overlap integrals:')

LCAO.Overlap(basis, ind1, R1, ind2, R2) - integral3(f, -Inf, Inf, -Inf, Inf, -Inf, Inf)
%%

%%%%%%%%%%%%%%%%%%%
% Fourier integrals

% set up the LCAO basis
degree = 1;
alpha = 4;
basis = LCAO(degree, alpha);
R1 = [0;0;0];
R2 = [1;0;0];
K = [1;0;0];


disp('Test Fourier integrals:')
Ref = zeros(basis.size);
Test = basis.Fourier(R1, R2, K);
for ind1 = 1:basis.size;
    for ind2 = 1:basis.size; 
        f = @(x,y,z) basis.Evaluation([x;y;z], R1, ind1) .* basis.Evaluation([x;y;z], R2, ind2) .* exp(-1i*dot(repmat(K, 1, length(x)), [x;y;z], 1));
        Ref(ind1, ind2) = integral3(f, -Inf, Inf, -Inf, Inf, -Inf, Inf);
    end
end

norm(Test-Ref)
%%
disp('Test 2D Fourier integrals:')

Test = basis.Fourier_2D(R1,R2, K);
Ref = zeros(basis.size);
for ind1 = 1:basis.size;
    for ind2 = 1:basis.size; 
        f = @(x,y) basis.Evaluation([x;y; ones(1, length(x))], R1, ind1) .* basis.Evaluation([x;y; ones(1, length(x))], R2, ind2) .* exp((basis.alpha(ind1)+basis.alpha(ind2))+1i*dot(repmat(K,1,length(x)), [x;y;zeros(1, length(x))], 1));
        Ref(ind1, ind2) = integral2(f, -Inf, Inf, -Inf, Inf);
    end
end
norm(Test-Ref)


%%
disp('Test LCAO.W integrals:')

% set up geometry
L = [1,0; 0,1]; % Plane layer lattice; columns are lattice vectors
Rot = eye(2);
T = Rot*L; % Perturbation lattice
z0 = 1.5;
Tr = 2*pi*T^-1;

% set up the LCAO basis
degree = 1;
alpha = 4;
basis = LCAO.Setup(degree, alpha);
R = [1;0;0];
K = [1;-1];
N = 30;
[Quad.z, Quad.w] = PerturbationPotential.hermquad(Info.N);

V = @(K,z) exp(-(z-z0).^2);

Test = zeros(basis.size);
Ref = zeros(basis.size);
for ind1 = 1:basis.size;
    for ind2 = 1:basis.size; 
        Test(ind1, ind2) = LCAO.W(Quad, basis, Tr, V, ind1, ind2, R, K);
        f = @(x,y,z) LCAO.Evaluation([x;y;z], basis, [0;0;0], ind1) .*LCAO.Evaluation([x;y;z], basis, R, ind2) .* V(-K,z) .* exp(-1i*dot(repmat(Tr*K,1,length(x)), [x;y], 1));
        Ref(ind1, ind2) = integral3(f, -Inf, Inf, -Inf, Inf, -Inf, Inf)/sqrt(abs(det(Tr)));
    end
end

norm(Test-Ref)

%%
%%%%%%%%%%%%%%%%%%%
% Nuclear integrals
% 
% set up the LCAO basis
degree = 1;
alpha = 4;
basis = LCAO.Setup(degree, alpha);

R1 = [1;0;0];
R2 = [0;1;0];
Z = 1; C = [0;0;0];
Test = zeros(basis.size);
Ref = zeros(basis.size);
for ind1 = 1:basis.size;
    for ind2 = 1:basis.size; 
        Test(ind1, ind2) = LCAO.Nuclear2(basis, ind1, R1, ind2, R2, Z, C);
        f = @(x,y,z) LCAO.Evaluation([x;y;z], basis, R1, ind1) .*LCAO.Evaluation([x;y;z], basis, R2, ind2) .* (-Z./sqrt(dot([x;y;z] - repmat(C,1,length(x)), [x;y;z] - repmat(C,1,length(x)), 1)));
        Ref(ind1, ind2) = integral3(f, -Inf, Inf, -Inf, Inf, -Inf, Inf);
        
    end
end

norm(Test-Ref)
