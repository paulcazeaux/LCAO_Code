L = eye(2); % Layer lattice
Lr = 2*pi*L^(-1);        % Reciprocal lattice

Rot = sqrt(2)*eye(2);
T = Rot*L; % Perturbation lattice
Tr = 2*pi*T^-1;
z0 = 1.5;
V = @(K,z) .25*all(abs(K) == 1)*exp(-(z-z0).^2); % Perturbing potential Fourier transformed in the (x,y) plane 

% set up the LCAO basis
degree = 1;
alpha = 10;

basis = LCAO.Setup(degree, alpha);
Info.N = 30;
[Info.Quad.z, Info.Quad.w] = PerturbationPotential.hermquad(Info.N);

Cutoff = 1; %  nearest neighbor distance cutoff
patch = geometry.Patch(L,Cutoff);

% gather relevant K_tilde's

K_tilde = [1 -1 -1 1; 1 1 -1 -1];

Nq = 10; % number of q points in discretization

% build the q' values that correspond to K_tilde's and the q values.
Q = zeros(2, Nq, Nq);
Q_p = zeros([size(K_tilde) Nq Nq]);
K_p = zeros([size(K_tilde) Nq Nq]);

for i = 1:Nq
    for j = 1:Nq
        Q(:,i,j) = Tr*([-1;-1] + 2/Nq*[i-1;j-1]);
        [Q_p(:,:,i,j), K_p(:,:,i,j)] = geometry.reciprocal_decomposition(Lr, Tr*K_tilde(1:2,:) + repmat(Q(:,i,j),1,size(K_tilde,2)));
    end
end

% build coefficients for calculating W floquet bloch terms.
W_co = perturbation.W_coeff(Info.Quad, basis, patch, Tr, V, K_tilde);

W_K = zeros(basis.size,basis.size,size(K_tilde,2),Nq,Nq);

for q1 = 1:Nq
    for q2 = 1:Nq
        W_K(:,:,:,q1,q2) = perturbation.W_floquet(basis, patch, W_co, K_tilde, Q_p(:,:,q1,q2));
    end
end

