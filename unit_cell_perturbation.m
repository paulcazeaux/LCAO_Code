clc;clear;
% fprintf('Importing Chebfun')
% addpath('chebfun')
% fprintf('   Done\n')

fprintf('==============\nSetup...')
% Set up the LCAO basis
degree = 1;
alpha = 10;
basis = LCAO(degree, alpha);

% Then the lattice and cutoff patch
L = Lattice(eye(2));
patch = Patch(L, 2);

% Now the perturbation lattice and potential
T = Lattice(eye(2));
V_K = cell(4,1);
K_list = T.Lr*[1 -1 -1 1; 1 1 -1 -1]; 
for i = 1:4
    V_K{i} = @(z) 40*exp(-(z-1).^2);
end
Nz = 30;
V = Potential(T, K_list, V_K, Nz); % Perturbing potential Fourier transformed in the (x,y) plane 
fprintf('   Done\n')

% Compute s, h
fprintf('==============\nComputing tight-binding matrices...')
Atomic_Cutoff = 2; Z = 6; % Atomic charge
[s,h] = sh(basis, patch, Atomic_Cutoff, Z);
fprintf('   Done\n')

% Compute D0
fprintf('==============\nComputing D0...')
gap_position = 1; Nq = 30;
D0 = PeriodicDensityMatrix(L, basis, patch, h, s, gap_position, Nq);
fprintf('   Done\n')

% Compute D1
fprintf('==============\nComputing D1...')
D1 = PerturbationDensityMatrix(D0, patch, V, Nq);
fprintf('   Done\n')

% Now the perturbation error tests
fprintf('==============\nComputing the perturbation error for different potential scalings...')
w = squeeze(sum(basis.W(V, patch),3));
N = 10;
t = 1000.^((0:N-1)/(N-1)-1);
err0 = zeros(1,N); err1 = zeros(1,N);
B0 = zeros(basis.size,basis.size,patch.size,patch.size);
B1 = zeros(basis.size,basis.size,patch.size,patch.size);
for r1 = 7:7
    for r2 = 1:patch.size
        B0(:,:,r1,r2) = D0.Block(patch.positions(:,r1), patch.positions(:,r2));
        B1(:,:,r1,r2) = D1.Block(patch.positions(:,r1), patch.positions(:,r2));
    end
end

for i = 1:N
    D = PeriodicDensityMatrix(L, basis, patch, h+t(i)*w, s, gap_position, Nq);
    for r1 = 7:7
        for r2 = 1:patch.size
            B = D.Block(patch.positions(:,r1), patch.positions(:,r2));
            err1(i) = max(err1(i), norm(B0(:,:,r1,r2) + t(i)*B1(:,:,r1,r2) - B));
            err0(i) = max(err0(i), norm(B0(:,:,r1,r2)-B));
        end
    end
end
fprintf('   Done\n')
figure(1); clf;
loglog(t, err0, 'b-'); hold on;
loglog(t, err1, 'r-');