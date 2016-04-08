clc;clear;

fprintf('==============\nLayer geometry...')
% Set up the LCAO basis
degree = 1;
alpha = 10;
basis = LCAO(degree, alpha, [0 1/2; 0 sqrt(3)/6; 0 0]);

% Then the lattice and cutoff patch
L = Lattice([1 cos(pi/3); 0 sin(pi/3)]);
patch = Patch(L, 2);
fprintf('   OK\n')

% Compute s, h
fprintf('==============\nComputing tight-binding matrices...')
Atomic_Cutoff = 2; Z = 6; % Atomic charge
[s,h] = sh(basis, patch, Atomic_Cutoff, Z);
fprintf('   OK\n')

% Compute D0
fprintf('==============\nComputing D0...')
gap_position = 4; Nq = 30;
D0 = PeriodicDensityMatrix(L, basis, patch, h, s, gap_position, Nq);
fprintf('   OK\n')

% %% Plot the band structure
% figure(1); clf;
% surf(squeeze(D0.Qgrid(1,:,:)), squeeze(D0.Qgrid(2,:,:)), squeeze(D0.eigenvalues_q(1,:,:)))
% hold on
% for i=2:basis.size
%     surf(squeeze(D0.Qgrid(1,:,:)), squeeze(D0.Qgrid(2,:,:)), squeeze(D0.eigenvalues_q(i,:,:)))
% end
% 


gap = min(reshape(D0.eigenvalues_q(gap_position+1,:,:),[],1)) - max(reshape(D0.eigenvalues_q(1:gap_position,:,:),[],1));
fprintf('Test band gap: gap = %f\n--------------\n\n', gap)

% Now the perturbation lattice and potential
fprintf('--------------\nSet up perturbation...')
p=5; q=16; n=19;
T = Lattice([2*p+q p-q; sqrt(3)*q sqrt(3)*(p+q)]/(2*n));
V_K = cell(4,1);
K_list = T.Lr*[1 -1 -1 1; 1 1 -1 -1]; 
for i = 1:4
    V_K{i} = @(z) 40*exp(-(z-1).^2);
end
Nz = 30;
V = Potential(T, K_list, V_K, Nz); % Perturbing potential Fourier transformed in the (x,y) plane 


% Here we construct the superlattice
S = SuperLattice(L, [p -q; q p+q]);
fprintf('   OK\n')

% Compute D1
fprintf('==============\nComputing D1...')
D1 = PerturbationDensityMatrix(D0, patch, V, Nq);
fprintf('   OK\n')


fprintf('==============\nSet up Supercell matrices...')
% Construct superlattice matrices
[super_basis,X] = S.Basis_Setup(basis);
[super_s, super_patch] = S.Change(basis, patch, s);
super_h = S.Change(basis, patch, h);
% Assuming that V is LS-Periodic
super_w = squeeze(sum(S.Change_W(basis, patch, V.K_list, basis.W(V, patch)),3));
%super_w_2 = squeeze(sum(super_basis.W(V, super_patch),3));
%norm(super_w(:) - super_w_2(:), 'inf')
fprintf('   OK\n--------------\n\n')


% Now the perturbation error tests
N = 10;
t = 100.^((0:N-1)/(N-1)-1);
O = [0;0];

fprintf('--------------\nComputing the zero and first-order perturbation blocks...')
B0 = S.Supercell_Block(D0, O, O);
B1 = S.Supercell_Block(D1, O, O);
fprintf('   OK\n')
%%
fprintf('==============\nComputing the perturbation error for different potential scalings...')
err0 = zeros(1,N); err1 = zeros(1,N);
parfor i = 1:N
    D = PeriodicDensityMatrix(S, super_basis, super_patch, super_h+t(i)*super_w, super_s, S.N*gap_position, 5);

    B{i} = D.Block(O,O);
    err1(i) = max(err1(i), norm(B{i} - B0 - t(i)*B1));
    err0(i) = max(err0(i), norm(B{i} - B0));
end
fprintf('   OK\n--------------\n\n')

fprintf('--------------\nSaving and plotting...')
save(sprintf('%ix%i_Workspace.mat',n ,n));

figure(1); clf;
loglog(t, err0, 'b-'); hold on;
loglog(t, err1, 'r-');

% %% Plot the resulting density matrix
% figure(3); clf;
% scatter3(super_basis.centers(1,:)', super_basis.centers(2,:)', abs(B1(:,5))/norm(B1, 'inf'), 'filled')
% ax = gca;
% ax.DataAspectRatio = [1 1 .01];