%% Lattice parameters and LCAO basis construction
%clc; close all; clear all;

%%%%%%%%%%%%%%%%%%%%%
% Set up the geometry
L = eye(2);
S = 3*L;

%%%%%%%%%%%%%%%%%%%%%%%
% set up the LCAO basis
degree = 2;
alpha = 4;
basis = LCAO.Setup(degree, alpha);


N = 7;
D = rand(basis.size, basis.size, N, N, N, N);

disp('Testing basis change...')
[new_basis, R, P] = supercell.Basis_Setup(basis, L, S);
new_D = supercell.Matrix_Change_To_Supercell(basis, L, S, D);
Test = supercell.Matrix_Change_From_Supercell(basis, L, S, new_D) - D;

if (norm(Test(:)) < 1e-10)
    disp('OK')
else
    disp('NOK! Watch out')
end