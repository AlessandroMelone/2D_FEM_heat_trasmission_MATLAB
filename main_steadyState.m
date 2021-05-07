clear, clc, close all
addpath('FEM_functions');
addpath('workspaces');
addpath('plot_functions');

%% Problem
% \div ( [D] \grad \theta(x,y)) + s(x,y) = 0
% q * n = - k \grad \theta * n_ = q0    on  \gamma_Q                      (Neumann Condition)
% q * n = - k \grad \theta * n_ = - h*( T - \theta)   on  \gamma_R      (Robin Condition)

% Model params
kw = 0.03; Large = 0.35;    A = 1;    ro = 1.20;  c_v = 717;    h = 8;  k = 0.025;
D = k/(ro*c_v); D = D*10;
Rtot = Large/(kw*A) + 1/(h*A);  Rtot = Rtot * 25;

%% Generate mesh
% % open pdetool and create a 4x3 rectangle, with corners coords readed counter
% % clock-wise: (0,0), (4,0), (4,3), (3,0)
% pderect([0 4 0 3]) 

% % now click on the boundary mode (boundary of omega symbol) and click
% % "Show Edge Labels" in "Boundary" menu, now you can read the label of each
% % segment that will be the in the 5-th row of the mesh matrix e.
% % You will see that the segments are numbered starting from 1 to 4
% % in counter clock-wise way

% % now click on the triangle symbol to triangulize the domain, then click 
% % on "Export Mesh" in the "Mesh" menu
%% Define the Domain 
load domain_toolbox_1
nodes_coordinates_full = p;
idxs_boundaries_full = e;
idxs_elements_full = t;

Lx = 4; Ly = 3; 

%% Forcement and boundary conditions

T = 20;  q0 = 0;
robin_segments_label = [1,2];       neumann_segments_label = [3,4];
s = @(x,y) 0;

%% Solve the problem

[L,K,M,f_omega,f_gammaR1,f_gammaQ] = solve_diffusiveEquation_FEM(D,k,Rtot,s,nodes_coordinates_full, idxs_boundaries_full, idxs_elements_full, robin_segments_label, neumann_segments_label);

f = [f_omega, f_gammaR1, f_gammaQ] * [1; T; q0];

THETA_fem_vector = L \ f;

%% Analytic solution and plots
plot_stationaryResults(THETA_fem_vector, nodes_coordinates_full);
