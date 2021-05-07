clear, clc, close all
addpath('FEM_functions');
addpath('workspaces');
addpath('plot_functions');

%% Problem
% \div ( [D] \grad \theta(x,y)) + s(x,y) = 0
% q * n = - k \grad \theta * n_ = q0    on  \gamma_Q                      (Neumann Condition)
% q * n = - k \grad \theta * n_ = - h*( T - \theta)   on  \gamma_R      (Robin Condition)

D = 1;   k = .5;     Rtot = 1;


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

%% Building the manufactured solution
T = 0;  q0 = 0;
robin_segments_label = [1,2,3,4];       neumann_segments_label = [];
p = @(x) x.*(1-1/Lx*x)/(k*Rtot);     g = @(y) y.*(1-1/Ly*y)/(k*Rtot);
p_dot = @(x) (1-2/Lx*x)/(k*Rtot);    g_dot = @(y) (1-2/Ly*y)/(k*Rtot);
p_ddot = @(x) (-2/Lx)/(k*Rtot);           g_ddot = @(y) (-2/Ly)/(k*Rtot);
s = @(x,y) -D*(p_ddot(x) + p_dot(x).^2 + g_ddot(y) + g_dot(y).^2) .* exp(p(x)) .* exp(g(y));
THETA_analyticFunction = @(x,y) exp(p(x)) .* exp(g(y));


%% Solve the problem

[L,K,M,f_omega,f_gammaR1,f_gammaQ] = solve_diffusiveEquation_FEM(D,k,Rtot,s,nodes_coordinates_full, idxs_boundaries_full, idxs_elements_full, robin_segments_label, neumann_segments_label);

f = [f_omega, f_gammaR1, f_gammaQ] * [1; T; q0];

THETA_fem_vector = L \ f;


%% Analytic solution and plots

THETA_analytic_vector = THETA_analyticFunction(nodes_coordinates_full(1,:), nodes_coordinates_full(2,:));
d_approx_anaylitic_vector = THETA_analytic_vector - THETA_fem_vector;   %%distance from the analytic solution

plot_stationaryResults_analyticSolution(THETA_fem_vector, THETA_analytic_vector, d_approx_anaylitic_vector, nodes_coordinates_full);
    
    
%% Error energy norm
errorNorm_L2 = compute_errorNormL2(THETA_fem_vector,THETA_analyticFunction, idxs_elements_full,nodes_coordinates_full,Lx,Ly);
errorNorm_inf = compute_errorNorm_inf(THETA_fem_vector, THETA_analytic_vector');

errorNorm_L2_1 = 0.175690537413827;
errorNorm_L2_2 = 0.045580224937643;
errorNorm_L2_3 = 0.011545437727128;
errorNorm_L2_4 = 0.002898269927076;
errorNorm_inf_1 = 0.179285520507033;
errorNorm_inf_2 = 0.072808785433235;
errorNorm_inf_3 = 0.025208772354777;
errorNorm_inf_4 = 0.008028336685381;

%p = log10(errorNorm_L2_2/errorNorm_L2_3) / log10(2);
p1 = log10(errorNorm_L2_1/errorNorm_L2_2) / log10(2);
p2 = log10(errorNorm_L2_2/errorNorm_L2_3) / log10(2);
p3 = log10(errorNorm_L2_3/errorNorm_L2_4) / log10(2);
figure, loglog([1 2 3 4],[errorNorm_L2_4 errorNorm_L2_3 errorNorm_L2_2 errorNorm_L2_1]), grid on,
hold on, scatter([1 2 3 4],[errorNorm_L2_4 errorNorm_L2_3 errorNorm_L2_2 errorNorm_L2_1],'b','filled');
xlabel('r'), ylabel('log (||e_h||_2)'), xticklabels({'1','','2','','3','','4'}), title('Error norm L2 as function of r');
figure, loglog([1 2 3 4],[errorNorm_inf_4 errorNorm_inf_3 errorNorm_inf_2 errorNorm_inf_1]), grid on,
hold on, scatter([1 2 3 4],[errorNorm_inf_4 errorNorm_inf_3 errorNorm_inf_2 errorNorm_inf_1],'b','filled');
xlabel('r'), ylabel('log (||e_h||_{inf})'), xticklabels({'1','','2','','3','','4'}), title('Error norm h inf as function of r');
