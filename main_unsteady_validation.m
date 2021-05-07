clear all, clc, close all
addpath('FEM_functions');
addpath('timeIntegration_functions');
addpath('workspaces');
addpath('plot_functions');

%% Problem
% \div ( [D] \grad \theta(x,y)) + s(x,y) * mu(t) = \der \theta
% q * n = - k \grad \theta * n_ = q0    on  \gamma_Q                      (Neumann Condition)
% q * n = - k \grad \theta * n_ = - 1/Rtot*( T - \theta)   on  \gamma_R      (Robin Condition)

%% Define the Domain 
load domain_toolbox_1
nodes_coordinates_full = p;
idxs_boundaries_full = e;
idxs_elements_full = t;
n_dof = size(nodes_coordinates_full,2);
Lx = 4; Ly = 3; 

%% Building the manufactured solution

% Model params
D = 1;   k = 1;     Rtot = 1;

%Boundary params
T = 0;  robin_segments_label = [1,2,3,4];   neumann_segments_label = [];    q0 = 0;

% Method of Manufactured solution to perform the analytic test
p = @(x) x.*(1-1/Lx*x)/(Rtot*k);     g = @(y) y.*(1-1/Ly*y)/(Rtot*k);
p_dot = @(x) (1-2/Lx*x)/(Rtot*k);    g_dot = @(y) (1-2/Ly*y)/(Rtot*k);
p_ddot = @(x) (-2/Lx)/(Rtot*k);           g_ddot = @(y) (-2/Ly)/(Rtot*k);
q = @(t) -t/20;                   q_dot = -1/20;
s = @(x,y) -(p_ddot(x) + p_dot(x)^2 + g_ddot(y) + g_dot(y)^2 - q_dot) .* exp(p(x)) .* exp(g(y));
THETA_analyticFunction = @(x,y,t) exp(p(x)) .* exp(g(y)) .* exp(q(t));
THETA_0 = (exp(p(nodes_coordinates_full(1,:))) .* exp(g(nodes_coordinates_full(2,:))) .* exp(q(0)))';



% Integration params
Dt = 0.002;
t_final = 30;

%% Solve FEM
% M \dot \theta + L \theta = f
% f = [f_omega f_gammaR1 f_gammaQ] * [u; T; q0];
% \dot \theta = - inv(M) * L \theta + inv(M) * f
% \dot \theta = A \theta + b                  (equivalent LTI system)
% where:    A = - inv(M) * L,   
% where: b = inv(M) * f = inv(M) * [f_omega f_gammaR1 f_gammaQ] * [mu; T; q0]
% where: b = B*u;

[L,K,M,f_omega,f_gammaR1,f_gammaQ] = solve_diffusiveEquation_FEM(D,k,Rtot,s,nodes_coordinates_full, idxs_boundaries_full, idxs_elements_full, robin_segments_label, neumann_segments_label);


%% We have obtained an LTI system that is rapresented by these matricies
global A_lti B_lti
A_lti = - M \ L;
B_lti = M \ [f_omega f_gammaR1 f_gammaQ];
C_lti = eye(size(A_lti));
D_lti = zeros(size(B_lti));
%u = [0; T; q0];                    %u = [mu; T; q0];


%% Integration over time
stabilityAnalysis_explEuler(A_lti, Dt);
THETA_fem_evolution_timeseries = solve_explicitEuler_analyticSolution(@input_fnc,q,THETA_0,T,q0, Dt,0,t_final);
%THETA_fem_evolution_timeseries = solve_crankNicolson_analyticSolution(@input_fnc,q,THETA_0,T,q0, Dt,0,t_final);
%THETA_fem_evolution_timeseries = solve_RungeKutta_analyticSolution(@input_fnc,q,THETA_0,T,q0, Dt,0,t_final);
%THETA_fem_evolution_timeseries = sim('simu_analyticSolution',t_final).THETA_evolution;

%% Plot
THETA_fem_evolution = THETA_fem_evolution_timeseries.Data;
time = THETA_fem_evolution_timeseries.Time;

figure,
plot(time, THETA_fem_evolution), grid on, ylabel({'$ \Theta$ '},'Interpreter','latex'), xlabel('Time (seconds)');
ylim([-1 21]); title('Nodal temperature evolutions');

THETA_analytic_evolution = THETA_analyticFunction(nodes_coordinates_full(1,:), nodes_coordinates_full(2,:), time);
difference_evolution = abs(THETA_fem_evolution - THETA_analytic_evolution);
plot_evolution_analyticSolution(THETA_fem_evolution,THETA_analytic_evolution,difference_evolution,time, nodes_coordinates_full, Lx,Ly, 500);


%% Error

THETA_fem_vector = (THETA_fem_evolution(end,:))';      %THETA(x,y; t_final);
THETA_analyticFunction_tfinal = @(x,y) THETA_analyticFunction(x,y,t_final);
errorNorm_L2 = compute_errorNormL2(THETA_fem_vector,THETA_analyticFunction_tfinal, idxs_elements_full,nodes_coordinates_full,Lx,Ly);


% errorNorm_L2 = 0.0021; %CrankNic: dT = 0.2, t_final = 40, domain2
% errorNorm_inf = 0.0049; %CrankNic: dT = 0.2, t_final = 40, domain2

% errorNorm_L2 = 9.2486e-04; %CrankNic: dT = 0.1, t_final = 40, domain2
% errorNorm_inf = 0.0031; %CrankNic: dT = 0.1, t_final = 40, domain2

% errorNorm_L2 = 0.0018; %CrankNic: dT = 0.2, t_final = 40, domain1
% errorNorm_inf = 0.0072; %CrankNic: dT = 0.2, t_final = 40, domain1

% errorNorm_L2 = 0.0016; %CrankNic: dT = 0.1, t_final = 40, domain1
% errorNorm_inf = 0.0054; %CrankNic: dT = 0.1, t_final = 40, domain1

%sol = Newton_Raphson(@observerdOrder_fnc,[10 10 2 2]', 10000, 1e-4);

