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

%% Params
    
% Model params
kw = 0.03; Large = 0.35;    A = 1;    ro = 1.20;  c_v = 717;    h = 8;  k = 0.025;
D = k/(ro*c_v); D = D*10;
Rtot = Large/(kw*A) + 1/(h*A);  Rtot = Rtot * 25;

% Integration params
Dt = 10;
days = 3;
t_final = 60*60*24 * days;


% Forcement params
global SOURCE_X_MIN SOURCE_X_MAX SOURCE_Y_MIN SOURCE_Y_MAX MU_MAX
global SOURCE_X_MIN_2 SOURCE_X_MAX_2 SOURCE_Y_MIN_2 SOURCE_Y_MAX_2
global RELE_THRESHOLD
SOURCE_X_MIN = 0.1;     SOURCE_X_MIN_2 = 1.6;
SOURCE_X_MAX = 0.3;     SOURCE_X_MAX_2 = 2.4;
SOURCE_Y_MIN = 1.1;     SOURCE_Y_MIN_2 = 0.1;
SOURCE_Y_MAX = 1.9;     SOURCE_Y_MAX_2 = 0.3;
SOURCE_AREA = (SOURCE_X_MAX - SOURCE_X_MIN) * (SOURCE_Y_MAX - SOURCE_Y_MIN);
SOURCE_POWER = 600;        %Power single radiator (watt)
SOURCE_POWER = SOURCE_POWER * 8/600;        %Scaling
MU_MAX = (SOURCE_POWER / SOURCE_AREA) / (ro*c_v);
RELE_THRESHOLD = 1.5;
s = @radiator_forcement;
G_s = tf(MU_MAX, [100*60/4.6 1]);         %low pass filter that model the forcement dynamic
G_z = c2d(G_s,Dt,'tustin');



% Boundary conditions
q0 = 0;    T = 10;    T_peak = 3.5;
freq_T = 1/(60*60*24) * 2*pi;   %rad/s
phase_T = +16/24 * 2*pi;    %rad
%q0 = 0;    T = 20;    T_peak = 0;
robin_segments_label = [1,2];   neumann_segments_label = [3,4];
THETA_0 = zeros(n_dof,1) + T;


%Controlled variable
sensor_coordinates = [2 , 1.5];
idx_sensor = find_nodeIndex(sensor_coordinates,nodes_coordinates_full);
reference = 20;


%% Solve FEM
% M \dot \theta + L \theta = f
% f = [f_omega f_gammaR1 f_gammaQ] * [mu; T; q0];
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
%%% In this integration solvers there is the controller
%THETA_fem_evolution_timeseries = solve_explicitEuler_withController(@input_fnc,THETA_0,T,freq_T,phase_T,q0,G_z,reference, idx_sensor, Dt,0,t_final);
%THETA_fem_evolution_timeseries = solve_crankNicolson_withController(@input_fnc,THETA_0,T,freq_T,phase_T,q0,G_z,reference, idx_sensor, Dt,0,t_final);
%THETA_fem_evolution_timeseries = solve_RungeKutta_withController(@input_fnc,THETA_0,T,freq_T,phase_T,q0,G_z,reference, idx_sensor, Dt,0,t_final);
THETA_fem_evolution_timeseries = sim('simu_unsteady',t_final).THETA_evolution;

%% Plot evolution
THETA_fem_evolution = THETA_fem_evolution_timeseries.Data;
time = THETA_fem_evolution_timeseries.Time;

figure, subplot(2,1,1), plot(time/3600, THETA_fem_evolution), grid on; ylabel({'$ \Theta$ '},'Interpreter','latex'), xlabel('Time (hours)');
ylim([-1 45]); title('Nodal temperature evolutions');
subplot(2,1,2), plot(time/3600, THETA_fem_evolution(:,idx_sensor)), grid on, ylabel('Controlled node'), xlabel('Time (hours)')

%% Plot 3d over time
plot_evolution(THETA_fem_evolution,time, nodes_coordinates_full, Lx,Ly, 10);



