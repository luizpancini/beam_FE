% 3D linear beam frame natural vibration response
clc
clear 
close all

% Comparison with results from 
% HJAJI & MOHAREB - Finite-element formulation for the linear steady-state response
% of asymmetric thin-walled members under harmonic forces - (2015)

%% Inputs
% Unit system
unit_sys = "SI-N.m-Hz";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "T";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
EXAMPLE = 2; % Choose example 1 or 2
if EXAMPLE == 1
    L0 = 4;
elseif EXAMPLE == 2
    L0 = 1;
end
E0 = 200e9; nu = 0.3; G0 =  E0/(2*(1+nu)); rho0 = 7.85e3;
A0 = 0.2e-2; Izz0 = 0.878e-6; Iyy0 = 3.723e-6; J0 = 0.5707e-7; 
Gamma0 = 0.8608e-9; e2_0 = 42.83e-3; e1_0 = 10.29e-3; Is0 = Iyy0+Izz0+A0*(e1_0^2+e2_0^2);
Ks0y = 0.3; Ks0z = 0.3; % One may change these to 0.1 to better match the reference results
% Geometry 
N_beams = 1;                                        % Number of beam segments (not the same as elements, e.g., one segment may have several elements)
L = L0*ones(N_beams,1);                             % Length of the beams
b_alpha = 0*ones(N_beams,1);                        % Angle of rotation of the beams about the z-axis [rad] 
b_beta = 0*ones(N_beams,1);                         % Angle of rotation the beams about the y'-axis [rad]
b_gamma = 0*ones(N_beams,1);                        % Angle of rotation the beams about the x''-axis [rad]
% Material and section properties - as functions of the beam local coordinate (x)
constitutive_model = "iso" + strings(N_beams,1);    % Constitutive model for beams (only isotropic = "iso" available for now)
E = @(x) E0*ones(N_beams,1);                        % Elastic modulus 
G = @(x) G0*ones(N_beams,1);                        % Shear modulus
rho = @(x) rho0*ones(N_beams,1);                    % Mass density
A = @(x) A0*ones(N_beams,1);                        % Cross-sectional area
J = @(x) J0*ones(N_beams,1);                        % Torsion constant 
Gamma = @(x) Gamma0*ones(N_beams,1);                % Warping constant
Is = @(x) Is0*ones(N_beams,1);                      % Polar moment of inertia 
Iyy = @(x) Iyy0*ones(N_beams,1);                    % Second moment of inertia about y-axis
Izz = @(x) Izz0*ones(N_beams,1);                    % Second moment of inertia about z-axis
Ksy = @(x) Ks0y*ones(N_beams,1);                    % Shear correction coefficient for y-direction
Ksz = @(x) Ks0z*ones(N_beams,1);                    % Shear correction coefficient for z-direction
e1 = @(x) e1_0*ones(N_beams,1);                     % z-axis distance from centroid to elastic axis / shear center
e2 = @(x) e2_0*ones(N_beams,1);                     % y-axis distance from centroid to elastic axis / shear center
% Set elements data
Ne_b = 10*ones(N_beams,1);          % Number of elements of each beam
element_order = "quadratic";        % Choose between "linear" or "quadratic"
elem_connect = "sequenced";         % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BC_nodes = {[1, "C"]};
% Option to consider rotational inertia (1 for true, 0 for false)
RI = 1;
% Number of modes sought
if EXAMPLE == 1
    N_modes = 8;
elseif EXAMPLE == 2
    N_modes = 6;
end

%% Call solver - Solve and plot results
LNV3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post EXAMPLE

%% Compare with literature
if EXAMPLE == 1
    ref_freqs = [3.68; 7.265; 22.35; 23.75; 40.37; 62.82; 73.86; 95.07]; 
elseif EXAMPLE == 2
    ref_freqs = [57.20; 72.32; 172.9; 264.1; 358.2; 550.7]; 
end
freq_rel_error = SOLdata.omega./ref_freqs-1