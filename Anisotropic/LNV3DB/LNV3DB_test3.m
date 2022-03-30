% 3D linear beam frame natural vibration response
clc
clear 
close all

% Comparison with results from 
% BERCIN & TANAKA - COUPLED FLEXURAL-TORSIONAL VIBRATIONS OF TIMOSHENKO BEAMS - (1997)

%% Inputs
% Unit system
unit_sys = "SI-N.m-Hz";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "T";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
EXAMPLE = 2;
if EXAMPLE == 1
    E0 = 200e9; nu = 0.25; G0 = E0/(2*(1+nu)); rho0 = 7.8e3;
    L0 = 0.82; A0 = 0.835/rho0; Izz0 = 6380/E0; Iyy0 = 0.94e-7; J0 = 43.46/G0; 
    Is0 = 5e-4/rho0; Gamma0 = 0.10473/E0; e1_0 = 0.0155; e2_0 = 0;
    Ks0y = 4.081e6/(G0*A0); Ks0z = 1;
elseif EXAMPLE == 2
    E0 = 30e9; nu = 0.18; G0 = E0/(2*(1+nu)); rho0 = 2.5e3; % For Timoshenko beams, the actual values of E and rho influence the solution
    L0 = 3.2; A0 = 225/rho0; Izz0 = 30.43e7/E0; Iyy0 = 0.94e-7; J0 = 97.83e4/G0; 
    Is0 = 56.87/rho0; Gamma0 = 81.58e5/E0; e1_0 = 0.336; e2_0 = 0;
    Ks0y = 52.82e7/(G0*A0); Ks0z = 1;
elseif EXAMPLE == 3
    E0 = 200e9; nu = 0.25; G0 = E0/(2*(1+nu)); rho0 = 7.8e3;
    L0 = 2.7; A0 = 4.256/rho0; Izz0 = 14.36e4/E0; Iyy0 = 0.94e-7; J0 = 346.71/G0; 
    Is0 = 3.17e-2/rho0; Gamma0 = 536.51/E0; e1_0 = 0.0735; e2_0 = 0;
    Ks0y = 20.81e6/(G0*A0); Ks0z = 1;
end
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
Ne_b = 3*ones(N_beams,1);         % Number of elements of each beam
element_order = "quadratic";      % Choose between "linear" or "quadratic"
elem_connect = "sequenced";       % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
if element_order == "linear", end_node = sum(Ne_b)+1; else, end_node = 2*sum(Ne_b)+1; end
BC_nodes = {[1, "C"]};
for i=2:end_node
    BC_nodes{i} = [i, "R135"];
end
% Option to consider rotational inertia (1 for true, 0 for false)
RI = 1;
% Number of modes sought
N_modes = 5;

%% Call solver - Solve and plot results
LNV3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post EXAMPLE

%% Compare with literature
if EXAMPLE == 1
    ref_freqs = [63.51; 137.39; 275.82; 481.10; 639.76];
elseif EXAMPLE == 2
    ref_freqs = [23.79; 78.26; 124.78; 295.26; 334.88];
elseif EXAMPLE == 3
    ref_freqs = [11.01; 38.93; 57.82; 150.51; 205.32];
end
freq_rel_error = SOLdata.omega./ref_freqs-1