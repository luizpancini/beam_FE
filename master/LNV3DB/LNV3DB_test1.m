% 3D linear beam frame natural vibration response
clc
clear 
close all

% Comparison with results from 
% TANAKA & BERCIN -  Finite element modelling of the coupled bending and torsional
% free vibration of uniform beams with an arbitrary cross-section - (1997)

%% Inputs
% Unit system
unit_sys = "SI-N.m-Hz";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "EB";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
E0 = 200e9; nu = 0.3; G0 = E0/(2*(1+nu)); rho0 = 7.8e3;
L0 = 1.5; A0 = 1.947/rho0; Iyy0 = 16680/E0; Izz0 = 73480/E0; J0 = 10.81/G0; 
Is0 = 76e-8; Gamma0 = 26.34/E0; e1_0 = 0.02625; e2_0 = 0.02316; Ks0y = 1; Ks0z = 1;
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
Ne_b = 12*ones(N_beams,1);          % Number of elements of each beam
element_order = "linear";          % Choose between "linear" or "quadratic"
elem_connect = "sequenced";        % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BCs = "C-F"; % Select between "C-F", "C-C", "S-S" or "F-F"
if element_order == "linear", end_node = sum(Ne_b)+1; else, end_node = 2*sum(Ne_b)+1; end
if BCs == "C-F"
    BC_nodes = {[1, "C"]};
elseif BCs == "C-C"
    BC_nodes = {[1, "C"]; [end_node, "C"]};
elseif BCs == "S-S"
    BC_nodes = {[1, "SS"]; [end_node, "SS"]};
elseif BCs == "F-F"
    BC_nodes = {[1, "R1"]; [end_node, "R1"]};    
end
% Option to consider rotational inertia (1 for true, 0 for false)
RI = 1;
% Number of modes sought
if BCs == "F-F"
    N_modes = 8;
else
    N_modes = 3;
end

%% Call solver - Solve and plot results
LNV3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post BCs

%% Compare with literature
omega = SOLdata.omega;
if BCs == "C-F"
    ref_freqs = [17.17; 27.31; 59.10];
elseif BCs == "C-C"
    ref_freqs = [98.71; 169.43; 270.89];
elseif BCs == "S-S"
    ref_freqs = [44.71; 75.14; 164.86];
elseif BCs == "F-F"
    omega = SOLdata.omega(end-2:end); % Discard rigid-body modes
    ref_freqs = [20.35; 101.26; 170.31];    
end  
freq_rel_error = omega./ref_freqs-1