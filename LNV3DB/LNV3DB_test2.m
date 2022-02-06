% 3D linear beam frame natural vibration response
clc
clear 
close all

% Comparison with results from BANARJEE et al. - Exact dynamic stiffness matrix 
% of a bending-torsion coupled beam including warping - (1996) and
% TANAKA & BERCIN -  Finite element modelling of the coupled bending and torsional
% free vibration of uniform beams with an arbitrary cross-section - (1997)

%% Inputs
% Unit system
unit_sys = "SI-N.m-Hz";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "EB";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
E0 = 215e9; nu = 0.25; G0 = E0/(2*(1+nu)); rho0 = 7.8e3;
L0 = 1.28; A0 = 2.095/rho0; Izz0 = 97400/E0; Iyy0 = 0.94e-7; J0 = 11.21/G0; 
Is0 = 0.725e-2/rho0; Gamma0 = 35.4/E0; e1_0 = 0.03771; e2_0 = 0;
Ks0y = 0.357; Ks0z = 0.376;
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
Ne_b = 5*ones(N_beams,1);          % Number of elements of each beam
element_order = "quadratic";       % Choose between "linear" or "quadratic"
elem_connect = "sequenced";        % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BCs = "C-F"; % Select between "C-F", "C-C" or "S-S"
if element_order == "linear", end_node = sum(Ne_b)+1; else, end_node = 2*sum(Ne_b)+1; end
if BCs == "C-F"
    BC_nodes = {[1, "C"]};
elseif BCs == "C-C"
    BC_nodes = {[1, "C"]; [end_node, "C"]};
elseif BCs == "S-S"
    BC_nodes = {[1, "SS"]; [end_node, "SS"]};
elseif BCs == "F-F"
    BC_nodes = {[1, "R13"]; [end_node, "R13"]};    
end
j = length(BC_nodes);
for i=2:end_node
    BC_nodes{j+i-1} = [i, "R135"];
end
% Option to consider rotational inertia (1 for true, 0 for false)
RI = 1;
% Number of modes sought
if BCs == "F-F"
    N_modes = 6;
else
    N_modes = 3;
end

%% Call solver - Solve and plot results
LNV3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post BCs

%% Compare with literature
omega = SOLdata.omega;
if BCs == "C-F"
    ref_freqs = [25.37; 98.57; 149.40];
elseif BCs == "C-C"
    ref_freqs = [149.39; 410.56; 624.60];
elseif BCs == "S-S"
    ref_freqs = [67.13; 263.67; 275.80];
elseif BCs == "F-F"
    omega = omega(4:end); % Discard rigid-body modes
    ref_freqs = [22.04; 152.08; 412.23];    
end  
freq_rel_error = omega./ref_freqs-1