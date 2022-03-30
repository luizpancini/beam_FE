% 3D linear beam frame static response
clc
clear 
close all

%% Inputs
% Unit system
unit_sys = "SI-N.m";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft
% Beam theory
beam_theory = "T";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
L0 = 1; E0 = 200e9; nu = 0.3; G0 = E0/(2*(1+nu));
b = 0.01; H = 0.01; A0 = b*H; Iyy0 = b*H^3/12; Izz0 = H*b^3/12; Ks0y = 5/6; Ks0z = 5/6;
I0 = Iyy0+Izz0; J0 = 0.025*A0^4/I0; Gamma0 = 0;
P0 = 10; 
% Geometry 
N_beams = 1;                                                % Number of beam segments (not the same as elements, e.g., one segment may have several elements)
L = L0*ones(N_beams,1);                                     % Length of the beams
b_alpha = 0*ones(N_beams,1);                                % Angle of rotation of the beams about the z-axis [rad] 
b_beta = 0*ones(N_beams,1);                                 % Angle of rotation the beams about the y'-axis [rad]
b_gamma = 0*ones(N_beams,1);                                % Angle of rotation the beams about the x''-axis [rad]
% Material and section properties - as functions of the beam local coordinate (x)
props.constitutive_model = "aniso" + strings(N_beams,1);    % Constitutive model for beams (only isotropic = "iso" available for now)
props.E = @(x) E0*ones(N_beams,1);                          % Elastic modulus 
props.G = @(x) G0*ones(N_beams,1);                          % Shear modulus
props.A = @(x) A0*ones(N_beams,1);                          % Cross-sectional area
props.J = @(x) J0*ones(N_beams,1);                          % Torsion constant 
props.Gamma = @(x) Gamma0*ones(N_beams,1);                  % Warping constant
props.Iyy = @(x) Iyy0*ones(N_beams,1);                      % Second moment of inertia about y-axis
props.Izz = @(x) Izz0*ones(N_beams,1);                      % Second moment of inertia about z-axis
props.Ksy = @(x) Ks0y*ones(N_beams,1);                      % Shear correction coefficient for y-direction
props.Ksz = @(x) Ks0z*ones(N_beams,1);                      % Shear correction coefficient for z-direction
EA = @(x) props.E(x).*props.A(x);                           % Axial stiffness
GAy = @(x) props.Ksy(x).*props.G(x).*props.A(x);            % Lateral shear stiffness
GAz = @(x) props.Ksz(x).*props.G(x).*props.A(x);            % Transverse shear stiffness
GJ = @(x) props.G(x).*props.J(x);                           % Torsional stiffness
EIy = @(x) props.E(x).*props.Iyy(x);                        % Lateral bending stiffness
EIz = @(x) props.E(x).*props.Izz(x);                        % Transverse bending stiffness
K_ax_oop = @(x) props.E(x).*props.Iyy(x)*10;                % Axial-out-of-plane-bending coupling stiffness
K_oop_tors = @(x) props.G(x).*props.J(x)*1/10;              % Torsion-out-of-plane-bending coupling stiffness
% Anisotropic constitutive matrix
props.C_aniso{1} = @(x) [EA(x)       0      0      0             K_ax_oop(x)        0;
                         0           GAy(x) 0      0             0                  0;
                         0           0      GAz(x) 0             0                  0;
                         0           0      0      GJ(x)         K_oop_tors(x)      0;
                         K_ax_oop(x) 0      0      K_oop_tors(x) EIy(x)             0;
                         0           0      0      0             0             EIz(x)];
% Distributed loads and sources - in the local beam coordinate (x)
% These are inputs in the format: load_of_x{beam} = @(x) ~load_expression~

% Concentraded loads and sources - in the local beam coordinate (x)
% These are inputs in the format: load{beam} = [vector of loads]; load_position{beam} = [vector of loads positions]
Pz{1} = P0;    apz{1} = L0; 
% Option for warping DOF (1 for true or 0 for false)
warp_DOF = 0;
% Set elements data
Ne_b = 10*ones(N_beams,1);          % Number of elements of each beam
element_order = "quadratic";          % Choose between "linear" or "quadratic"
elem_connect = "sequenced";        % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BC_nodes = {[1, "C"]}; 
% Scale for deformed structure plot
scale = 1;

%% Call solver - Solve and plot results
LS3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post