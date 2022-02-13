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
E0 = 200e9; nu0 = 0.3; G0 = E0/(2*(1+nu0)); 
L0 = 1; b = 0.01; H = 0.01; A0 = b*H; Izz0 = H*b^3/12; Iyy0 = b*H^3/12; I0 = Iyy0+Izz0; J0 = 0.025*A0^4/I0; Gamma0 = 10/E0; Ks0y = 5/6; Ks0z = 5/6; 
P0 = 10; q0 = 10; T0 = 10; k0 = 1e-0*E0*I0;
% Geometry 
N_beams = 3;                                        % Number of beam segments (not the same as elements, e.g., one segment may have several elements)
L = L0*ones(N_beams,1);                             % Length of the beams
b_alpha = [0 pi/2 0];                               % Angle of rotation of the beams about the z-axis [rad] 
b_beta = [0 0 -pi/2];                               % Angle of rotation the beams about the y'-axis [rad]
b_gamma = 0*ones(N_beams,1);                        % Angle of rotation the beams about the x''-axis [rad]
% Material and section properties - as functions of the beam local coordinate (x)
constitutive_model = "iso" + strings(N_beams,1);    % Constitutive model for beams (only isotropic = "iso" available for now)
E = @(x) E0*ones(N_beams,1);                        % Elastic modulus 
G = @(x) G0*ones(N_beams,1);                        % Shear modulus
A = @(x) A0*ones(N_beams,1);                        % Cross-sectional area
J = @(x) J0*ones(N_beams,1);                        % Torsion constant 
Gamma = @(x) Gamma0*ones(N_beams,1);                % Warping constant
Iyy = @(x) Iyy0*ones(N_beams,1);                    % Second moment of inertia about y-axis
Izz = @(x) Izz0*ones(N_beams,1);                    % Second moment of inertia about z-axis
Ksy = @(x) Ks0y*ones(N_beams,1);                    % Shear correction coefficient for y-direction
Ksz = @(x) Ks0z*ones(N_beams,1);                    % Shear correction coefficient for z-direction
% Distributed loads and sources - in the local beam coordinate (x)
% These are inputs in the format: load_of_x{beam} = @(x) ~load_expression~
tq_of_x{3} = @(x) q0*x.^2;
ql_of_x{1} = @(x) q0*(L0-x); 
qz_of_x{2} = @(x) q0;
cft_of_x{2} = @(x) k0;
% Concentraded loads and sources - in the local beam coordinate (x)
% These are inputs in the format: load{beam} = [vector of loads]; load_position{beam} = [vector of loads positions]
CSu{1} = [k0];          aku{1} = [L0];
Pa{3} = [-P0];          apa{3} = [L0];
Pt{3} = [-P0/4; P0/2];  apt{3} = [L0; L0/2];
Tq{1} = [T0];           atq{1} = [L0/2];
Mz{2} = [T0];           amz{2} = [.85*L0];
% Option for warping DOF (1 for true or 0 for false)
warp_DOF = 1;
% Set elements data
Ne_b = 5*ones(N_beams,1);          % Number of elements of each beam
element_order = "quadratic";       % Choose between "linear" or "quadratic"
elem_connect = "sequenced";        % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BC_nodes = {[1, "C"]}; 
% Scale for deformed structure plot
scale = 1;

%% Call solver - Solve and plot results
LS3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post 