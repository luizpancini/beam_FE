% 3D linear beam frame natural vibration response
clc
clear 
close all

% Example 3.9 from PETYT - Introduction to Finite Element Vibration Analysis - [2nd Ed.] (2010)

%% Inputs
% Unit system
unit_sys = "SI-N.m-Hz";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "EB";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
E0 = 219.9e9; nu = 0.25; G0 = E0/(2*(1+nu)); rho0 = 7.9e3; L0 = 1;
ba = 0.05; Ha = 0.05; Aa = ba*Ha; Izza = Ha*ba^3/12; Iyya = ba*Ha^3/12; Ia = Iyya+Izza; Ja = 0.025*Aa^4/Ia; 
bb = 0.15; Hb = 0.05; Ab = bb*Hb; Izzb = Hb*bb^3/12; Iyyb = bb*Hb^3/12; Ib = Iyyb+Izzb; Jb = 0.025*Ab^4/Ib;
Gamma0 = 0; e1_0 = 0; e2_0 = 0; Isa = Iyya+Izza+Aa*(e1_0^2+e2_0^2); Isb = Iyyb+Izzb+Ab*(e1_0^2+e2_0^2); Ks0y = 1; Ks0z = 1;
% Geometry 
N_beams = 16;                                       % Number of beam segments (not the same as elements, e.g., one segment may have several elements)
L = L0*ones(N_beams,1);                             % Length of the beams
b_beta = -[pi/2; pi/2; 0; -pi/2; -pi/2; 0;
           0; 0; 0; 0;
           -pi/2; pi/2; 0; -pi/2; -pi/2; 0];        % Angle of rotation of the beams about the z-axis [rad] 
b_alpha = [0; 0; pi/2; 0; 0; pi/2;
          0; 0; 0; 0;
          0; 0; pi/2; 0; 0; pi/2];                  % Angle of rotation the beams about the y'-axis [rad]
b_gamma = 0*ones(N_beams,1);                        % Angle of rotation the beams about the x''-axis [rad]
% Material and section properties - as functions of the beam local coordinate (x)
constitutive_model = "iso" + strings(N_beams,1);    % Constitutive model for beams (only isotropic = "iso" available for now)
E = @(x) E0*ones(N_beams,1);                        % Elastic modulus 
G = @(x) G0*ones(N_beams,1);                        % Shear modulus
rho = @(x) rho0*ones(N_beams,1);                    % Mass density
A = @(x) [Aa; Aa; Ab; Aa; Aa; Ab;
          Ab; Ab; Ab; Ab;
          Aa; Aa; Ab; Aa; Aa; Ab];                  % Cross-sectional area
J = @(x) [Ja; Ja; Jb; Ja; Ja; Jb;
          Jb; Jb; Jb; Jb;
          Ja; Ja; Jb; Ja; Ja; Jb];                  % Torsion constant 
Gamma = @(x) Gamma0*ones(N_beams,1);                % Warping constant
Is = @(x) [Isa; Isa; Isb; Isa; Isa; Isb;
          Isb; Isb; Isb; Isb;
          Isa; Isa; Isb; Isa; Isa; Isb];            % Polar moment of inertia 
Iyy = @(x) [Iyya; Iyya; Iyyb; Iyya; Iyya; Iyyb;
            Iyyb; Iyyb; Iyyb; Iyyb;
            Iyya; Iyya; Iyyb; Iyya; Iyya; Iyyb];    % Second moment of inertia about y-axis
Izz = @(x) [Izza; Izza; Izzb; Izza; Izza; Izzb;
            Izzb; Izzb; Izzb; Izzb;
            Izza; Izza; Izzb; Izza; Izza; Izzb];    % Second moment of inertia about z-axis
Ksy = @(x) Ks0y*ones(N_beams,1);                    % Shear correction coefficient for y-direction
Ksz = @(x) Ks0z*ones(N_beams,1);                    % Shear correction coefficient for z-direction
e1 = @(x) e1_0*ones(N_beams,1);                     % z-axis distance from centroid to elastic axis / shear center
e2 = @(x) e2_0*ones(N_beams,1);                     % y-axis distance from centroid to elastic axis / shear center
% Set elements data
Ne_b = 2*ones(N_beams,1);          % Number of elements of each beam
element_order = "quadratic";       % Choose between "linear" or "quadratic"
elem_connect = "unsequenced";      % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
elem_nodes = [1 2 3;    % Left 1st vertical upwards beam
              3 4 5;
              5 6 7;    % Left 2nd vertical upwards beam 
              7 8 9
              9 10 11;  % Left horizontal high beam
              11 12 13;
              13 14 15; % Left 1st vertical downwards beam
              15 16 17;
              17 18 19; % Left 2nd vertical downwards beam
              19 20 21;
              5 22 23;  % Left horizontal low beam
              23 24 17;
              5 25 26;  % Front low beam
              26 27 28;
              17 29 30; % Back low beam
              30 31 32;
              9 33 34;  % Front high beam
              34 35 36;
              13 37 38; % Back high beam
              38 39 40;
              28 41 42; % Right 1st vertical upwards beam
              42 43 44;
              28 45 46; % Right 2nd vertical upwards beam
              46 47 36;
              36 48 49; % Right horizontal high beam
              49 50 40;
              40 51 52; % Right 1st vertical downwards beam
              52 53 32;
              32 54 55; % Right 2nd vertical downwards beam
              55 56 57;
              28 58 59; % Right horizontal low beam 
              59 60 32];
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BC_nodes = {[1, "C"], [21, "C"], [44, "C"], [57, "C"]};
% Option to consider rotational inertia (1 for true, 0 for false)
RI = 0;
% Number of modes sought
N_modes = 4;

%% Call solver - Solve and plot results
LNV3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post 

%% Compare with literature
ref_freqs = [11.8; 34.1];
freq_rel_error = [SOLdata.omega(1); SOLdata.omega(4)]./ref_freqs-1