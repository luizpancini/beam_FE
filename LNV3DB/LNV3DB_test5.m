% 3D linear beam frame natural vibration response
clc
clear 
close all

% Comparison with results from 
% DENNIS & JONES - Flexural-torsional vibration of a tapered C-section beam - (2017)

%% Inputs
% Unit system
unit_sys = "SI-N.m-Hz";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "T";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters 
case_select = "1s";
if case_select == "1s" 
    h = 12.7e-3; b = 25.4e-3; t = 0.635e-3; J_root = 1/3*(h*t^3+2*b*t^3); Gamma_root = 1/6*h^2*b^3*t; Izz_root = t*h^3/12+b*t*h^2/2;
    E0 = 205e9; G0 = 79.3e9; rho0 = 7.855e3;
    L0 = 0.914; A0 = 39.52e-6; Izz0 = 2643e-12; Iyy0 = 1254e-12; J0 = 5.54e-12; 
    Gamma0 = 74997e-18; e2_0 = (10.43+11.25)*1e-3; e1_0 = 0; Is0 = Iyy0+Izz0+A0*e2_0^2; Ks0y = 1; Ks0z = 1;
elseif case_select == "2s"
    E0 = 205e9; G0 = 79.3e9; rho0 = 7.855e3;
    L0 = 0.914; A0 = 40.53e-6; Izz0 = 1857e-12; Iyy0 = 4126e-12; J0 = 5.69e-12; 
    Gamma0 = 176698e-18; e2_0 = (6.79+8.2)*1e-3; e1_0 = 0; Is0 = Iyy0+Izz0+A0*e2_0^2; Ks0y = 1; Ks0z = 1;
elseif case_select == "3a"
    E0 = 68.3e9; G0 = 26.2e9; rho0 = 2.7e3;
    L0 = 0.914; A0 = 27.42e-6; Izz0 = 973e-12; Iyy0 = 205e-12; J0 = 3.79e-12; 
    Gamma0 = 6818e-18; e2_0 = (8.44+8.47)*1e-3; e1_0 = 0; Is0 = Iyy0+Izz0+A0*e2_0^2; Ks0y = 1; Ks0z = 1;
elseif case_select == "4s"
    E0 = 205e9; G0 = 79.3e9; rho0 = 7.855e3; L0 = 0.914; Ks0y = 1; Ks0z = 1;
    A_root = 0.395E-4; Iyy_root = 0.354e-8; Izz_root = 0.181e-8; J_root = 0.532e-11; 
    Gamma_root = 0.151e-12; e2_root = (6.95+8.30)*1e-3; e1_root = 0; Is_root = Iyy_root+Izz_root+A_root*e2_root^2;
    A_tip = 0.395E-4; Iyy_tip = 0.125e-8; Izz_tip = 0.264e-8; J_tip = 0.532e-11; 
    Gamma_tip = 0.751e-13; e2_tip = (10.43+11.25)*1e-3; e1_tip = 0; Is_tip = Iyy_tip+Izz_tip+A_tip*e2_tip^2;
    A_mid = 0.395E-4; Iyy_mid = 0.224e-8; Izz_mid = 0.230e-8; J_mid = 0.532e-11; 
    Gamma_mid = 0.117e-12; e2_mid = (8.60+9.76)*1e-3; e1_mid = 0; Is_mid = Iyy_mid+Izz_mid+A_mid*e2_mid^2;
elseif case_select == "5s"
    E0 = 205e9; G0 = 79.3e9; rho0 = 7.855e3; L0 = 0.914; Ks0y = 1; Ks0z = 1;
    A_root = 0.102E-3; Izz_root = 0.739e-8; Iyy_root = 0.428e-7; J_root = 0.341e-10; 
    Gamma_root = 0.321e-11; e2_root = (7.21+9.43)*1e-3; e1_root = 0; Is_root = Iyy_root+Izz_root+A_root*e2_root^2;
    A_tip = 0.641e-4; Izz_tip = 0.465e-8; Iyy_tip = 0.193e-8; J_tip = 0.214e-10; 
    Gamma_tip = 0.127e-12; e2_tip = (11.21+11.58)*1e-3; e1_tip = 0; Is_tip = Iyy_tip+Izz_tip+A_tip*e2_tip^2; 
    A_mid = 0.832e-4; Izz_mid = 0.633e-8; Iyy_mid = 0.148e-7; J_mid = 0.277e-10; 
    Gamma_mid = 0.106e-11; e2_mid = (8.75+10.43)*1e-3; e1_mid = 0; Is_mid = Iyy_mid+Izz_mid+A_mid*e2_mid^2; 
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
Ksy = @(x) Ks0y*ones(N_beams,1);                    % Shear correction coefficient for y-direction
Ksz = @(x) Ks0z*ones(N_beams,1);                    % Shear correction coefficient for z-direction
if case_select == "1s" || case_select == "2s" || case_select == "3a"
    A = @(x) A0;            % Cross-sectional area
    Iyy = @(x) Iyy0;        % Second moment of inertia about y-axis
    Izz = @(x) Izz0;        % Second moment of inertia about z-axis
    J = @(x) J0;            % Torsion constant 
    Is = @(x) Is0;          % Polar moment of inertia 
    Gamma = @(x) Gamma0;    % Warping constant
    e1 = @(x) e1_0;         % z-axis distance from centroid to elastic axis / shear center
    e2 = @(x) e2_0;         % y-axis distance from centroid to elastic axis / shear center
else
    A = @(x) (A_root+x/(L/2)*(A_mid-A_root)) .* (x <= L/2) + (A_mid+(x-L/2)/(L/2)*(A_tip-A_mid)) .* (x > L/2);                              % Cross-sectional area
    Iyy = @(x) (Iyy_root+x/(L/2)*(Iyy_mid-Iyy_root)) .* (x <= L/2) + (Iyy_mid+(x-L/2)/(L/2)*(Iyy_tip-Iyy_mid)) .* (x > L/2);                % Second moment of inertia about y-axis
    Izz = @(x) (Izz_root+x/(L/2)*(Izz_mid-Izz_root)) .* (x <= L/2) + (Izz_mid+(x-L/2)/(L/2)*(Izz_tip-Izz_mid)) .* (x > L/2);                % Second moment of inertia about z-axis
    J = @(x) (J_root+x/(L/2)*(J_mid-J_root)) .* (x <= L/2) + (J_mid+(x-L/2)/(L/2)*(J_tip-J_mid)) .* (x > L/2);                              % Torsion constant 
    Is = @(x) (Is_root+x/(L/2)*(Is_mid-Is_root)) .* (x <= L/2) + (Is_mid+(x-L/2)/(L/2)*(Is_tip-Is_mid)) .* (x > L/2);                       % Polar moment of inertia 
    Gamma = @(x) (Gamma_root+x/(L/2)*(Gamma_mid-Gamma_root)) .* (x <= L/2) + (Gamma_mid+(x-L/2)/(L/2)*(Gamma_tip-Gamma_mid)) .* (x > L/2);  % Warping constant
    e1 = @(x) (e1_root+x/(L/2)*(e1_mid-e1_root)) .* (x <= L/2) + (e1_mid+(x-L/2)/(L/2)*(e1_tip-e1_mid)) .* (x > L/2);                       % z-axis distance from centroid to elastic axis / shear center
    e2 = @(x) (e2_root+x/(L/2)*(e2_mid-e2_root)) .* (x <= L/2) + (e2_mid+(x-L/2)/(L/2)*(e2_tip-e2_mid)) .* (x > L/2);                       % y-axis distance from centroid to elastic axis / shear center
end
% Set elements data
Ne_b = 10*ones(N_beams,1);         % Number of elements of each beam
element_order = "quadratic";       % Choose between "linear" or "quadratic"
elem_connect = "sequenced";        % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BC_nodes = {[1, "C"]};
% Option to consider rotational inertia (1 for true, 0 for false)
RI = 1;
% Number of modes sought
if case_select == "5s"
    N_modes = 9;
else
    N_modes = 8;
end

%% Call solver - Solve and plot results
LNV3DB_main;
clearvars -except edata FEMdata SOLdata time_total time_pre time_solve time_post case_select

%% Compare with literature
if case_select == "1s" 
    ref_freqs = [13.91; 28.80; 52.58; 69.06; 123.1; 180.6; 217.5; 339.0];
elseif case_select == "2s" 
    ref_freqs = [20.23; 23.91; 65.31; 87.66; 151.3; 211.8; 374.1; 386.6];
elseif case_select == "3a" 
    ref_freqs = [10.08; 20.31; 37.50; 64.84; 88.59; 128.4; 143.4; 197.5];
elseif case_select == "4s"   
    ref_freqs = [20.16; 25.63; 60.16; 83.05; 157.2; 165.9; 284.1; 356.1];
elseif case_select == "5s"   
    ref_freqs = [35.69; 43.63; 99.56; 151.9; 200.3; 252.6; 441.0; 525.8; 530.7];    
end  
freq_rel_error = SOLdata.omega./ref_freqs-1