% 3D linear beam frame natural vibration response
clc
clear 
close all

% Exercise 3.13 of HODGES & PIERCE -  Introduction to Structural Dynamics and Aeroelasticity - (2011)

%% Inputs
% Unit system
unit_sys = "SI-N.m-rads";  % Choose between SI-N.m, SI-kN.m, E-lb.in, E-lb.ft -Hz or -rads
% Beam theory
beam_theory = "EB";    % Choose between "EB" (Euler-Bernoulli) or "T" (Timoshenko)
% Problem parameters - I pick very large values of EA, EIzz, GJ, rhoIs so
% that the lowest frequencies are on bending about the y-axis (EIyy
% relatively low), without the need to BC every node
E0 = 200e9; nu = 0.3; G0 = E0/(2*(1+nu)); rho0 = 7.8e3;
L0 = 10; A0 = 1e9/E0; Iyy0 = 1e2/E0; Izz0 = 1e6*Iyy0; J0 = 1e6/G0; 
Is0 = 1e3/rho0; Gamma0 = 0/E0; e1_0 = 0; e2_0 = 0; Ks0y = 5/6; Ks0z = 5/6;
kappa = linspace(0,30,31); % Vector of kappa values                     
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
Ne_b = 4*ones(N_beams,1);          % Number of elements of each beam
element_order = "quadratic";       % Choose between "linear" or "quadratic"
elem_connect = "sequenced";        % Choose between "sequenced" (beams are connected one after the other) or "unsequenced". If unsequenced, specify the matrix "elem_nodes" below, containing as many rows as elements and as many columns as nodes/element
% Specify geometrical global frame BCs at nodes as a node-BCs tuple (e.g., BC_nodes = {[1, "C"], [2, "SS"], [3, "R1"], [4, "R25", 1e-3, 2e-3], specifies node 1 to be clamped, node 2 to be "simply supported", node 3 to have DOF 1 equal to zero, and node 4 to have its DOFs 2 and 5 set to the values of 1e-3 and 2e-3})
BC_nodes = {[1, "R124567"]};       % Sliding at left end, free at right end
% Option to consider rotary inertia (1 for true, 0 for false)
RI = 1;
% Number of modes sought
N_modes = 3;

%% Call solver 
norm_omega = zeros(length(kappa),N_modes);
for i=1:length(kappa)
    % Set the spring stiffness value
    k0 = kappa(i)*E0*Iyy0/L0^3;
    CSz{1} = k0; akz{1} = 0;
    % Run the program
    disp("Running kappa value " + num2str(i) + " out of " + num2str(length(kappa)));
    LNV3DB_main; close all; 
    % Save normalized natural frequencies to replicate plot of Figure 3.49
    norm_omega(i,:) = SOLdata.omega/sqrt((E0*Iyy0)/(rho0*A0*L0^4));
end

%% Compare with literature
% Replication of Figure 3.49
axes_size = 20; lw = 1; ms = 5;
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontName','times new roman','FontSize',axes_size);
axes1.XLim = [kappa(1) kappa(end)]; axes1.YLim = [0 32];
hold(axes1,'on');
for m=1:N_modes
    plot(kappa,norm_omega(:,m),'k-o','MarkerSize',ms,'Parent',axes1);
end
xlabel('$\kappa$','FontWeight','normal','FontSize',axes_size);
ylabel('$(\alpha_i l)^2$','FontWeight','normal','FontSize',axes_size);
grid on
% Comments: Figure 3.49 is exactly reproduced
