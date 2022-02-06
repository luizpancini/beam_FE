function FEMdata = LNV3DB_interp_funs(FEMdata)

%% Interpolation functions
if FEMdata.element_order == "linear"
    % zeta - Lagrange linear polynomials for axial displacement (and torsional rotations in case of no warping DOF)
    zeta_1 = @(xi) 1/2*(1-xi);
    zeta_2 = @(xi) 1/2*(1+xi);
    zeta = {zeta_1; zeta_2};
    dzeta_1 = @(xi) -1/2*ones(length(xi),1); % *ones(length(xi),1) is just a trick to output dzeta as a vector in case xi is input as a vector
    dzeta_2 = @(xi) 1/2*ones(length(xi),1);
    dzeta = {dzeta_1; dzeta_2};
    % psi - generalized Hermite cubic polynomials for bending displacements (and torsional rotation in case of warping DOF)
    psi_1 = @(xi,h_e,Gam) 1./(1+12*Gam).*((1+12*Gam)-6*Gam.*(xi+1)+1/4.*(xi-2).*(xi+1).^2);
    psi_2 = @(xi,h_e,Gam) h_e./(8*(1+12*Gam)).*(xi.^2-1).*(12*Gam-xi+1);
    psi_3 = @(xi,h_e,Gam) 1./(1+12*Gam).*(6*Gam.*(xi+1)-1/4.*(xi-2).*(xi+1).^2);
    psi_4 = @(xi,h_e,Gam) -h_e./(8*(1+12*Gam)).*(xi.^2-1).*(12*Gam+xi+1);
    psi = {psi_1; psi_2; psi_3; psi_4};
    dpsi_1 = @(xi,h_e,Gam) -1./(4*(1+12*Gam)).*(-3*xi.^2+24*Gam+3);
    dpsi_2 = @(xi,h_e,Gam) h_e./(8*(1+12*Gam)).*(-3*xi.^2+2*xi+1+24*Gam.*xi);
    dpsi_3 = @(xi,h_e,Gam) 1./(4*(1+12*Gam)).*(-3*xi.^2+24*Gam+3);
    dpsi_4 = @(xi,h_e,Gam) -h_e./(8*(1+12*Gam)).*(3*xi.^2+2*xi+24*Gam.*xi-1);
    dpsi = {dpsi_1; dpsi_2; dpsi_3; dpsi_4};
    d2psi_1 = @(xi,h_e,Gam) 3/2*xi./(12*Gam+1);
    d2psi_2 = @(xi,h_e,Gam) h_e/4*(12*Gam-3*xi+1)./(12*Gam+1);
    d2psi_3 = @(xi,h_e,Gam) -3/2*xi./(12*Gam+1);
    d2psi_4 = @(xi,h_e,Gam) -h_e/4*(12*Gam+3*xi+1)./(12*Gam+1);
    d2psi = {d2psi_1; d2psi_2; d2psi_3; d2psi_4};
    d3psi_1 = @(xi,h_e,Gam) 3/2./(12*Gam+1)*ones(length(xi),1);
    d3psi_2 = @(xi,h_e,Gam) -3/4*h_e./(12*Gam+1)*ones(length(xi),1);
    d3psi_3 = @(xi,h_e,Gam) -3/2./(12*Gam+1)*ones(length(xi),1);
    d3psi_4 = @(xi,h_e,Gam) -3/4*h_e./(12*Gam+1)*ones(length(xi),1);
    d3psi = {d3psi_1; d3psi_2; d3psi_3; d3psi_4};
    % phi - interdependent generalized Hermite quadratic polynomials for bending rotations 
    phi_1 = @(xi,h_e,Gam) -3/2*(xi.^2-1)./(h_e*(1+12*Gam));
    phi_2 = @(xi,h_e,Gam) -1./(4*(1+12*Gam)).*(-3*(xi+1).^2+8*xi-4*(1+12*Gam)+24*Gam.*(xi+1)+8);
    phi_3 = @(xi,h_e,Gam) 3/2*(xi.^2-1)./(h_e*(1+12*Gam));
    phi_4 = @(xi,h_e,Gam) 1./(4*(1+12*Gam)).*(xi+1).*(3*xi-1+24*Gam);
    phi = {phi_1; phi_2; phi_3; phi_4};
    dphi_1 = @(xi,h_e,Gam) -3*xi./(h_e*(1+12*Gam));
    dphi_2 = @(xi,h_e,Gam) -1./(2*(1+12*Gam)).*(-3*xi+1+12*Gam);
    dphi_3 = @(xi,h_e,Gam) 3*xi./(h_e*(1+12*Gam));
    dphi_4 = @(xi,h_e,Gam) 1./(2*(1+12*Gam)).*(3*xi+1+12*Gam);
    dphi = {dphi_1; dphi_2; dphi_3; dphi_4};
    d2phi_1 = @(xi,h_e,Gam) -3./(h_e*(1+12*Gam));
    d2phi_2 = @(xi,h_e,Gam) 3./(2*(1+12*Gam));
    d2phi_3 = @(xi,h_e,Gam) 3./(h_e*(1+12*Gam));
    d2phi_4 = @(xi,h_e,Gam) 3./(2*(1+12*Gam));
    d2phi = {d2phi_1; d2phi_2; d2phi_3; d2phi_4};
else
    % zeta - Lagrange quadratic polynomials for axial displacement (and torsional rotations in case of no warping DOF)
    zeta_1 = @(xi) xi/2.*(xi-1);
    zeta_2 = @(xi) 1-xi.^2;
    zeta_3 = @(xi) xi/2.*(xi+1);
    zeta = {zeta_1; zeta_2; zeta_3};
    dzeta_1 = @(xi) xi-1/2;
    dzeta_2 = @(xi) -2*xi;
    dzeta_3 = @(xi) xi+1/2;
    dzeta = {dzeta_1; dzeta_2; dzeta_3};
    % psi - generalized Hermite quintic polynomials for transverse displacements (and torsional rotations in case of warping DOF)
    psi_1 = @(xi,h_e,Gam) -(xi.*(xi - 1).*(11520*Gam.^2.*xi + 5760*Gam.^2 - 144*Gam.*xi.^3 - 24*Gam.*xi.^2 + 456*Gam.*xi + 120*Gam - 3*xi.^3 - xi.^2 + 4*xi))./(4*(1+108*Gam+2880*Gam.^2));
    psi_2 = @(xi,h_e,Gam) (h_e*xi.*(xi.^2 - 1).*(- 92160*Gam.^3 + 1152*Gam.^2.*xi.^2 - 768*Gam.^2 - 24*Gam.*xi.^2 + 60*Gam.*xi + 24*Gam - xi.^2 + xi))./(8*(1+108*Gam+2880*Gam.^2));
    psi_3 = @(xi,h_e,Gam) -((xi.^2 - 1).*(- xi.^2 + 48*Gam + 1))./(48*Gam + 1);
    psi_4 = @(xi,h_e,Gam) (h_e*xi.*(xi.^2 - 1).*(960*Gam.^2 - 12*Gam.*xi.^2 + 108*Gam - xi.^2 + 1))./(120*Gam + 2);
    psi_5 = @(xi,h_e,Gam) (xi.*(xi + 1).*(11520*Gam.^2.*xi - 5760*Gam.^2 - 144*Gam.*xi.^3 + 24*Gam.*xi.^2 + 456*Gam.*xi - 120*Gam - 3*xi.^3 + xi.^2 + 4*xi))./(4*(1+108*Gam+2880*Gam.^2));
    psi_6 = @(xi,h_e,Gam) -(h_e*xi.*(xi.^2 - 1).*(92160*Gam.^3 - 1152*Gam.^2.*xi.^2 + 768*Gam.^2 + 24*Gam.*xi.^2 + 60*Gam.*xi - 24*Gam + xi.^2 + xi))./(8*(1+108*Gam+2880*Gam.^2));
    psi = {psi_1; psi_2; psi_3; psi_4; psi_5; psi_6};
    dpsi_1 = @(xi,h_e,Gam) (- 34560*Gam.^2.*xi.^2 + 11520*Gam.^2.*xi + 5760*Gam.^2 + 720*Gam.*xi.^4 - 480*Gam.*xi.^3 - 1440*Gam.*xi.^2 + 672*Gam.*xi + 120*Gam + 15*xi.^4 - 8*xi.^3 - 15*xi.^2 + 8*xi)./(4*(1+108*Gam+2880*Gam.^2));
    dpsi_2 = @(xi,h_e,Gam) -(h_e*(276480*Gam.^3.*xi.^2 - 92160*Gam.^3 - 5760*Gam.^2.*xi.^4 + 5760*Gam.^2.*xi.^2 - 768*Gam.^2 + 120*Gam.*xi.^4 - 240*Gam.*xi.^3 - 144*Gam.*xi.^2 + 120*Gam.*xi + 24*Gam + 5*xi.^4 - 4*xi.^3 - 3*xi.^2 + 2*xi))./(8*(1+108*Gam+2880*Gam.^2));
    dpsi_3 = @(xi,h_e,Gam) -(4*xi.*(- xi.^2 + 24*Gam + 1))./(48*Gam + 1);
    dpsi_4 = @(xi,h_e,Gam) -(h_e*(- 2880*Gam.^2.*xi.^2 + 960*Gam.^2 + 60*Gam.*xi.^4 - 360*Gam.*xi.^2 + 108*Gam + 5*xi.^4 - 6*xi.^2 + 1))./(2*(60*Gam + 1));
    dpsi_5 = @(xi,h_e,Gam) -(- 34560*Gam.^2.*xi.^2 - 11520*Gam.^2.*xi + 5760*Gam.^2 + 720*Gam.*xi.^4 + 480*Gam.*xi.^3 - 1440*Gam.*xi.^2 - 672*Gam.*xi + 120*Gam + 15*xi.^4 + 8*xi.^3 - 15*xi.^2 - 8*xi)./(4*(1+108*Gam+2880*Gam.^2));
    dpsi_6 = @(xi,h_e,Gam) -(h_e*(276480*Gam.^3.*xi.^2 - 92160*Gam.^3 - 5760*Gam.^2.*xi.^4 + 5760*Gam.^2.*xi.^2 - 768*Gam.^2 + 120*Gam.*xi.^4 + 240*Gam.*xi.^3 - 144*Gam.*xi.^2 - 120*Gam.*xi + 24*Gam + 5*xi.^4 + 4*xi.^3 - 3*xi.^2 - 2*xi))./(8*(1+108*Gam+2880*Gam.^2));
    dpsi = {dpsi_1; dpsi_2; dpsi_3; dpsi_4; dpsi_5; dpsi_6};
    d2psi_1 = @(xi,h_e,Gam) (Gam*336-xi*15-Gam.*xi*1440-Gam.*xi.^2*720-Gam.^2.*xi*34560+Gam.*xi.^3*1440+Gam.^2*5760-xi.^2*12+xi.^3*30+4)./(2*(1+108*Gam+2880*Gam.^2));
    d2psi_2 = @(xi,h_e,Gam) -h_e*(Gam.*120-xi.*6-Gam.*xi.*288-Gam.*xi.^2.*720+Gam.^2.*xi.*11520+Gam.*xi.^3.*480+Gam.^3.*xi.*552960-xi.^2.*12+xi.^3.*20-Gam.^2.*xi.^3.*23040+2)./(8*(1+108*Gam+2880*Gam.^2));
    d2psi_3 = @(xi,h_e,Gam) -4*(Gam.*24-xi.^2.*3+1)./(48*Gam+1);
    d2psi_4 = @(xi,h_e,Gam) (h_e*xi.*(Gam.*180-Gam.*xi.^2.*60+Gam.^2.*1440-xi.^2.*5+3).*2)./(60*Gam + 1);
    d2psi_5 = @(xi,h_e,Gam) (Gam.*336+xi.*15+Gam.*xi.*1440-Gam.*xi.^2.*720+Gam.^2.*xi.*34560-Gam.*xi.^3.*1440+Gam.^2.*5760-xi.^2.*12-xi.^3.*30+4)./(2*(1+108*Gam+2880*Gam.^2));
    d2psi_6 = @(xi,h_e,Gam) -h_e.*(Gam.*-120-xi.*6-Gam.*xi.*288+Gam.*xi.^2.*720+Gam.^2.*xi.*11520+Gam.*xi.^3.*480+Gam.^3.*xi.*552960+xi.^2.*12+xi.^3.*20-Gam.^2.*xi.^3.*23040-2)./(8*(1+108*Gam+2880*Gam.^2));
    d2psi = {d2psi_1; d2psi_2; d2psi_3; d2psi_4; d2psi_5; d2psi_6};
    d3psi_1 = @(xi,h_e,Gam) -(Gam.*1440+xi.*24+Gam.*xi.*1440-Gam.*xi.^2.*4320+Gam.^2.*34560-xi.^2.*90+15)./(2*(1+108*Gam+2880*Gam.^2));
    d3psi_2 = @(xi,h_e,Gam) (h_e.*(Gam.*288+xi.*24+Gam.*xi.*1440-Gam.*xi.^2.*1440-Gam.^2.*11520-Gam.^3.*552960-xi.^2.*60+Gam.^2.*xi.^2.*69120+6))./(8*(1+108*Gam+2880*Gam.^2));
    d3psi_3 = @(xi,h_e,Gam) 24*xi./(48*Gam+1);
    d3psi_4 = @(xi,h_e,Gam) 6*h_e.*(Gam.*60-Gam.*xi.^2.*60+Gam.^2.*480-xi.^2.*5+1)./(60*Gam+1);
    d3psi_5 = @(xi,h_e,Gam) -(Gam.*-1440+xi.*24+Gam.*xi.*1440+Gam.*xi.^2.*4320-Gam.^2.*34560+xi.^2.*90-15)./(2*(1+108*Gam+2880*Gam.^2));
    d3psi_6 = @(xi,h_e,Gam) -h_e.*(Gam.*-288+xi.*24+Gam.*xi.*1440+Gam.*xi.^2.*1440+Gam.^2.*11520+Gam.^3.*552960+xi.^2.*60-Gam.^2.*xi.^2.*69120-6)./(8*(1+108*Gam+2880*Gam.^2));
    d3psi = {d3psi_1; d3psi_2; d3psi_3; d3psi_4; d3psi_5; d3psi_6};
    % phi - interdependent generalized Hermite quartic polynomials for bending rotations 
    phi_1 = @(xi,h_e,Gam) (xi.*(xi - 1).*(xi + 1).*(480*Gam - 15*xi - 720*Gam.*xi + 8))./(2*h_e*(1+108*Gam+2880*Gam.^2));
    phi_2 = @(xi,h_e,Gam) -(xi.*(xi - 1).*(5760*Gam.^2.*xi.^2 + 5760*Gam.^2.*xi - 5760*Gam.^2 - 120*Gam.*xi.^2 + 120*Gam.*xi + 24*Gam - 5*xi.^2 - xi + 2))./(4*(1+108*Gam+2880*Gam.^2));
    phi_3 = @(xi,h_e,Gam) -(8*xi.*(xi.^2 - 1))./(h_e*(48*Gam + 1));
    phi_4 = @(xi,h_e,Gam) -((xi.^2 - 1).*(60*Gam - 60*Gam.*xi.^2 - 5*xi.^2 + 1))./(60*Gam + 1);
    phi_5 = @(xi,h_e,Gam) (xi.*(xi - 1).*(xi + 1).*(480*Gam + 15*xi + 720*Gam.*xi + 8))./(2*h_e*(1+108*Gam+2880*Gam.^2));
    phi_6 = @(xi,h_e,Gam) (xi.*(xi + 1).*(- 5760*Gam.^2.*xi.^2 + 5760*Gam.^2.*xi + 5760*Gam.^2 + 120*Gam.*xi.^2 + 120*Gam.*xi - 24*Gam + 5*xi.^2 - xi - 2))./(4*(1+108*Gam+2880*Gam.^2));
    phi = {phi_1; phi_2; phi_3; phi_4; phi_5; phi_6};
    dphi_1 = @(xi,h_e,Gam) -(240*Gam - 15*xi - 720*Gam.*xi - 720*Gam.*xi.^2 + 1440*Gam*xi.^3 - 12*xi.^2 + 30*xi.^3 + 4)./(h_e*(1+108*Gam+2880*Gam.^2));
    dphi_2 = @(xi,h_e,Gam) (- 11520*Gam.^2.*xi.^3 + 11520*Gam.^2.*xi - 2880*Gam.^2 + 240*Gam.*xi.^3 - 360*Gam.*xi.^2 + 96*Gam.*xi + 12*Gam + 10*xi.^3 - 6*xi.^2 - 3*xi + 1)./(2*(1+108*Gam+2880*Gam.^2));
    dphi_3 = @(xi,h_e,Gam) -(8*(3*xi.^2 - 1))./(h_e*(48*Gam + 1));
    dphi_4 = @(xi,h_e,Gam) -(4*xi.*(60*Gam - 60*Gam.*xi.^2 - 5*xi.^2 + 3))./(60*Gam + 1);
    dphi_5 = @(xi,h_e,Gam) -(240*Gam + 15*xi + 720*Gam.*xi - 720*Gam.*xi.^2 - 1440*Gam.*xi.^3 - 12*xi.^2 - 30*xi.^3 + 4)./(h_e*(1+108*Gam+2880*Gam.^2));
    dphi_6 = @(xi,h_e,Gam) (- 11520*Gam.^2.*xi.^3 + 11520*Gam.^2.*xi + 2880*Gam.^2 + 240*Gam.*xi.^3 + 360*Gam.*xi.^2 + 96*Gam.*xi - 12*Gam + 10*xi.^3 + 6*xi.^2 - 3*xi - 1)./(2*(1+108*Gam+2880*Gam.^2));
    dphi = {dphi_1; dphi_2; dphi_3; dphi_4; dphi_5; dphi_6};
    d2phi_1 = @(xi,h_e,Gam) (720*Gam+24*xi+1440*Gam.*xi-4320*Gam.*xi.^2-90*xi.^2+15)./(h_e*(1+108*Gam+2880*Gam.^2));
    d2phi_2 = @(xi,h_e,Gam) (96*Gam-12*xi-720*Gam.*xi+720*Gam.*xi.^2+11520*Gam.^2+30*xi.^2-34560*Gam.^2.*xi.^2-3)./(2*(1+108*Gam+2880*Gam.^2));
    d2phi_3 = @(xi,h_e,Gam) (-48*xi)./(h_e*(48*Gam + 1));
    d2phi_4 = @(xi,h_e,Gam) (-12*(20*Gam-60*Gam.*xi.^2-5*xi.^2+1))./(60*Gam + 1);
    d2phi_5 = @(xi,h_e,Gam) (-720*Gam+24*xi+1440*Gam.*xi+4320*Gam.*xi.^2+90*xi.^2-15)./(h_e*(1+108*Gam+2880*Gam.^2));
    d2phi_6 = @(xi,h_e,Gam) (96*Gam+12*xi+720*Gam.*xi+720*Gam.*xi.^2+11520*Gam.^2+30*xi.^2-34560*Gam.^2.*xi.^2-3)./(2*(1+108*Gam+2880*Gam.^2));
    d2phi = {d2phi_1; d2phi_2; d2phi_3; d2phi_4; d2phi_5; d2phi_6};
end

%% Interpolation function matrices
% Strain-displacement (B) matrices
B = @(i,xi,h_e,J,Gamx,Gamy,Gamz) [J^-1*dzeta{i}(xi) 0                                                              0                                                              0                                   0                                                           0                                                                                              0                                              
                                  0                 phi{2*i-1}(xi,h_e,Gamy(xi))+J^-1*dpsi{2*i-1}(xi,h_e,Gamy(xi))  0                                                              0                                   0                                                           phi{2*i}(xi,h_e,Gamy(xi))+J^-1*dpsi{2*i}(xi,h_e,Gamy(xi))                                      0                      
                                  0                 0                                                              phi{2*i-1}(xi,h_e,Gamz(xi))+J^-1*dpsi{2*i-1}(xi,h_e,Gamz(xi))  0                                   phi{2*i}(xi,h_e,Gamz(xi))+J^-1*dpsi{2*i}(xi,h_e,Gamz(xi))   0                                                                                              0                                            
                                  0                 0                                                              0                                                              J^-1*dpsi{2*i-1}(xi,h_e,Gamx(xi))   0                                                           0                                                                J^-1*dpsi{2*i}(xi,h_e,Gamx(xi))                        
                                  0                 0                                                              J^-1*dphi{2*i-1}(xi,h_e,Gamz(xi))                              0                                   J^-1*dphi{2*i}(xi,h_e,Gamz(xi))                             0                                                                                              0
                                  0                 J^-1*dphi{2*i-1}(xi,h_e,Gamy(xi))                              0                                                              0                                   0                                                           J^-1*dphi{2*i}(xi,h_e,Gamy(xi))                                                                0
                                  0                 0                                                              0                                                              J^-2*d2psi{2*i-1}(xi,h_e,Gamx(xi))  0                                                           0                                                               J^-2*d2psi{2*i}(xi,h_e,Gamx(xi))];
B_cf = @(i,xi,h_e,J,Gamy,Gamz) [0 psi{2*i-1}(xi,h_e,Gamy(xi))   0                           0 0                         psi{2*i}(xi,h_e,Gamy(xi)) 0
                                0 0                             psi{2*i-1}(xi,h_e,Gamz(xi)) 0 psi{2*i}(xi,h_e,Gamz(xi)) 0                         0];
% Translational inertia shape function matrix
B_Mt = @(i,xi,h_e,J,Gamx,Gamy,Gamz) [zeta{i}(xi)       0                            0                            0                             0                           0                                                     0                                              
                                     0                 psi{2*i-1}(xi,h_e,Gamy(xi))  0                            0                             0                           psi{2*i}(xi,h_e,Gamy(xi))                             0                      
                                     0                 0                            psi{2*i-1}(xi,h_e,Gamz(xi))  0                             psi{2*i}(xi,h_e,Gamz(xi))   0                                                     0                                            
                                     0                 0                            0                            psi{2*i-1}(xi,h_e,Gamx(xi))   0                           0                             psi{2*i}(xi,h_e,Gamx(xi))];
% Rotational inertia shape function matrix
B_Mr = @(i,xi,h_e,J,Gamx,Gamy,Gamz) [0 phi{2*i-1}(xi,h_e,Gamy(xi))  0                            0                             0                           phi{2*i}(xi,h_e,Gamy(xi))                             0                      
                                     0 0                            phi{2*i-1}(xi,h_e,Gamz(xi))  0                             phi{2*i}(xi,h_e,Gamz(xi))   0                                                     0                                            
                                     0 0                            0                            phi{2*i-1}(xi,h_e,Gamx(xi))   0                           0                             phi{2*i}(xi,h_e,Gamx(xi))];
% Shape functions for interpolation of concentraded sources
H_psi = @(i,xi,h_e,Gam) [psi{2*i-1}(xi,h_e,Gam(xi)); psi{2*i}(xi,h_e,Gam(xi))];
H_phi = @(i,xi,h_e,Gam) [phi{2*i-1}(xi,h_e,Gam(xi)); phi{2*i}(xi,h_e,Gam(xi))];
H_zeta = @(i,xi) zeta{i}(xi);

%% Add variables to FEMdata
FEMdata.funs.psi = psi;
FEMdata.funs.dpsi = dpsi;
FEMdata.funs.d2psi = d2psi;
FEMdata.funs.d3psi = d3psi;
FEMdata.funs.phi = phi;
FEMdata.funs.dphi = dphi;
FEMdata.funs.d2phi = d2phi;
FEMdata.funs.zeta = zeta;
FEMdata.funs.dzeta = dzeta;
FEMdata.funs.B = B;
FEMdata.funs.B_cf = B_cf;
FEMdata.funs.B_Mt = B_Mt;
FEMdata.funs.B_Mr = B_Mr;
FEMdata.funs.H_psi = H_psi;
FEMdata.funs.H_phi = H_phi;
FEMdata.funs.H_zeta = H_zeta;

