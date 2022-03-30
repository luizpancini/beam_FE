function [K_CS_g,K_CS_l,M_CS_l] = LNV3DB_get_con_sources(FEMdata,edata,e)

% Unpack FEM and element data 
[edof,enn,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,H_psi,H_phi,H_zeta] = unpack_FEMdata(FEMdata,'concentraded_sources');
[L0,Gamx,Gamy,Gamz,CSx,CSy,CSz,CSu,CSv,CSw,CSt,CSbl,CSbt,Ma,MaRix,MaRiy,MaRiz,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt,xi_akbl,xi_akbt,xi_ama] = unpack_edata(edata,e,'concentraded_sources'); 

% DOF matrices
DOF_bendl = [DOF_v; DOF_phiz];
DOF_bendt = [DOF_w; DOF_phiy];
DOF_tors = [DOF_phix; DOF_dphix];

% Initialize outputs:
K_CS_g = zeros(edof);   % Matrix of global frame concentraded springs' stiffnesses
K_CS_l = zeros(edof);   % Matrix of local frame concentraded springs' stiffnesses
M_CS_l = zeros(edof);   % Matrix of local frame concentraded masses' inertias

%% Add concentraded stiffness sources contributions
% Global frame
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSx,xi_akx,H_zeta,DOF_u);                % x-direction springs
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSy,xi_aky,H_psi,DOF_bendl,L0,Gamy);     % y-direction springs
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSz,xi_akz,H_psi,DOF_bendt,L0,Gamz);     % z-direction springs

% Local frame 
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSu,xi_aku,H_zeta,DOF_u);                % Axial springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSv,xi_akv,H_psi,DOF_bendl,L0,Gamy);     % Lateral springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSw,xi_akw,H_psi,DOF_bendt,L0,Gamz);     % Transverse springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSt,xi_akt,H_psi,DOF_tors,L0,Gamx);      % Torsional springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSbl,xi_akbl,H_phi,DOF_bendl,L0,Gamy);   % Lateral rotation springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSbt,xi_akbt,H_phi,DOF_bendt,L0,Gamz);   % Transverse rotation springs

%% Add concentraded inertia sources contributions
% Local frame 
M_CS_l = apply_concentraded_sources(M_CS_l,enn,Ma,xi_ama,H_zeta,DOF_u);                 % Mass' axial translational inertia
M_CS_l = apply_concentraded_sources(M_CS_l,enn,Ma,xi_ama,H_psi,DOF_bendl,L0,Gamy);      % Mass' lateral translational inertia
M_CS_l = apply_concentraded_sources(M_CS_l,enn,Ma,xi_ama,H_psi,DOF_bendt,L0,Gamz);      % Mass' transverse translational inertia
M_CS_l = apply_concentraded_sources(M_CS_l,enn,MaRix,xi_ama,H_psi,DOF_tors,L0,Gamx);    % Mass' torsional "translational" inertia
M_CS_l = apply_concentraded_sources(M_CS_l,enn,MaRiz,xi_ama,H_phi,DOF_bendl,L0,Gamy);   % Mass' lateral bending rotary inertia inertia
M_CS_l = apply_concentraded_sources(M_CS_l,enn,MaRiy,xi_ama,H_phi,DOF_bendt,L0,Gamz);   % Mass' transverse bending rotary inertia

%% Nested functions    
    % Update the concentraded sources matrix
    function K = apply_concentraded_sources(K,enn,K_CS,xi_a,H_interp,DOFs,varargin)
        indexat = @(expr, index) expr(index);   % To get specific vector function index
        for j=1:length(xi_a)                    % Loop over applied sources
            for ni=1:enn                        % Loop over element's nodes
                for nj=ni:enn
                    % Distribute the sources across the nodes according to the interpolation function
                    K(DOFs(:,ni),DOFs(:,nj)) = K(DOFs(:,ni),DOFs(:,nj)) + indexat(K_CS(xi_a(j)),j) * H_interp(ni,xi_a(j),varargin{:}) * (H_interp(nj,xi_a(j),varargin{:}))';
                end
            end
        end        
    end
    
end