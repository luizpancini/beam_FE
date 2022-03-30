function [K_CS_g,K_CS_l,F_con_g,F_con_l] = LS3DB_get_con_loads_sources(FEMdata,edata,e)

% Unpack FEM and element data 
[warp_DOF,edof,enn,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,H_psi,H_phi,H_zeta] = unpack_FEMdata(FEMdata,'concentraded_loads');
[L0,Gamx,Gamy,Gamz,Px,Py,Pz,Mx,My,Mz,Pa,Pl,Pt,Ml,Mt,Tq,Bm,CSx,CSy,CSz,CSu,CSv,CSw,CSt,CSbl,CSbt,xi_apx,xi_apy,xi_apz,xi_amx,xi_amy,xi_amz,xi_apa,xi_apl,xi_apt,xi_atq,xi_aml,xi_amt,xi_abm,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt,xi_akbl,xi_akbt] = unpack_edata(edata,e,'concentraded_loads'); 

% DOF matrices
DOF_bendl = [DOF_v; DOF_phiz];
DOF_bendt = [DOF_w; DOF_phiy];
if warp_DOF 
    DOF_tors = [DOF_phix; DOF_dphix]; 
    H_tors = H_psi; H_mod_tors = [1; 1];
    H_tors_argin = {L0,Gamx};
else
    DOF_tors = DOF_phix; 
    H_tors = H_zeta; H_mod_tors = 1;
    H_tors_argin = {};
end

% Initialize outputs:
K_CS_g = zeros(edof);            % Matrix of global frame concentraded springs' stiffnesses
K_CS_l = zeros(edof);            % Matrix of local frame concentraded springs' stiffnesses
F_con_g = zeros(edof,1);         % Vector of concentraded forces applied in the global frame
F_con_l = zeros(edof,1);         % Vector of concentraded forces applied in the local frame

%% Add concentraded loads contributions
% Global frame 
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_apx,Px,H_zeta,1,DOF_u);                           % x-direction forces
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_apy,Py,H_psi,[1; -1],DOF_bendl,L0,Gamy);          % y-direction forces
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_apz,Pz,H_psi,[1; 1],DOF_bendt,L0,Gamz);           % z-direction forces
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amx,Mx,H_tors,H_mod_tors,DOF_tors,H_tors_argin);  % x-direction moments
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amy,My,H_phi,[1; 1],DOF_bendt,L0,Gamz);           % y-direction moments
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amz,Mz,H_phi,[-1; 1],DOF_bendl,L0,Gamy);          % z-direction moments

% Local frame 
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_apa,Pa,H_zeta,1,DOF_u);                             % Axial forces
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_apl,Pl,H_psi,[1; -1],DOF_bendl,L0,Gamy);            % Lateral forces
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_apt,Pt,H_psi,[1; 1],DOF_bendt,L0,Gamy);             % Transverse forces
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_atq,Tq,H_tors,H_mod_tors,DOF_tors,H_tors_argin);    % Twisting moments
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_aml,Ml,H_phi,[1; 1],DOF_bendt,L0,Gamz);             % Lateral moments
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_amt,Mt,H_phi,[-1; 1],DOF_bendl,L0,Gamy);            % Transverse moments
if warp_DOF
    F_con_l = apply_concentraded_loads(F_con_l,enn,xi_abm,Bm,H_phi,H_mod_tors,DOF_tors,H_tors_argin); % Bimoments
end  

%% Add concentraded sources contributions
% Global frame 
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSx,xi_akx,H_zeta,DOF_u);                 % x-direction translational springs
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSy,xi_aky,H_psi,DOF_bendl,L0,Gamy);      % y-direction translational springs
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSz,xi_akz,H_psi,DOF_bendt,L0,Gamz);      % z-direction translational springs

% Local frame 
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSu,xi_aku,H_zeta,DOF_u);                 % Axial springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSv,xi_akv,H_psi,DOF_bendl,L0,Gamy);      % Lateral springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSw,xi_akw,H_psi,DOF_bendt,L0,Gamz);      % Transverse springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSt,xi_akt,H_tors,DOF_tors,H_tors_argin); % Torsional springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSbl,xi_akbl,H_phi,DOF_bendl,L0,Gamy);    % Lateral rotation springs
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSbt,xi_akbt,H_phi,DOF_bendt,L0,Gamz);    % Transverse rotation springs

%% Nested functions  
    % Update the concentraded loads vector
    function F = apply_concentraded_loads(F,enn,xi_a,P,H_interp,H_mod,DOFs,varargin) 
        if ~isempty(varargin), if iscell(varargin{1}), varargin = varargin{1}; end; end
        indexat = @(expr, index) expr(index);   % To get specific vector function index
        for j=1:length(xi_a)                    % Loop over applied forces
            for n=1:enn                         % Loop over element's nodes
                % Distribute the generalized force across the nodes according to the interpolation function
                F(DOFs(:,n)) = F(DOFs(:,n)) + indexat(P(xi_a(j)),j) * H_interp(n,xi_a(j),varargin{:}) .* H_mod;
            end
        end
    end

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