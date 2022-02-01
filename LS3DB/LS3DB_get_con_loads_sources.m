function [K_CS_g,K_CS_l,F_con_g,F_con_l] = LS3DB_get_con_loads_sources(FEMdata,edata,e)

% Unpack FEM and element data 
[warp_DOF,edof,enn,DOF_u,DOF_phix,DOF_bendl,DOF_bendt,H_psi,H_phi,H_zeta] = unpack_FEMdata(FEMdata,'concentraded_loads');
[L0,Gamx,Gamy,Gamz,Px,Py,Pz,Mx,My,Mz,Pa,Pl,Pt,Ml,Mt,Tq,Bm,CSx,CSy,CSz,CSu,CSv,CSw,CSt,xi_apx,xi_apy,xi_apz,xi_amx,xi_amy,xi_amz,xi_apa,xi_apl,xi_apt,xi_atq,xi_aml,xi_amt,xi_abm,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt] = unpack_edata(edata,e,'concentraded_loads'); 

% Reshape DOF matrices
DOF_bendl = reshape(DOF_bendl,[2,enn]);
DOF_bendt = reshape(DOF_bendt,[2,enn]);
if warp_DOF, DOF_phix_warp = reshape(DOF_phix,[2,enn]); end

% Initialize outputs:
K_CS_g = zeros(edof);            % Matrix of global frame concentraded springs' stiffnesses
K_CS_l = zeros(edof);            % Matrix of local frame concentraded springs' stiffnesses
F_con_g = zeros(edof,1);         % Vector of concentraded forces applied in the global frame
F_con_l = zeros(edof,1);         % Vector of concentraded forces applied in the local frame

%% Add concentraded loads/sources contributions
% Global frame x-direction forces:
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_apx,Px,H_zeta,1,DOF_u);
% Global frame y-direction forces:
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_apy,Py,H_psi,[1; -1],DOF_bendl,L0,Gamy);
% Global frame z-direction forces:
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_apz,Pz,H_psi,[1; 1],DOF_bendt,L0,Gamz);
% Global frame x-direction moments:
if warp_DOF
    F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amx,Mx,H_psi,[1; 1],DOF_phix_warp,L0,Gamx);
else
    F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amx,Mx,H_zeta,1,DOF_phix);
end
% Global frame y-direction moments:
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amy,My,H_phi,[1; 1],DOF_bendt,L0,Gamz);
% Global frame z-direction moments:
F_con_g = apply_concentraded_loads(F_con_g,enn,xi_amz,Mz,H_phi,[-1; 1],DOF_bendl,L0,Gamy);
% Global frame x-direction springs:
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSx,xi_akx,H_zeta,DOF_u);
% Global frame y-direction springs:
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSy,xi_aky,H_psi,DOF_bendl,L0,Gamy);
% Global frame z-direction springs:
K_CS_g = apply_concentraded_sources(K_CS_g,enn,CSz,xi_akz,H_psi,DOF_bendt,L0,Gamz);
% Local frame axial forces:
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_apa,Pa,H_zeta,1,DOF_u);
% Local frame lateral forces:
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_apl,Pl,H_psi,[1; -1],DOF_bendl,L0,Gamy);
% Local frame transverse forces:
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_apt,Pt,H_psi,[1; 1],DOF_bendt,L0,Gamy);
% Local frame twisting torque:
if warp_DOF
    F_con_l = apply_concentraded_loads(F_con_l,enn,xi_atq,Tq,H_psi,[1; 1],DOF_phix_warp,L0,Gamx);
else
    F_con_l = apply_concentraded_loads(F_con_l,enn,xi_atq,Tq,H_zeta,1,DOF_phix);
end
% Local frame lateral moments:
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_aml,Ml,H_phi,[1; 1],DOF_bendt,L0,Gamz);
% Local frame transverse moments:
F_con_l = apply_concentraded_loads(F_con_l,enn,xi_amt,Mt,H_phi,[-1; 1],DOF_bendl,L0,Gamy);
% Local frame bimoments:
if warp_DOF 
    F_con_l = apply_concentraded_loads(F_con_l,enn,xi_abm,Bm,H_phi,[1; 1],DOF_phix_warp,L0,Gamx);
end
% Local frame axial springs:
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSu,xi_aku,H_zeta,DOF_u);
% Local frame lateral springs:
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSv,xi_akv,H_psi,DOF_bendl,L0,Gamy);
% Local frame transverse springs:
K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSw,xi_akw,H_psi,DOF_bendt,L0,Gamz);
% Local frame torsional springs:
if warp_DOF
    K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSt,xi_akt,H_psi,DOF_phix_warp,L0,Gamx);
else
    K_CS_l = apply_concentraded_sources(K_CS_l,enn,CSt,xi_akt,H_zeta,DOF_phix);
end

%% Nested functions
    
    % Update the concentraded loads vector
    function F = apply_concentraded_loads(F,enn,xi_a,P,H_interp,H_mod,DOFs,varargin)  
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