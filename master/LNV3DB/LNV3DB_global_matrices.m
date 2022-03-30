function FEMdata = LNV3DB_global_matrices(FEMdata,edata)

%% Global matrices
% Unpack FEM data
[Ne,Ndof,edof,enn,EDOFs,NGP,RI,B,B_cf,B_Mt,B_Mr,D,D_cf,D_Mt,D_Mr] = unpack_FEMdata(FEMdata,'global_matrices');

% Initialize outputs: global stiffness matrix and global mass/inertia matrix
K = zeros(Ndof); M = zeros(Ndof); 

% Loop over elements
for e=1:Ne         
    %% Element matrices
    % Unpack element data on parent (xi) coordinate system
    [e_dof_range,L0,beam,T0,Jac,E,G,rho,A,J,Gamma,Is,Iyy,Izz,Ksy,Ksz,e1,e2,Gamx,Gamy,Gamz,cfl,cft,cf_on_element] = unpack_edata(edata,e,'global_matrices');
    % Element stiffness and mass matrices
    K_e = zeros(edof); M_e = zeros(edof); K_DS_l = zeros(edof); 
    % Loop over nodes
    for ni=1:enn
        for nj=ni:enn
            % Element stiffness/inertia integrands 
            int_K_fun = @(xi) Jac * (B(ni,xi,L0,Jac,Gamx,Gamy,Gamz))' * D{beam}(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) * B(nj,xi,L0,Jac,Gamx,Gamy,Gamz);
            int_Mt_fun = @(xi) Jac * (B_Mt(ni,xi,L0,Jac,Gamx,Gamy,Gamz))' * D_Mt{beam}(xi,rho,A,Is,e1,e2) * B_Mt(nj,xi,L0,Jac,Gamx,Gamy,Gamz);
            % Integration over element domain
            K_e(EDOFs{ni},EDOFs{nj}) = gauss_legendre(int_K_fun,-1,1,NGP);
            if RI % Consider rotary inertia
                int_Mr_fun = @(xi) Jac * (B_Mr(ni,xi,L0,Jac,Gamx,Gamy,Gamz))' * D_Mr{beam}(xi,rho,Iyy,Izz,Gamma) * B_Mr(nj,xi,L0,Jac,Gamx,Gamy,Gamz);
                int_Mtr_fun = @(xi) int_Mt_fun(xi) + int_Mr_fun(xi);
                M_e(EDOFs{ni},EDOFs{nj}) = gauss_legendre(int_Mtr_fun,-1,1,NGP);
            else
                M_e(EDOFs{ni},EDOFs{nj}) = gauss_legendre(int_Mt_fun,-1,1,NGP);
            end
            % Contribution from distributed elastic foundations 
            if cf_on_element 
                int_Kcf_fun = @(xi) Jac * (B_cf(ni,xi,L0,Jac,Gamy,Gamz))' * D_cf{beam}(xi,cfl,cft) * B_cf(nj,xi,L0,Jac,Gamy,Gamz); 
                K_DS_l(EDOFs{ni},EDOFs{nj}) = gauss_legendre(int_Kcf_fun,-1,1,NGP+2);
            end
        end       
    end
    % Get global and local concentrated sources
    [K_CS_g,K_CS_l,M_CS_l] = LNV3DB_get_con_sources(FEMdata,edata,e);
    % Add stiffness contributions from nodes nj < ni
    K_e = mirror_matrix(K_e);
    M_e = mirror_matrix(M_e);
    K_DS_l = mirror_matrix(K_DS_l);
    K_CS_l = mirror_matrix(K_CS_l);
    K_CS_g = mirror_matrix(K_CS_g);
    M_CS_l = mirror_matrix(M_CS_l);
    % Add local frame sources
    K_e = K_e + K_DS_l + K_CS_l;
    M_e = M_e + M_CS_l;
    % Rotate from local to global frame
    K_e = T0*K_e*T0';
    M_e = T0*M_e*T0';
    % Add global frame sources
    K_e = K_e + K_CS_g;      
    
    %% Global assembly 
    % Add element contribution to global stiffness matrix 
    K(e_dof_range,e_dof_range) = K(e_dof_range,e_dof_range) + K_e;
    % Add element contribution to global mass/inertia matrix 
    M(e_dof_range,e_dof_range) = M(e_dof_range,e_dof_range) + M_e;   
end
% Add to FEM data
FEMdata.K = K;
FEMdata.M = M;

%% Nested functions
    % Mirror the upper triangular part of a matrix to make it symmetric
    function K = mirror_matrix(K)              
        % Keep only the upper entries
        K = triu(K);
        % Mirror
        K = K + K' - diag(diag(K));       
    end
end