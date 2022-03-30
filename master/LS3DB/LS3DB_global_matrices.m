function FEMdata = LS3DB_global_matrices(FEMdata,edata)

%% Global matrices
% Unpack FEM data
[Ne,Ndof,edof,enn,EDOFs,NGP,B,B_cf,H_dist,f_dist_l,f_dist_g,D,D_cf] = unpack_FEMdata(FEMdata,'global_matrices');

% Initialize outputs: global stiffness matrix and global force vector
K = zeros(Ndof); F = zeros(Ndof,1); 

% Loop over elements
for e=1:Ne         
    %% Element matrices
    % Unpack element data on parent (xi) coordinate system
    [e_dof_range,L0,beam,T0,Jac,E,G,A,J,Gamma,Iyy,Izz,Ksy,Ksz,Gamx,Gamy,Gamz,fx,qy,qz,mx,fa,ql,qt,tq,bm,cfl,cft,cf_on_element] = unpack_edata(edata,e,'global_matrices');
    % Element stiffness matrix, and distributed/concentraded, local/global force vectors
    K_e = zeros(edof); K_DS_l = zeros(edof); F_dist_l = zeros(edof,1); F_dist_g = zeros(edof,1); 
    % Loop over nodes
    for ni=1:enn
        for nj=ni:enn
            % Element stiffness integrand 
            int_K_fun = @(xi) Jac * (B(ni,xi,L0,Jac,Gamx,Gamy,Gamz))' * D{beam}(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) * B(nj,xi,L0,Jac,Gamx,Gamy,Gamz);
            % Integration over element domain
            K_e(EDOFs{ni},EDOFs{nj}) = gauss_legendre(int_K_fun,-1,1,NGP);
            % Contribution from distributed elastic foundations 
            if cf_on_element 
                int_Kcf_fun = @(xi) Jac * (B_cf(ni,xi,L0,Jac,Gamy,Gamz))' * D_cf{beam}(xi,cfl,cft) * B_cf(nj,xi,L0,Jac,Gamy,Gamz); 
                K_DS_l(EDOFs{ni},EDOFs{nj}) = gauss_legendre(int_Kcf_fun,-1,1,NGP+2);
            end
        end       
        % Distributed forces integrands
        int_F_dist_l_fun = @(xi) Jac * f_dist_l(xi,fa,ql,qt,tq,bm) .* H_dist(ni,xi,L0,Gamx,Gamy,Gamz);  
        int_F_dist_g_fun = @(xi) Jac * f_dist_g(xi,fx,qy,qz,mx) .* H_dist(ni,xi,L0,Gamx,Gamy,Gamz);
        % Integration over element domain
        F_dist_l(EDOFs{ni}) = gauss_legendre(int_F_dist_l_fun,-1,1,NGP);
        F_dist_g(EDOFs{ni}) = gauss_legendre(int_F_dist_g_fun,-1,1,NGP);
    end
    % Get global and local concentrated loads/sources
    [K_CS_g,K_CS_l,F_con_g,F_con_l] = LS3DB_get_con_loads_sources(FEMdata,edata,e);
    % Add stiffness contributions from nodes nj < ni
    K_e = mirror_matrix(K_e);
    K_DS_l = mirror_matrix(K_DS_l);
    K_CS_l = mirror_matrix(K_CS_l);
    K_CS_g = mirror_matrix(K_CS_g);
    % Add local frame sources
    F_e = F_dist_l + F_con_l;
    K_e = K_e + K_DS_l + K_CS_l;
    % Rotate from local to global frame
    F_e = T0*F_e;
    K_e = T0*K_e*T0';
    % Add global frame sources
    F_e = F_e + F_con_g + F_dist_g; 
    K_e = K_e + K_CS_g;      

    %% Global assembly 
    % Add element contribution to global stiffness matrix 
    K(e_dof_range,e_dof_range) = K(e_dof_range,e_dof_range) + K_e;
    % Add element contribution to global external forces vector 
    F(e_dof_range) = F(e_dof_range) + F_e;
    
end

% Add to FEM data
FEMdata.K = K;
FEMdata.F = F;

%% Nested functions

    % Mirror the upper triangular part of a matrix to make it symmetric
    function K = mirror_matrix(K)              
        % Keep only the upper entries
        K = triu(K);
        % Mirror
        K = K + K' - diag(diag(K));       
    end

end