function SOLdata = LNV3DB_get_outputs(edata,FEMdata,SOLdata)

% Unpack FEM data
[N_modes,Ne,enn,n_div,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,psi,dpsi,phi,zeta] = unpack_FEMdata(FEMdata,'outputs');

% Initialize outputs: generalized local displacements, and global mode shapes
u_interp = cell(Ne,N_modes); v_interp = cell(Ne,N_modes); w_interp = cell(Ne,N_modes); phix_interp = cell(Ne,N_modes); phiy_interp = cell(Ne,N_modes); phiz_interp = cell(Ne,N_modes); dphix_interp = cell(Ne,N_modes); 
u_nodes = cell(Ne,N_modes); v_nodes = cell(Ne,N_modes); w_nodes = cell(Ne,N_modes); phix_nodes = cell(Ne,N_modes); phiy_nodes = cell(Ne,N_modes); phiz_nodes = cell(Ne,N_modes); dphix_nodes = cell(Ne,N_modes); 
x_def_nodes = cell(Ne,N_modes); y_def_nodes = cell(Ne,N_modes); z_def_nodes = cell(Ne,N_modes); x_def_interp = cell(Ne,N_modes); y_def_interp = cell(Ne,N_modes); z_def_interp = cell(Ne,N_modes);

% Loop over modes
for m=1:N_modes
    % Loop over elements
    for e=1:Ne        
        %% Initialize element interpolation functions
        u_interp_fun = @(x) 0; v_interp_fun = @(x) 0; w_interp_fun = @(x) 0; phix_interp_fun = @(x) 0; phiy_interp_fun = @(x) 0; phiz_interp_fun = @(x) 0; dphix_interp_fun = @(x) 0;
        
        %% Generalized displacements - as functions of the local beam coordinate x
        % Unpack element data in the local coordinate x
        [e_dof_range,L0,x1,R0,T0,Jac,x_vec_nodes,x_vec_interp,x0,y0,z0,Gamx,Gamy,Gamz] = unpack_edata(edata,e,'outputs');
        % Define parent coordinate xi as a function of x
        xi = @(x) Jac^-1*(x-x1)-1;
        % Element generalized displacements in local frame
        u_local = T0'*SOLdata.W(e_dof_range,m);
        % Loop over nodes - get functions of x
        for n=1:enn
            % Local displacements
            u_interp_fun = @(x) u_interp_fun(x) + u_local(DOF_u(n))*zeta{n}(xi(x));
            v_interp_fun = @(x) v_interp_fun(x) + u_local(DOF_v(n))*psi{2*n-1}(xi(x),L0,Gamy(x)) + u_local(DOF_phiz(n))*psi{2*n}(xi(x),L0,Gamy(x));
            w_interp_fun = @(x) w_interp_fun(x) + u_local(DOF_w(n))*psi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*psi{2*n}(xi(x),L0,Gamz(x));
            phix_interp_fun = @(x) phix_interp_fun(x) + u_local(DOF_phix(n))*psi{2*n-1}(xi(x),L0,Gamx(x)) + u_local(DOF_dphix(n))*psi{2*n}(xi(x),L0,Gamx(x));
            phiy_interp_fun = @(x) phiy_interp_fun(x) + (u_local(DOF_w(n))*phi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*phi{2*n}(xi(x),L0,Gamz(x)));
            phiz_interp_fun = @(x) phiz_interp_fun(x) + (u_local(DOF_v(n))*phi{2*n-1}(xi(x),L0,Gamy(x)) + u_local(DOF_phiz(n))*phi{2*n}(xi(x),L0,Gamy(x)));
            dphix_interp_fun = @(x) dphix_interp_fun(x) + Jac^-1*(u_local(DOF_phix(n))*dpsi{2*n-1}(xi(x),L0,Gamx(x)) + u_local(DOF_dphix(n))*dpsi{2*n}(xi(x),L0,Gamx(x)));
        end
        
        %% Generalized displacements - local numerical values
        % Nodal values
        u_nodes{e,m} = u_interp_fun(x_vec_nodes);
        v_nodes{e,m} = v_interp_fun(x_vec_nodes);
        w_nodes{e,m} = w_interp_fun(x_vec_nodes);
        phix_nodes{e,m} = phix_interp_fun(x_vec_nodes);
        phiy_nodes{e,m} = phiy_interp_fun(x_vec_nodes);
        phiz_nodes{e,m} = phiz_interp_fun(x_vec_nodes);
        dphix_nodes{e,m} = dphix_interp_fun(x_vec_nodes);
        % Interpolated values
        u_interp{e,m} = u_interp_fun(x_vec_interp);
        v_interp{e,m} = v_interp_fun(x_vec_interp);
        w_interp{e,m} = w_interp_fun(x_vec_interp);
        phix_interp{e,m} = phix_interp_fun(x_vec_interp);
        phiy_interp{e,m} = phiy_interp_fun(x_vec_interp);
        phiz_interp{e,m} = phiz_interp_fun(x_vec_interp);
        dphix_interp{e,m} = dphix_interp_fun(x_vec_interp);     
        
        %% Generalized displacements - global numerical values
        % Nodal values
        vec_def_nodal = zeros(3,enn); vec_defr_nodal = vec_def_nodal;
        x_def_nodal_l = x_vec_nodes + u_nodes{e,m};
        y_def_nodal_l = v_nodes{e,m};
        z_def_nodal_l = w_nodes{e,m};
        for n=1:enn
            vec_def_nodal(:,n) = [x_def_nodal_l(n); y_def_nodal_l(n); z_def_nodal_l(n)];
            vec_defr_nodal(:,n) = R0*vec_def_nodal(:,n);
            x_def_nodal_l(n) = vec_defr_nodal(1,n);
            y_def_nodal_l(n) = vec_defr_nodal(2,n);
            z_def_nodal_l(n) = vec_defr_nodal(3,n);
        end
        x_def_nodes{e,m} = x_def_nodal_l + x0;
        y_def_nodes{e,m} = y_def_nodal_l + y0;
        z_def_nodes{e,m} = z_def_nodal_l + z0;
        % Interpolated values
        vec_def_interp = zeros(3,n_div); vec_defr_cont = vec_def_interp;
        x_def_interp_l = x_vec_interp + u_interp{e,m};
        y_def_interp_l = v_interp{e,m};
        z_def_interp_l = w_interp{e,m};
        for n=1:n_div
            vec_def_interp(:,n) = [x_def_interp_l(n); y_def_interp_l(n); z_def_interp_l(n)];
            vec_defr_cont(:,n) = R0*vec_def_interp(:,n);
            x_def_interp_l(n) = vec_defr_cont(1,n);
            y_def_interp_l(n) = vec_defr_cont(2,n);
            z_def_interp_l(n) = vec_defr_cont(3,n);
        end
        x_def_interp{e,m} = x_def_interp_l + x0;
        y_def_interp{e,m} = y_def_interp_l + y0;
        z_def_interp{e,m} = z_def_interp_l + z0;
    end
end

% Add to solution structure
SOLdata.u_nodes = u_nodes;
SOLdata.v_nodes = v_nodes;
SOLdata.w_nodes = w_nodes;
SOLdata.phix_nodes = phix_nodes;
SOLdata.phiy_nodes = phiy_nodes;
SOLdata.phiz_nodes = phiz_nodes;
SOLdata.dphix_nodes = dphix_nodes;
SOLdata.u_interp = u_interp;
SOLdata.v_interp = v_interp;
SOLdata.w_interp = w_interp;
SOLdata.phix_interp = phix_interp;
SOLdata.phiy_interp = phiy_interp;
SOLdata.phiz_interp = phiz_interp;
SOLdata.dphix_interp = dphix_interp;
SOLdata.x_def_nodes = x_def_nodes;
SOLdata.y_def_nodes = y_def_nodes;
SOLdata.z_def_nodes = z_def_nodes;
SOLdata.x_def_interp = x_def_interp;
SOLdata.y_def_interp = y_def_interp;
SOLdata.z_def_interp = z_def_interp;
