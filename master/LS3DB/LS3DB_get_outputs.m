function SOLdata = LS3DB_get_outputs(edata,FEMdata,SOLdata)

% Unpack FEM data
[beam_theory,warp_DOF,constitutive_model,scale,Ne,enn,n_div,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,psi,dpsi,d2psi,d3psi,phi,dphi,d2phi,zeta,dzeta] = unpack_FEMdata(FEMdata,'outputs');

% Initialize outputs: generalized displacements and internal forces, and deformed structure position
u_local_interp = cell(Ne,1); v_local_interp = cell(Ne,1); w_local_interp = cell(Ne,1); phix_local_interp = cell(Ne,1); phiy_local_interp = cell(Ne,1); phiz_local_interp = cell(Ne,1); dphix_local_interp = cell(Ne,1); Nx_interp = cell(Ne,1); My_interp = cell(Ne,1); Mz_interp = cell(Ne,1); Vy_interp = cell(Ne,1); Vz_interp = cell(Ne,1); Tq_interp = cell(Ne,1); Bm_interp = cell(Ne,1);
u_local_nodes = cell(Ne,1); v_local_nodes = cell(Ne,1); w_local_nodes = cell(Ne,1); phix_local_nodes = cell(Ne,1); phiy_local_nodes = cell(Ne,1); phiz_local_nodes = cell(Ne,1); dphix_local_nodes = cell(Ne,1); Nx_nodes = cell(Ne,1); My_nodes = cell(Ne,1); Mz_nodes = cell(Ne,1); Vy_nodes = cell(Ne,1); Vz_nodes = cell(Ne,1); Tq_nodes = cell(Ne,1); Bm_nodes = cell(Ne,1);
x_def_nodes = cell(Ne,1); y_def_nodes = cell(Ne,1); z_def_nodes = cell(Ne,1); x_def_interp = cell(Ne,1); y_def_interp = cell(Ne,1); z_def_interp = cell(Ne,1);

% Loop over elements
for e=1:Ne
    
    %% Initialize element interpolation functions 
    u_interp_fun = @(x) 0; v_interp_fun = @(x) 0; w_interp_fun = @(x) 0; 
    du_interp_fun = @(x) 0; dv_interp_fun = @(x) 0; dw_interp_fun = @(x) 0;
    phix_interp_fun = @(x) 0; phiy_interp_fun = @(x) 0; phiz_interp_fun = @(x) 0; 
    dphix_interp_fun = @(x) 0; dphiy_interp_fun = @(x) 0; dphiz_interp_fun = @(x) 0;
    d2phix_interp_fun = @(x) 0; d3phix_interp_fun = @(x) 0; 
    d2v_interp_fun = @(x) 0; d3v_interp_fun = @(x) 0; 
    d2w_interp_fun = @(x) 0; d3w_interp_fun = @(x) 0;
    
    %% Displacements and internal forces - as functions of the local beam coordinate x
    % Unpack element data in the local coordinate x    
    [e_dof_range,L0,x1,beam,R0,T0,Jac,x_vec_nodes,x_vec_interp,x0,y0,z0,E,G,A,J,Gamma,Iyy,Izz,Ksy,Ksz,Gamx,Gamy,Gamz] = unpack_edata(edata,e,'outputs');  
    % Define parent coordinate xi as a function of x
    xi = @(x) Jac^-1*(x-x1)-1;
    % Element generalized displacements in local frame
    u_local = T0'*SOLdata.u_global(e_dof_range);
    SOLdata.u_local{e} = u_local;
    % Loop over nodes - get functions of x
    for n=1:enn        
        % Local displacements
        u_interp_fun = @(x) u_interp_fun(x) + u_local(DOF_u(n))*zeta{n}(xi(x));
        v_interp_fun = @(x) v_interp_fun(x) + u_local(DOF_v(n))*psi{2*n-1}(xi(x),L0,Gamy(x)) + (-1)*u_local(DOF_phiz(n))*psi{2*n}(xi(x),L0,Gamy(x));
        w_interp_fun = @(x) w_interp_fun(x) + u_local(DOF_w(n))*psi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*psi{2*n}(xi(x),L0,Gamz(x));
        if warp_DOF
            phix_interp_fun = @(x) phix_interp_fun(x) + u_local(DOF_phix(n))*psi{2*n-1}(xi(x),L0,Gamx(x)) + u_local(DOF_dphix(n))*psi{2*n}(xi(x),L0,Gamx(x));
        else
            phix_interp_fun = @(x) phix_interp_fun(x) + u_local(DOF_phix(n))*zeta{n}(xi(x));
        end
        phiy_interp_fun = @(x) phiy_interp_fun(x) + (u_local(DOF_w(n))*phi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*phi{2*n}(xi(x),L0,Gamz(x)));
        phiz_interp_fun = @(x) phiz_interp_fun(x) - (u_local(DOF_v(n))*phi{2*n-1}(xi(x),L0,Gamy(x)) + (-1)*u_local(DOF_phiz(n))*phi{2*n}(xi(x),L0,Gamy(x)));
        % Local derivatives
        du_interp_fun = @(x) du_interp_fun(x) + Jac^-1*(u_local(DOF_u(n))*dzeta{n}(xi(x)));
        dv_interp_fun = @(x) dv_interp_fun(x) + Jac^-1*(u_local(DOF_v(n))*dpsi{2*n-1}(xi(x),L0,Gamy(x)) + (-1)*u_local(DOF_phiz(n))*dpsi{2*n}(xi(x),L0,Gamy(x)));
        dw_interp_fun = @(x) dw_interp_fun(x) + Jac^-1*(u_local(DOF_w(n))*dpsi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*dpsi{2*n}(xi(x),L0,Gamz(x)));
        d2v_interp_fun = @(x) d2v_interp_fun(x) + Jac^-2*(u_local(DOF_v(n))*d2psi{2*n-1}(xi(x),L0,Gamy(x)) + (-1)*u_local(DOF_phiz(n))*d2psi{2*n}(xi(x),L0,Gamy(x)));
        d2w_interp_fun = @(x) d2w_interp_fun(x) + Jac^-2*(u_local(DOF_w(n))*d2psi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*d2psi{2*n}(xi(x),L0,Gamz(x)));
        d3v_interp_fun = @(x) d3v_interp_fun(x) + Jac^-3*(u_local(DOF_v(n))*d3psi{2*n-1}(xi(x),L0,Gamy(x)) + (-1)*u_local(DOF_phiz(n))*d3psi{2*n}(xi(x),L0,Gamy(x)));
        d3w_interp_fun = @(x) d3w_interp_fun(x) + Jac^-3*(u_local(DOF_w(n))*d3psi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*d3psi{2*n}(xi(x),L0,Gamz(x)));
        dphiy_interp_fun = @(x) dphiy_interp_fun(x) - Jac^-1*(u_local(DOF_w(n))*dphi{2*n-1}(xi(x),L0,Gamz(x)) + u_local(DOF_phiy(n))*dphi{2*n}(xi(x),L0,Gamz(x)));
        dphiz_interp_fun = @(x) dphiz_interp_fun(x) - Jac^-1*(u_local(DOF_v(n))*dphi{2*n-1}(xi(x),L0,Gamy(x)) + (-1)*u_local(DOF_phiz(n))*dphi{2*n}(xi(x),L0,Gamy(x)));
        if warp_DOF
            dphix_interp_fun = @(x) dphix_interp_fun(x) - (u_local(DOF_phix(n))*phi{2*n-1}(xi(x),L0,Gamx(x)) + u_local(DOF_dphix(n))*phi{2*n}(xi(x),L0,Gamx(x)));
            d2phix_interp_fun = @(x) d2phix_interp_fun(x) - Jac^-1*(u_local(DOF_phix(n))*dphi{2*n-1}(xi(x),L0,Gamx(x)) + u_local(DOF_dphix(n))*dphi{2*n}(xi(x),L0,Gamx(x)));
            d3phix_interp_fun = @(x) d3phix_interp_fun(x) - Jac^-2*(u_local(DOF_phix(n))*d2phi{2*n-1}(xi(x),L0,Gamx(x)) + u_local(DOF_dphix(n))*d2phi{2*n}(xi(x),L0,Gamx(x)));
        else
            dphix_interp_fun = @(x) dphix_interp_fun(x) + Jac^-1*(u_local(DOF_phix(n))*dzeta{n}(xi(x)));
        end
        % Local internal loads
        if constitutive_model(beam) == "iso"
            % Nx = E*A*du/dx, Vy = G*A*Ksy*(phiz+dv/dx) or E*Izz*d3v/dx3, Vz = G*A*Ksz*(phiy+dw/dx) or E*Iyy*d3w/dx3, Tq = G*J*dphi/dx - E*Gamma*d3phix/dx3, My = E*Iyy*dphiy/dx or E*Iyy*d2w/dx2 , Mz = E*Izz*dphiz/dx or or E*Izz*d2v/dx2, Bm = E*Gamma*d2phix/dx2
            Nx_interp_fun = @(x) E(x).*A(x).*du_interp_fun(x);
            if beam_theory == "EB"
                Vy_interp_fun = @(x) E(x).*Izz(x).*d3v_interp_fun(x);
                Vz_interp_fun = @(x) E(x).*Iyy(x).*d3w_interp_fun(x);
                My_interp_fun = @(x) E(x).*Iyy(x).*d2w_interp_fun(x);
                Mz_interp_fun = @(x) E(x).*Izz(x).*d2v_interp_fun(x);
            elseif beam_theory == "T"
                Vy_interp_fun = @(x) Ksy(x).*G(x).*A(x).*(phiz_interp_fun(x) - dv_interp_fun(x));
                Vz_interp_fun = @(x) -Ksz(x).*G(x).*A(x).*(phiy_interp_fun(x) + dw_interp_fun(x));
                My_interp_fun = @(x) E(x).*Iyy(x).*dphiy_interp_fun(x);
                Mz_interp_fun = @(x) E(x).*Izz(x).*dphiz_interp_fun(x);
            end
            Tq_interp_fun = @(x) -(G(x).*J(x).*dphix_interp_fun(x) - E(x).*Gamma(x).*d3phix_interp_fun(x));
            Bm_interp_fun = @(x) -E(x).*Gamma(x).*d2phix_interp_fun(x);
        end
    end
    
    %% Displacements and internal forces - numerical values
    % Nodal values   
    u_local_nodes{e} = u_interp_fun(x_vec_nodes);
    v_local_nodes{e} = v_interp_fun(x_vec_nodes);
    w_local_nodes{e} = w_interp_fun(x_vec_nodes);
    phix_local_nodes{e} = phix_interp_fun(x_vec_nodes);
    phiy_local_nodes{e} = phiy_interp_fun(x_vec_nodes);
    phiz_local_nodes{e} = phiz_interp_fun(x_vec_nodes);
    dphix_local_nodes{e} = dphix_interp_fun(x_vec_nodes);
    Nx_nodes{e} = Nx_interp_fun(x_vec_nodes);
    Vy_nodes{e} = Vy_interp_fun(x_vec_nodes);
    Vz_nodes{e} = Vz_interp_fun(x_vec_nodes);
    Tq_nodes{e} = Tq_interp_fun(x_vec_nodes);
    My_nodes{e} = My_interp_fun(x_vec_nodes);
    Mz_nodes{e} = Mz_interp_fun(x_vec_nodes);
    Bm_nodes{e} = Bm_interp_fun(x_vec_nodes);
    % Interpolated values
    u_local_interp{e} = u_interp_fun(x_vec_interp);
    v_local_interp{e} = v_interp_fun(x_vec_interp);
    w_local_interp{e} = w_interp_fun(x_vec_interp);
    phix_local_interp{e} = phix_interp_fun(x_vec_interp);
    phiy_local_interp{e} = phiy_interp_fun(x_vec_interp);
    phiz_local_interp{e} = phiz_interp_fun(x_vec_interp);
    dphix_local_interp{e} = dphix_interp_fun(x_vec_interp);
    Nx_interp{e} = Nx_interp_fun(x_vec_interp);
    Vy_interp{e} = Vy_interp_fun(x_vec_interp);
    Vz_interp{e} = Vz_interp_fun(x_vec_interp);
    Tq_interp{e} = Tq_interp_fun(x_vec_interp);
    My_interp{e} = My_interp_fun(x_vec_interp);
    Mz_interp{e} = Mz_interp_fun(x_vec_interp);
    Bm_interp{e} = Bm_interp_fun(x_vec_interp);
    
    %% Deformed structure data
    % Nodal values    
    vec_def_nodal = zeros(3,enn); vec_defr_nodal = vec_def_nodal;    
    x_def_nodal_l = x_vec_nodes + scale*u_local_nodes{e};
    y_def_nodal_l = scale*v_local_nodes{e};
    z_def_nodal_l = scale*w_local_nodes{e};
    for n=1:enn
        vec_def_nodal(:,n) = [x_def_nodal_l(n); y_def_nodal_l(n); z_def_nodal_l(n)];
        vec_defr_nodal(:,n) = R0*vec_def_nodal(:,n);
        x_def_nodal_l(n) = vec_defr_nodal(1,n); 
        y_def_nodal_l(n) = vec_defr_nodal(2,n); 
        z_def_nodal_l(n) = vec_defr_nodal(3,n);
    end
    x_def_nodes{e} = x_def_nodal_l + x0;
    y_def_nodes{e} = y_def_nodal_l + y0; 
    z_def_nodes{e} = z_def_nodal_l + z0;
    % Interpolated values
    vec_def_interp = zeros(3,n_div); vec_defr_cont = vec_def_interp; 
    x_def_interp_l = x_vec_interp + scale*u_local_interp{e};
    y_def_interp_l = scale*v_local_interp{e};
    z_def_interp_l = scale*w_local_interp{e};
    for n=1:n_div
        vec_def_interp(:,n) = [x_def_interp_l(n); y_def_interp_l(n); z_def_interp_l(n)];
        vec_defr_cont(:,n) = R0*vec_def_interp(:,n);
        x_def_interp_l(n) = vec_defr_cont(1,n);
        y_def_interp_l(n) = vec_defr_cont(2,n); 
        z_def_interp_l(n) = vec_defr_cont(3,n);
    end
    x_def_interp{e} = x_def_interp_l + x0; 
    y_def_interp{e} = y_def_interp_l + y0; 
    z_def_interp{e} = z_def_interp_l + z0; 
        
end

% Add to solution structure
SOLdata.u_local_nodes = u_local_nodes;
SOLdata.v_local_nodes = v_local_nodes;
SOLdata.w_local_nodes = w_local_nodes;
SOLdata.phix_local_nodes = phix_local_nodes;
SOLdata.phiy_local_nodes = phiy_local_nodes;
SOLdata.phiz_local_nodes = phiz_local_nodes;
SOLdata.dphix_local_nodes = dphix_local_nodes;
SOLdata.Nx_nodes = Nx_nodes;
SOLdata.Vy_nodes = Vy_nodes;
SOLdata.Vz_nodes = Vz_nodes;
SOLdata.Tq_nodes = Tq_nodes;
SOLdata.My_nodes = My_nodes;
SOLdata.Mz_nodes = Mz_nodes;
SOLdata.Bm_nodes = Bm_nodes;
SOLdata.u_local_interp = u_local_interp;
SOLdata.v_local_interp = v_local_interp;
SOLdata.w_local_interp = w_local_interp;
SOLdata.phix_local_interp = phix_local_interp;
SOLdata.phiy_local_interp = phiy_local_interp;
SOLdata.phiz_local_interp = phiz_local_interp;
SOLdata.dphix_local_interp = dphix_local_interp;
SOLdata.Nx_interp = Nx_interp;
SOLdata.Vy_interp = Vy_interp;
SOLdata.Vz_interp = Vz_interp;
SOLdata.Tq_interp = Tq_interp;
SOLdata.My_interp = My_interp;
SOLdata.Mz_interp = Mz_interp;
SOLdata.Bm_interp = Bm_interp;
SOLdata.x_def_nodes = x_def_nodes;
SOLdata.y_def_nodes = y_def_nodes;
SOLdata.z_def_nodes = z_def_nodes;
SOLdata.x_def_interp = x_def_interp;
SOLdata.y_def_interp = y_def_interp;
SOLdata.z_def_interp = z_def_interp;
