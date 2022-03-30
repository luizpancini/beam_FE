function varargout = unpack_edata(edata,e,varargin)

%% Get outputs
if isempty(varargin) % Unpack all data
    % Geometry-related
    e_dof_range = edata.e_dof_range{e};
    e_node_range = edata.e_node_range{e};
    L0 = edata.L0{e};
    x1 = edata.x1{e};
    beam = edata.beam{e};
    e_on_beam = edata.e_on_beam{e};
    R0 = edata.R0{e};
    T0 = edata.T0{e};
    Jac = edata.Jac{e};
    x_vec_nodes = edata.x_vec_nodes{e};
    x_vec_interp = edata.x_vec_interp{e};
    x0 = edata.x0{e};
    y0 = edata.y0{e};
    z0 = edata.z0{e};
    x_und_nodes = edata.x_und_nodes{e};
    y_und_nodes = edata.y_und_nodes{e};
    z_und_nodes = edata.z_und_nodes{e};
    x_und_cont = edata.x_und_cont{e};
    y_und_cont = edata.y_und_cont{e};
    z_und_cont = edata.z_und_cont{e};
    % Constitutive properties (x)
    E_of_x = edata.E_of_x{e};
    A_of_x = edata.A_of_x{e};
    rho_of_x = edata.rho_of_x{e};
    G_of_x = edata.G_of_x{e};
    J_of_x = edata.J_of_x{e};
    Gamma_of_x = edata.Gamma_of_x{e};
    Is_of_x = edata.Is_of_x{e};
    Iyy_of_x = edata.Iyy_of_x{e};
    Izz_of_x = edata.Izz_of_x{e};
    Ksy_of_x = edata.Ksy_of_x{e};
    Ksz_of_x = edata.Ksz_of_x{e};
    e1_of_x = edata.e1_of_x{e};
    e2_of_x = edata.e2_of_x{e};
    Gamx_of_x = edata.Gamx_of_x{e};
    Gamy_of_x = edata.Gamy_of_x{e};
    Gamz_of_x = edata.Gamz_of_x{e};
    % Constitutive properties (xi)
    E_of_xi = edata.E_of_xi{e};
    A_of_xi = edata.A_of_xi{e};
    rho_of_xi = edata.rho_of_xi{e};
    G_of_xi = edata.G_of_xi{e};
    J_of_xi = edata.J_of_xi{e};
    Gamma_of_xi = edata.Gamma_of_xi{e};
    Is_of_xi = edata.Is_of_xi{e};
    Iyy_of_xi = edata.Iyy_of_xi{e};
    Izz_of_xi = edata.Izz_of_xi{e};
    Ksy_of_xi = edata.Ksy_of_xi{e};
    Ksz_of_xi = edata.Ksz_of_xi{e};
    e1_of_xi = edata.e1_of_xi{e};
    e2_of_xi = edata.e2_of_xi{e};
    Gamx_of_xi = edata.Gamx_of_xi{e};
    Gamy_of_xi = edata.Gamy_of_xi{e};
    Gamz_of_xi = edata.Gamz_of_xi{e};
    % Distributed sources (xi)
    cfl_of_xi = edata.cfl_of_xi{e};
    cft_of_xi = edata.cft_of_xi{e};
    cf_on_element = edata.cf_on_element{e};
    % Concentraded sources (xi)
    CSx_of_xi = edata.CSx_of_xi{e};
    CSy_of_xi = edata.CSy_of_xi{e};
    CSz_of_xi = edata.CSz_of_xi{e};
    CSu_of_xi = edata.CSu_of_xi{e};
    CSv_of_xi = edata.CSv_of_xi{e};
    CSw_of_xi = edata.CSw_of_xi{e};
    CSt_of_xi = edata.CSt_of_xi{e};
    CSbl_of_xi = edata.CSbl_of_xi{e};
    CSbt_of_xi = edata.CSbt_of_xi{e};
    Ma_of_xi = edata.Ma_of_xi{e};
    MaRix_of_xi = edata.MaRix_of_xi{e};
    MaRiy_of_xi = edata.MaRiy_of_xi{e};
    MaRiz_of_xi = edata.MaRiz_of_xi{e};
    % Concentraded sources positions (xi)
    xi_akx = edata.xi_akx{e};
    xi_aky = edata.xi_aky{e};
    xi_akz = edata.xi_akz{e};
    xi_aku = edata.xi_aku{e};
    xi_akv = edata.xi_akv{e};
    xi_akw = edata.xi_akw{e};
    xi_akt = edata.xi_akt{e};
    xi_akbl = edata.xi_akbl{e};
    xi_akbt = edata.xi_akbt{e};
    xi_ama = edata.xi_ama{e};
    % Outputs:
    varargout = {e_dof_range,e_node_range,L0,x1,beam,e_on_beam,R0,T0,Jac,x_vec_nodes,x_vec_interp,x0,y0,z0,x_und_nodes,y_und_nodes,z_und_nodes,x_und_cont,y_und_cont,z_und_cont,E_of_x,G_of_x,rho_of_x,A_of_x,J_of_x,Gamma_of_x,Is_of_x,Iyy_of_x,Izz_of_x,Ksy_of_x,Ksz_of_x,e1_of_x,e2_of_x,Gamx_of_x,Gamy_of_x,Gamz_of_x,E_of_xi,G_of_xi,rho_of_xi,A_of_xi,J_of_xi,Gamma_of_xi,Is_of_xi,Iyy_of_xi,Izz_of_xi,Ksy_of_xi,Ksz_of_xi,e1_of_xi,e2_of_xi,Gamx_of_xi,Gamy_of_xi,Gamz_of_xi,cfl_of_xi,cft_of_xi,cf_on_element,CSx_of_xi,CSy_of_xi,CSz_of_xi,CSu_of_xi,CSv_of_xi,CSw_of_xi,CSt_of_xi,CSbl_of_xi,CSbt_of_xi,Ma_of_xi,MaRix_of_xi,MaRiy_of_xi,MaRiz_of_xi,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt,xi_akbl,xi_akbt,xi_ama};
else % Unpack specific data for function
    switch lower(varargin{1})
        case 'global_matrices'
            % Geometry-related
            e_dof_range = edata.e_dof_range{e};
            L0 = edata.L0{e};
            beam = edata.beam{e};
            T0 = edata.T0{e};
            Jac = edata.Jac{e};
            % Constitutive properties (xi)
            E_of_xi = edata.E_of_xi{e};
            A_of_xi = edata.A_of_xi{e};
            rho_of_xi = edata.rho_of_xi{e};
            G_of_xi = edata.G_of_xi{e};
            J_of_xi = edata.J_of_xi{e};
            Gamma_of_xi = edata.Gamma_of_xi{e};
            Is_of_xi = edata.Is_of_xi{e};
            Iyy_of_xi = edata.Iyy_of_xi{e};
            Izz_of_xi = edata.Izz_of_xi{e};
            Ksy_of_xi = edata.Ksy_of_xi{e};
            Ksz_of_xi = edata.Ksz_of_xi{e};
            e1_of_xi = edata.e1_of_xi{e};
            e2_of_xi = edata.e2_of_xi{e};
            Gamx_of_xi = edata.Gamx_of_xi{e};
            Gamy_of_xi = edata.Gamy_of_xi{e};
            Gamz_of_xi = edata.Gamz_of_xi{e};
            % Distributed sources (xi)
            cfl_of_xi = edata.cfl_of_xi{e};
            cft_of_xi = edata.cft_of_xi{e};
            cf_on_element = edata.cf_on_element{e};
            % Outputs:
            varargout = {e_dof_range,L0,beam,T0,Jac,E_of_xi,G_of_xi,rho_of_xi,A_of_xi,J_of_xi,Gamma_of_xi,Is_of_xi,Iyy_of_xi,Izz_of_xi,Ksy_of_xi,Ksz_of_xi,e1_of_xi,e2_of_xi,Gamx_of_xi,Gamy_of_xi,Gamz_of_xi,cfl_of_xi,cft_of_xi,cf_on_element};     
        case 'concentraded_sources'
            % Geometry-related
            L0 = edata.L0{e};
            % Constitutive properties (xi)
            Gamx_of_xi = edata.Gamx_of_xi{e};
            Gamy_of_xi = edata.Gamy_of_xi{e};
            Gamz_of_xi = edata.Gamz_of_xi{e}; 
            % Concentraded sources (xi)
            CSx_of_xi = edata.CSx_of_xi{e};
            CSy_of_xi = edata.CSy_of_xi{e};
            CSz_of_xi = edata.CSz_of_xi{e};
            CSu_of_xi = edata.CSu_of_xi{e};
            CSv_of_xi = edata.CSv_of_xi{e};
            CSw_of_xi = edata.CSw_of_xi{e};
            CSt_of_xi = edata.CSt_of_xi{e};
            CSbl_of_xi = edata.CSbl_of_xi{e};
            CSbt_of_xi = edata.CSbt_of_xi{e};
            Ma_of_xi = edata.Ma_of_xi{e};
            MaRix_of_xi = edata.MaRix_of_xi{e};
            MaRiy_of_xi = edata.MaRiy_of_xi{e};
            MaRiz_of_xi = edata.MaRiz_of_xi{e};
            % Concentraded sources positions (xi)
            xi_akx = edata.xi_akx{e};
            xi_aky = edata.xi_aky{e};
            xi_akz = edata.xi_akz{e};
            xi_aku = edata.xi_aku{e};
            xi_akv = edata.xi_akv{e};
            xi_akw = edata.xi_akw{e};
            xi_akt = edata.xi_akt{e};
            xi_akbl = edata.xi_akbl{e};
            xi_akbt = edata.xi_akbt{e};
            xi_ama = edata.xi_ama{e};
            % Outputs:
            varargout = {L0,Gamx_of_xi,Gamy_of_xi,Gamz_of_xi,CSx_of_xi,CSy_of_xi,CSz_of_xi,CSu_of_xi,CSv_of_xi,CSw_of_xi,CSt_of_xi,CSbl_of_xi,CSbt_of_xi,Ma_of_xi,MaRix_of_xi,MaRiy_of_xi,MaRiz_of_xi,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt,xi_akbl,xi_akbt,xi_ama};
        case 'outputs'
            % Geometry-related
            e_dof_range = edata.e_dof_range{e};
            L0 = edata.L0{e};
            x1 = edata.x1{e};
            R0 = edata.R0{e};
            T0 = edata.T0{e};
            Jac = edata.Jac{e};
            x_vec_nodes = edata.x_vec_nodes{e};
            x_vec_interp = edata.x_vec_interp{e};
            x0 = edata.x0{e};
            y0 = edata.y0{e};
            z0 = edata.z0{e};
            % Constitutive properties (x)
            Gamx_of_x = edata.Gamx_of_x{e};
            Gamy_of_x = edata.Gamy_of_x{e};
            Gamz_of_x = edata.Gamz_of_x{e};
            % Outputs:
            varargout = {e_dof_range,L0,x1,R0,T0,Jac,x_vec_nodes,x_vec_interp,x0,y0,z0,Gamx_of_x,Gamy_of_x,Gamz_of_x};
        case 'plot_structure'
            % Geometry-related
            e_node_range = edata.e_node_range{e};
            beam = edata.beam{e};
            R0 = edata.R0{e};
            x_vec_interp = edata.x_vec_interp{e};
            x0 = edata.x0{e};
            y0 = edata.y0{e};
            z0 = edata.z0{e};
            x_und_nodes = edata.x_und_nodes{e};
            y_und_nodes = edata.y_und_nodes{e};
            z_und_nodes = edata.z_und_nodes{e};
            x_und_cont = edata.x_und_cont{e};
            y_und_cont = edata.y_und_cont{e};
            z_und_cont = edata.z_und_cont{e};
            % Outputs:
            varargout = {e_node_range,beam,R0,x_vec_interp,x0,y0,z0,x_und_nodes,y_und_nodes,z_und_nodes,x_und_cont,y_und_cont,z_und_cont};
        case 'plot_gen_outs'
            % Geometry-related
            x_vec_nodes = edata.x_vec_nodes{e};
            x_vec_interp = edata.x_vec_interp{e};
            % Outputs:
            varargout = {x_vec_nodes,x_vec_interp};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
end
