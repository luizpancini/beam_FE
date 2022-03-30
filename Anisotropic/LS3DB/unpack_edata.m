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
    G_of_x = edata.G_of_x{e};
    J_of_x = edata.J_of_x{e};
    Gamma_of_x = edata.Gamma_of_x{e};
    Iyy_of_x = edata.Iyy_of_x{e};
    Izz_of_x = edata.Izz_of_x{e};
    Ksy_of_x = edata.Ksy_of_x{e};
    Ksz_of_x = edata.Ksz_of_x{e};
    Gamx_of_x = edata.Gamx_of_x{e};
    Gamy_of_x = edata.Gamy_of_x{e};
    Gamz_of_x = edata.Gamz_of_x{e};
    % Constitutive properties (xi)
    E_of_xi = edata.E_of_xi{e};
    A_of_xi = edata.A_of_xi{e};
    G_of_xi = edata.G_of_xi{e};
    J_of_xi = edata.J_of_xi{e};
    Gamma_of_xi = edata.Gamma_of_xi{e};
    Iyy_of_xi = edata.Iyy_of_xi{e};
    Izz_of_xi = edata.Izz_of_xi{e};
    Ksy_of_xi = edata.Ksy_of_xi{e};
    Ksz_of_xi = edata.Ksz_of_xi{e};
    Gamx_of_xi = edata.Gamx_of_xi{e};
    Gamy_of_xi = edata.Gamy_of_xi{e};
    Gamz_of_xi = edata.Gamz_of_xi{e};
    % Distributed loads (xi)
    fx_of_xi = edata.fx_of_xi{e};
    qy_of_xi = edata.qy_of_xi{e};
    qz_of_xi = edata.qz_of_xi{e};
    mx_of_xi = edata.mx_of_xi{e};
    fa_of_xi = edata.fa_of_xi{e};
    ql_of_xi = edata.ql_of_xi{e};
    qt_of_xi = edata.qt_of_xi{e};
    tq_of_xi = edata.tq_of_xi{e};
    bm_of_xi = edata.bm_of_xi{e};
    % Distributed sources (xi)
    cfl_of_xi = edata.cfl_of_xi{e};
    cft_of_xi = edata.cft_of_xi{e};
    cf_on_element = edata.cf_on_element{e};
    % Concentraded loads (xi)
    Px_of_xi = edata.Px_of_xi{e};
    Py_of_xi = edata.Py_of_xi{e};
    Pz_of_xi = edata.Pz_of_xi{e};
    Mx_of_xi = edata.Mx_of_xi{e};
    My_of_xi = edata.My_of_xi{e};
    Mz_of_xi = edata.Mz_of_xi{e};
    Pa_of_xi = edata.Pa_of_xi{e};
    Pl_of_xi = edata.Pl_of_xi{e};
    Pt_of_xi = edata.Pt_of_xi{e};
    Tq_of_xi = edata.Tq_of_xi{e};
    Ml_of_xi = edata.Ml_of_xi{e};
    Mt_of_xi = edata.Mt_of_xi{e};
    Bm_of_xi = edata.Bm_of_xi{e};
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
    % Concentraded loads positions (xi)
    xi_apx = edata.xi_apx{e};
    xi_apy = edata.xi_apy{e};
    xi_apz = edata.xi_apz{e};
    xi_amx = edata.xi_amx{e};
    xi_amy = edata.xi_amy{e};
    xi_amz = edata.xi_amz{e};
    xi_apa = edata.xi_apa{e};
    xi_apl = edata.xi_apl{e};
    xi_apt = edata.xi_apt{e};
    xi_atq = edata.xi_atq{e};
    xi_aml = edata.xi_aml{e};
    xi_amt = edata.xi_amt{e};
    xi_abm = edata.xi_abm{e};
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
    % Outputs:
    varargout = {e_dof_range,e_node_range,L0,x1,beam,e_on_beam,R0,T0,Jac,x_vec_nodes,x_vec_interp,x0,y0,z0,x_und_nodes,y_und_nodes,z_und_nodes,x_und_cont,y_und_cont,z_und_cont,E_of_x,G_of_x,A_of_x,J_of_x,Gamma_of_x,Iyy_of_x,Izz_of_x,Ksy_of_x,Ksz_of_x,Gamx_of_x,Gamy_of_x,Gamz_of_x,E_of_xi,G_of_xi,A_of_xi,J_of_xi,Gamma_of_xi,Iyy_of_xi,Izz_of_xi,Ksy_of_xi,Ksz_of_xi,Gamx_of_xi,Gamy_of_xi,Gamz_of_xi,fx_of_xi,qy_of_xi,qz_of_xi,mx_of_xi,fa_of_xi,ql_of_xi,qt_of_xi,tq_of_xi,bm_of_xi,cfl_of_xi,cft_of_xi,cf_on_element,Px_of_xi,Py_of_xi,Pz_of_xi,Mx_of_xi,My_of_xi,Mz_of_xi,Pa_of_xi,Pl_of_xi,Pt_of_xi,Ml_of_xi,Mt_of_xi,Tq_of_xi,Bm_of_xi,CSx_of_xi,CSy_of_xi,CSz_of_xi,CSu_of_xi,CSv_of_xi,CSw_of_xi,CSt_of_xi,CSbl_of_xi,CSbt_of_xi,xi_apx,xi_apy,xi_apz,xi_amx,xi_amy,xi_amz,xi_apa,xi_apl,xi_apt,xi_atq,xi_aml,xi_amt,xi_abm,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt,xi_akbl,xi_akbt};
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
            props = edata.e_props{e};
            Gamx_of_xi = edata.e_props{e}.Gamx_of_xi;
            Gamy_of_xi = edata.e_props{e}.Gamy_of_xi;
            Gamz_of_xi = edata.e_props{e}.Gamz_of_xi;
            % Distributed loads and sources (xi)
            fx_of_xi = edata.fx_of_xi{e};
            qy_of_xi = edata.qy_of_xi{e};
            qz_of_xi = edata.qz_of_xi{e};
            mx_of_xi = edata.mx_of_xi{e};
            fa_of_xi = edata.fa_of_xi{e};
            ql_of_xi = edata.ql_of_xi{e};
            qt_of_xi = edata.qt_of_xi{e};
            tq_of_xi = edata.tq_of_xi{e};
            bm_of_xi = edata.bm_of_xi{e};
            cfl_of_xi = edata.cfl_of_xi{e};
            cft_of_xi = edata.cft_of_xi{e};
            cf_on_element = edata.cf_on_element{e};
            % Outputs:
            varargout = {e_dof_range,L0,beam,T0,Jac,props,Gamx_of_xi,Gamy_of_xi,Gamz_of_xi,fx_of_xi,qy_of_xi,qz_of_xi,mx_of_xi,fa_of_xi,ql_of_xi,qt_of_xi,tq_of_xi,bm_of_xi,cfl_of_xi,cft_of_xi,cf_on_element};     
        case 'concentraded_loads'
            % Geometry-related
            L0 = edata.L0{e};
            % Constitutive properties (xi)
            Gamx_of_xi = edata.e_props{e}.Gamx_of_xi;
            Gamy_of_xi = edata.e_props{e}.Gamy_of_xi;
            Gamz_of_xi = edata.e_props{e}.Gamz_of_xi;
            % Concentraded loads (xi)
            Px_of_xi = edata.Px_of_xi{e};
            Py_of_xi = edata.Py_of_xi{e};
            Pz_of_xi = edata.Pz_of_xi{e};
            Mx_of_xi = edata.Mx_of_xi{e};
            My_of_xi = edata.My_of_xi{e};
            Mz_of_xi = edata.Mz_of_xi{e};
            Pa_of_xi = edata.Pa_of_xi{e};
            Pl_of_xi = edata.Pl_of_xi{e};
            Pt_of_xi = edata.Pt_of_xi{e};
            Tq_of_xi = edata.Tq_of_xi{e};
            Ml_of_xi = edata.Ml_of_xi{e};
            Mt_of_xi = edata.Mt_of_xi{e};
            Bm_of_xi = edata.Bm_of_xi{e};
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
            % Concentraded loads positions (xi)
            xi_apx = edata.xi_apx{e};
            xi_apy = edata.xi_apy{e};
            xi_apz = edata.xi_apz{e};
            xi_amx = edata.xi_amx{e};
            xi_amy = edata.xi_amy{e};
            xi_amz = edata.xi_amz{e};
            xi_apa = edata.xi_apa{e};
            xi_apl = edata.xi_apl{e};
            xi_apt = edata.xi_apt{e};
            xi_atq = edata.xi_atq{e};
            xi_aml = edata.xi_aml{e};
            xi_amt = edata.xi_amt{e};
            xi_abm = edata.xi_abm{e};
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
            % Outputs:
            varargout = {L0,Gamx_of_xi,Gamy_of_xi,Gamz_of_xi,Px_of_xi,Py_of_xi,Pz_of_xi,Mx_of_xi,My_of_xi,Mz_of_xi,Pa_of_xi,Pl_of_xi,Pt_of_xi,Ml_of_xi,Mt_of_xi,Tq_of_xi,Bm_of_xi,CSx_of_xi,CSy_of_xi,CSz_of_xi,CSu_of_xi,CSv_of_xi,CSw_of_xi,CSt_of_xi,CSbl_of_xi,CSbt_of_xi,xi_apx,xi_apy,xi_apz,xi_amx,xi_amy,xi_amz,xi_apa,xi_apl,xi_apt,xi_atq,xi_aml,xi_amt,xi_abm,xi_akx,xi_aky,xi_akz,xi_aku,xi_akv,xi_akw,xi_akt,xi_akbl,xi_akbt};
        case 'outputs'
            % Geometry-related
            e_dof_range = edata.e_dof_range{e};
            L0 = edata.L0{e};
            x1 = edata.x1{e};
            beam = edata.beam{e};
            R0 = edata.R0{e};
            T0 = edata.T0{e};
            Jac = edata.Jac{e};
            x_vec_nodes = edata.x_vec_nodes{e};
            x_vec_interp = edata.x_vec_interp{e};
            x0 = edata.x0{e};
            y0 = edata.y0{e};
            z0 = edata.z0{e};
            % Constitutive properties (x)
            props = edata.e_props{e};
            Gamx_of_x = edata.e_props{e}.Gamx_of_x;
            Gamy_of_x = edata.e_props{e}.Gamy_of_x;
            Gamz_of_x = edata.e_props{e}.Gamz_of_x;
            % Outputs:
            varargout = {e_dof_range,L0,x1,beam,R0,T0,Jac,x_vec_nodes,x_vec_interp,x0,y0,z0,props,Gamx_of_x,Gamy_of_x,Gamz_of_x};
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
            beam = edata.beam{e};
            x_vec_nodes = edata.x_vec_nodes{e};
            x_vec_interp = edata.x_vec_interp{e};
            % Outputs:
            varargout = {beam,x_vec_nodes,x_vec_interp};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
end
