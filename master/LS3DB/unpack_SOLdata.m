function varargout = unpack_SOLdata(SOLdata,e,varargin)

%% Get outputs
if isempty(varargin) % Unpack all data
    u_global = SOLdata.u_global;
    u_local = SOLdata.u_local{e};
    u_local_nodes = SOLdata.u_local_nodes{e};
    v_local_nodes = SOLdata.v_local_nodes{e};
    w_local_nodes = SOLdata.w_local_nodes{e};
    phix_local_nodes = SOLdata.phix_local_nodes{e};
    phiy_local_nodes = SOLdata.phiy_local_nodes{e};
    phiz_local_nodes = SOLdata.phiz_local_nodes{e};
    Nx_nodes = SOLdata.Nx_nodes{e};
    Vy_nodes = SOLdata.Vy_nodes{e};
    Vz_nodes = SOLdata.Vz_nodes{e};
    Tq_nodes = SOLdata.Tq_nodes{e};
    My_nodes = SOLdata.My_nodes{e};
    Mz_nodes = SOLdata.Mz_nodes{e};
    Bm_nodes = SOLdata.Bm_nodes{e};
    u_local_interp = SOLdata.u_local_interp{e};
    v_local_interp = SOLdata.v_local_interp{e};
    w_local_interp = SOLdata.w_local_interp{e};
    phix_local_interp = SOLdata.phix_local_interp{e};
    phiy_local_interp = SOLdata.phiy_local_interp{e};
    phiz_local_interp = SOLdata.phiz_local_interp{e};
    Nx_interp = SOLdata.Nx_interp{e};
    Vy_interp = SOLdata.Vy_interp{e};
    Vz_interp = SOLdata.Vz_interp{e};
    Tq_interp = SOLdata.Tq_interp{e};
    My_interp = SOLdata.My_interp{e};
    Mz_interp = SOLdata.Mz_interp{e};
    Bm_interp = SOLdata.Bm_interp{e};
    x_def_nodes = SOLdata.x_def_nodes{e};
    y_def_nodes = SOLdata.y_def_nodes{e};
    z_def_nodes = SOLdata.z_def_nodes{e};
    x_def_interp = SOLdata.x_def_interp{e};
    y_def_interp = SOLdata.y_def_interp{e};
    z_def_interp = SOLdata.z_def_interp{e};
    
    varargout = {u_global,u_local,u_local_nodes,v_local_nodes,w_local_nodes,phix_local_nodes,phiy_local_nodes,phiz_local_nodes,Nx_nodes,Vy_nodes,Vz_nodes,Tq_nodes,My_nodes,Mz_nodes,Bm_nodes,u_local_interp,v_local_interp,w_local_interp,phix_local_interp,phiy_local_interp,phiz_local_interp,Nx_interp,Vy_interp,Vz_interp,Tq_interp,My_interp,Mz_interp,Bm_interp,x_def_nodes,y_def_nodes,z_def_nodes,x_def_interp,y_def_interp,z_def_interp};
else % Unpack specific data for function
    switch lower(varargin{1})
        case 'plot_structure'
            x_def_nodes = SOLdata.x_def_nodes{e};
            y_def_nodes = SOLdata.y_def_nodes{e};
            z_def_nodes = SOLdata.z_def_nodes{e};
            x_def_interp = SOLdata.x_def_interp{e};
            y_def_interp = SOLdata.y_def_interp{e};
            z_def_interp = SOLdata.z_def_interp{e};
            
            varargout = {x_def_nodes,y_def_nodes,z_def_nodes,x_def_interp,y_def_interp,z_def_interp};
        case 'plot_u'
            u_local_nodes = SOLdata.u_local_nodes{e};
            u_local_interp = SOLdata.u_local_interp{e};
            
            varargout = {u_local_nodes,u_local_interp};            
        case 'plot_v'
            v_local_nodes = SOLdata.v_local_nodes{e};
            v_local_interp = SOLdata.v_local_interp{e};
            
            varargout = {v_local_nodes,v_local_interp};    
        case 'plot_w'
            w_local_nodes = SOLdata.w_local_nodes{e};
            w_local_interp = SOLdata.w_local_interp{e};
            
            varargout = {w_local_nodes,w_local_interp};    
        case 'plot_phix'
            phix_local_nodes = SOLdata.phix_local_nodes{e};
            phix_local_interp = SOLdata.phix_local_interp{e};
            
            varargout = {phix_local_nodes,phix_local_interp};    
        case 'plot_phiy'
            phiy_local_nodes = SOLdata.phiy_local_nodes{e};
            phiy_local_interp = SOLdata.phiy_local_interp{e};
            
            varargout = {phiy_local_nodes,phiy_local_interp};    
        case 'plot_phiz'
            phiz_local_nodes = SOLdata.phiz_local_nodes{e};
            phiz_local_interp = SOLdata.phiz_local_interp{e};
            
            varargout = {phiz_local_nodes,phiz_local_interp};    
        case 'plot_dphix'
            dphix_local_nodes = SOLdata.dphix_local_nodes{e};
            dphix_local_interp = SOLdata.dphix_local_interp{e};
            
            varargout = {dphix_local_nodes,dphix_local_interp};    
        case 'plot_nx'
            Nx_nodes = SOLdata.Nx_nodes{e};
            Nx_interp = SOLdata.Nx_interp{e};
            
            varargout = {Nx_nodes,Nx_interp};    
        case 'plot_vy'
            Vy_nodes = SOLdata.Vy_nodes{e};
            Vy_interp = SOLdata.Vy_interp{e};
            
            varargout = {Vy_nodes,Vy_interp};    
        case 'plot_vz'
            Vz_nodes = SOLdata.Vz_nodes{e};
            Vz_interp = SOLdata.Vz_interp{e};
            
            varargout = {Vz_nodes,Vz_interp};    
        case 'plot_tq'
            Tq_nodes = SOLdata.Tq_nodes{e};
            Tq_interp = SOLdata.Tq_interp{e};
            
            varargout = {Tq_nodes,Tq_interp};    
        case 'plot_my'
            My_nodes = SOLdata.My_nodes{e};
            My_interp = SOLdata.My_interp{e};
            
            varargout = {My_nodes,My_interp};    
        case 'plot_mz'
            Mz_nodes = SOLdata.Mz_nodes{e};
            Mz_interp = SOLdata.Mz_interp{e};
            
            varargout = {Mz_nodes,Mz_interp};    
        case 'plot_bm'
            Bm_nodes = SOLdata.Bm_nodes{e};
            Bm_interp = SOLdata.Bm_interp{e};
            
            varargout = {Bm_nodes,Bm_interp};    
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    
end