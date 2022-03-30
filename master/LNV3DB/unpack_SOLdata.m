function varargout = unpack_SOLdata(SOLdata,e,varargin)

%% Get outputs
if isempty(varargin) % Unpack all data
    
    N_modes = SOLdata.N_modes;
    omega = SOLdata.omega;
    W = SOLdata.W;
    Wr = SOLdata.Wr; 
    lambda = SOLdata.lambda;
    u_nodes = SOLdata.u_nodes{e,:};
    v_nodes = SOLdata.v_nodes{e,:};
    w_nodes = SOLdata.w_nodes{e,:};
    phix_nodes = SOLdata.phix_nodes{e,:};
    phiy_nodes = SOLdata.phiy_nodes{e,:};
    phiz_nodes = SOLdata.phiz_nodes{e,:};
    dphix_nodes = SOLdata.dphix_nodes{e,:};
    u_interp = SOLdata.u_interp{e,:};
    v_interp = SOLdata.v_interp{e,:};
    w_interp = SOLdata.w_interp{e,:};
    phix_interp = SOLdata.phix_interp{e,:};
    phiy_interp = SOLdata.phiy_interp{e,:};
    phiz_interp = SOLdata.phiz_interp{e,:};
    dphix_interp = SOLdata.dphix_interp{e,:};
    x_def_nodes = SOLdata.x_def_nodes{e,:};
    y_def_nodes = SOLdata.y_def_nodes{e,:};
    z_def_nodes = SOLdata.z_def_nodes{e,:};
    x_def_interp = SOLdata.x_def_interp{e,:};
    y_def_interp = SOLdata.y_def_interp{e,:};
    z_def_interp = SOLdata.z_def_interp{e,:};
    
    varargout = {N_modes,omega,W,Wr,lambda,u_nodes,v_nodes,w_nodes,phix_nodes,phiy_nodes,phiz_nodes,dphix_nodes,u_interp,v_interp,w_interp,phix_interp,phiy_interp,phiz_interp,dphix_interp,x_def_nodes,y_def_nodes,z_def_nodes,x_def_interp,y_def_interp,z_def_interp};
else % Unpack specific data for function
    mode = varargin{2};
    switch lower(varargin{1})
        case 'plot_structure'
            omega = SOLdata.omega;
            x_def_nodes = SOLdata.x_def_nodes{e,mode};
            y_def_nodes = SOLdata.y_def_nodes{e,mode};
            z_def_nodes = SOLdata.z_def_nodes{e,mode};
            x_def_interp = SOLdata.x_def_interp{e,mode};
            y_def_interp = SOLdata.y_def_interp{e,mode};
            z_def_interp = SOLdata.z_def_interp{e,mode};
            
            varargout = {omega,x_def_nodes,y_def_nodes,z_def_nodes,x_def_interp,y_def_interp,z_def_interp};
        case 'plot_u'
            u_nodes = SOLdata.u_nodes{e,mode};
            u_interp = SOLdata.u_interp{e,mode};
            
            varargout = {u_nodes,u_interp};            
        case 'plot_v'
            v_nodes = SOLdata.v_nodes{e,mode};
            v_interp = SOLdata.v_interp{e,mode};
            
            varargout = {v_nodes,v_interp};    
        case 'plot_w'
            w_nodes = SOLdata.w_nodes{e,mode};
            w_interp = SOLdata.w_interp{e,mode};
            
            varargout = {w_nodes,w_interp};    
        case 'plot_phix'
            phix_nodes = SOLdata.phix_nodes{e,mode};
            phix_interp = SOLdata.phix_interp{e,mode};
            
            varargout = {phix_nodes,phix_interp};    
        case 'plot_phiy'
            phiy_nodes = SOLdata.phiy_nodes{e,mode};
            phiy_interp = SOLdata.phiy_interp{e,mode};
            
            varargout = {phiy_nodes,phiy_interp};    
        case 'plot_phiz'
            phiz_nodes = SOLdata.phiz_nodes{e,mode};
            phiz_interp = SOLdata.phiz_interp{e,mode};
            
            varargout = {phiz_nodes,phiz_interp};    
        case 'plot_dphix'
            dphix_nodes = SOLdata.dphix_nodes{e,mode};
            dphix_interp = SOLdata.dphix_interp{e,mode};
            
            varargout = {dphix_nodes,dphix_interp};     
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    
end