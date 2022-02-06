function varargout = unpack_FEMdata(FEMdata,varargin)

%% Get outputs
if isempty(varargin) % Unpack all data
    % General data
    unit_sys = FEMdata.unit_sys;
    beam_theory = FEMdata.beam_theory;
    N_modes = FEMdata.N_modes;
    RI = FEMdata.RI;
    N_beams = FEMdata.N_beams;
    L = FEMdata.L;
    b_alpha = FEMdata.b_alpha;
    b_beta = FEMdata.b_beta;
    b_gamma = FEMdata.b_gamma;
    constitutive_model = FEMdata.constitutive_model;
    Ne_b = FEMdata.Ne_b;
    element_order = FEMdata.element_order;
    elem_connect = FEMdata.elem_connect;
    BC_nodes = FEMdata.BC_nodes;
    Ne = FEMdata.Ne;
    N_nodes = FEMdata.N_nodes;
    Ndof = FEMdata.Ndof;
    edof = FEMdata.edof;
    enn = FEMdata.enn;
    ndof = FEMdata.ndof;
    elem_nodes = FEMdata.elem_nodes;
    nodes_coords = FEMdata.nodes_coords;
    n_div = FEMdata.n_div;
    EDOFs = FEMdata.EDOFs;
    DOF_u = FEMdata.DOF_u;
    DOF_v = FEMdata.DOF_v;
    DOF_w = FEMdata.DOF_w;
    DOF_phix = FEMdata.DOF_phix;
    DOF_phiy = FEMdata.DOF_phiy;
    DOF_phiz = FEMdata.DOF_phiz;
    DOF_dphix = FEMdata.DOF_dphix;
    NGP = FEMdata.NGP;
    dirac_delta = FEMdata.dirac_delta;  
    % Concentraded sources - numerical values
    CSx = FEMdata.sources.CSx;
    CSy = FEMdata.sources.CSy;
    CSz = FEMdata.sources.CSz;
    CSu = FEMdata.sources.CSu;
    CSv = FEMdata.sources.CSv;
    CSw = FEMdata.sources.CSw;
    CSt = FEMdata.sources.CSt;
    CSbl = FEMdata.sources.CSbl;
    CSbt = FEMdata.sources.CSbt;
    Ma = FEMdata.sources.Ma;
    MaRix = FEMdata.sources.MaRix;
    MaRiy = FEMdata.sources.MaRiy;
    MaRiz = FEMdata.sources.MaRiz;
    % Concentraded sources - as functions of x
    CSx_of_x = FEMdata.sources.CSx_of_x;
    CSy_of_x = FEMdata.sources.CSy_of_x;
    CSz_of_x = FEMdata.sources.CSz_of_x;
    CSu_of_x = FEMdata.sources.CSu_of_x;
    CSv_of_x = FEMdata.sources.CSv_of_x;
    CSw_of_x = FEMdata.sources.CSw_of_x;
    CSt_of_x = FEMdata.sources.CSt_of_x;
    CSbl_of_x = FEMdata.sources.CSbl_of_x;
    CSbt_of_x = FEMdata.sources.CSbt_of_x;
    Ma_of_x = FEMdata.sources.Ma_of_x;
    MaRix_of_x = FEMdata.sources.MaRix_of_x;
    MaRiy_of_x = FEMdata.sources.MaRiy_of_x;
    MaRiz_of_x = FEMdata.sources.MaRiz_of_x;
    % Concentraded sources positions
    akx = FEMdata.sources.akx;
    aky = FEMdata.sources.aky;
    akz = FEMdata.sources.akz;
    aku = FEMdata.sources.aku;
    akv = FEMdata.sources.akv;
    akw = FEMdata.sources.akw;
    akt = FEMdata.sources.akt;
    akbl = FEMdata.sources.akbl;
    akbt = FEMdata.sources.akbt;
    ama = FEMdata.sources.ama;
    % Distributed sources
    cfl_of_x = FEMdata.sources.cfl_of_x;
    cft_of_x = FEMdata.sources.cft_of_x;
    % Interpolation, constitutive, force and strain-displacement function matrices
    psi = FEMdata.funs.psi;
    dpsi = FEMdata.funs.dpsi;
    d2psi = FEMdata.funs.d2psi;
    d3psi = FEMdata.funs.d3psi;
    phi = FEMdata.funs.phi;
    dphi = FEMdata.funs.dphi;
    d2phi = FEMdata.funs.d2phi;
    zeta = FEMdata.funs.zeta;
    dzeta = FEMdata.funs.dzeta;
    B = FEMdata.funs.B;
    B_cf = FEMdata.funs.B_cf;
    B_Mt = FEMdata.funs.B_Mt;
    B_Mr = FEMdata.funs.B_Mr;
    H_psi = FEMdata.funs.H_psi;
    H_phi = FEMdata.funs.H_phi;
    H_zeta = FEMdata.funs.H_zeta;
    D = FEMdata.funs.D;
    D_cf = FEMdata.funs.D_cf;
    D_Mt = FEMdata.funs.D_Mt;
    D_Mr = FEMdata.funs.D_Mr;
    % Outputs:
    varargout = {unit_sys,beam_theory,N_modes,RI,N_beams,L,b_alpha,b_beta,b_gamma,constitutive_model,Ne_b,element_order,elem_connect,BC_nodes,Ne,N_nodes,Ndof,edof,enn,ndof,elem_nodes,nodes_coords,n_div,EDOFs,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,NGP,dirac_delta,CSx,CSy,CSz,CSu,CSv,CSw,CSt,CSbl,CSbt,Ma,MaRix,MaRiy,MaRiz,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,CSbl_of_x,CSbt_of_x,Ma_of_x,MaRix_of_x,MaRiy_of_x,MaRiz_of_x,akx,aky,akz,aku,akv,akw,akt,akbl,akbt,ama,cfl_of_x,cft_of_x,psi,dpsi,d2psi,d3psi,phi,dphi,d2phi,zeta,dzeta,B,B_cf,B_Mt,B_Mr,H_psi,H_phi,H_zeta,D,D_cf,D_Mt,D_Mr};
else % Unpack specific data for function
    switch lower(varargin{1})
          case 'global_matrices'
              % General data
              Ne = FEMdata.Ne;
              Ndof = FEMdata.Ndof;
              edof = FEMdata.edof;
              enn = FEMdata.enn;
              EDOFs = FEMdata.EDOFs;
              NGP = FEMdata.NGP;
              RI = FEMdata.RI;
              % Interpolation, constitutive, force and strain-displacement function matrices
              B = FEMdata.funs.B;
              B_cf = FEMdata.funs.B_cf;
              B_Mt = FEMdata.funs.B_Mt;
              B_Mr = FEMdata.funs.B_Mr;
              D = FEMdata.funs.D;
              D_cf = FEMdata.funs.D_cf;
              D_Mt = FEMdata.funs.D_Mt;
              D_Mr = FEMdata.funs.D_Mr;
              % Outputs:
              varargout = {Ne,Ndof,edof,enn,EDOFs,NGP,RI,B,B_cf,B_Mt,B_Mr,D,D_cf,D_Mt,D_Mr};        
        case 'concentraded_sources'
              % General data
              edof = FEMdata.edof;
              enn = FEMdata.enn;
              DOF_u = FEMdata.DOF_u;
              DOF_v = FEMdata.DOF_v;
              DOF_w = FEMdata.DOF_w;
              DOF_phix = FEMdata.DOF_phix;
              DOF_phiy = FEMdata.DOF_phiy;
              DOF_phiz = FEMdata.DOF_phiz;
              DOF_dphix = FEMdata.DOF_dphix;
              % Interpolation function matrices
              H_psi = FEMdata.funs.H_psi;
              H_phi = FEMdata.funs.H_phi;
              H_zeta = FEMdata.funs.H_zeta;
              % Outputs:
              varargout = {edof,enn,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,H_psi,H_phi,H_zeta};
        case 'bcs'
              % General data
              BC_nodes = FEMdata.BC_nodes;
              N_nodes = FEMdata.N_nodes;
              ndof = FEMdata.ndof;
              elem_nodes = FEMdata.elem_nodes;
              EDOFs = FEMdata.EDOFs;
              % Outputs:
              varargout = {BC_nodes,N_nodes,ndof,elem_nodes,EDOFs};
         case 'process_modes'
              % General data
              unit_sys = FEMdata.unit_sys;
              N_modes = FEMdata.N_modes;
              Ndof = FEMdata.Ndof;
              acc_modes = FEMdata.acc_modes;
              remove_list = FEMdata.remove_list;
              % Outputs:
              varargout = {unit_sys,N_modes,Ndof,acc_modes,remove_list};
         case 'outputs'
              % General data
              N_modes = FEMdata.N_modes;
              Ne = FEMdata.Ne;
              enn = FEMdata.enn;
              n_div = FEMdata.n_div;
              DOF_u = FEMdata.DOF_u;
              DOF_v = FEMdata.DOF_v;
              DOF_w = FEMdata.DOF_w;
              DOF_phix = FEMdata.DOF_phix;
              DOF_phiy = FEMdata.DOF_phiy;
              DOF_phiz = FEMdata.DOF_phiz;
              DOF_dphix = FEMdata.DOF_dphix;
              % Interpolation functions
              psi = FEMdata.funs.psi;
              dpsi = FEMdata.funs.dpsi;
              phi = FEMdata.funs.phi;
              zeta = FEMdata.funs.zeta;
              % Outputs:
              varargout = {N_modes,Ne,enn,n_div,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_dphix,psi,dpsi,phi,zeta};         
        case 'plots'
            % General data
            unit_sys = FEMdata.unit_sys;
            L = FEMdata.L;
            b_alpha = FEMdata.b_alpha;
            b_beta = FEMdata.b_beta;
            b_gamma = FEMdata.b_gamma;
            BC_nodes = FEMdata.BC_nodes;
            Ne = FEMdata.Ne;
            N_modes = FEMdata.N_modes;
            % Concentraded sources positions
            akx = FEMdata.sources.akx;
            aky = FEMdata.sources.aky;
            akz = FEMdata.sources.akz;
            aku = FEMdata.sources.aku;
            akv = FEMdata.sources.akv;
            akw = FEMdata.sources.akw;
            akt = FEMdata.sources.akt;
            akbl = FEMdata.sources.akbl;
            akbt = FEMdata.sources.akbt;
            ama = FEMdata.sources.ama;
            % Distributed sources
            cfl_of_x = FEMdata.sources.cfl_of_x;
            cft_of_x = FEMdata.sources.cft_of_x;
            % Outputs:
            varargout = {unit_sys,L,b_alpha,b_beta,b_gamma,BC_nodes,Ne,N_modes,akx,aky,akz,aku,akv,akw,akt,akbl,akbt,ama,cfl_of_x,cft_of_x};
        otherwise
              error(['Unexpected option: ' varargin{1}])
    end
end

