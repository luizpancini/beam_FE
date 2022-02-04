function varargout = unpack_FEMdata(FEMdata,varargin)

%% Get outputs
if isempty(varargin) % Unpack all data
    % General data
    unit_sys = FEMdata.unit_sys;
    beam_theory = FEMdata.beam_theory;
    warp_DOF = FEMdata.warp_DOF;
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
    scale = FEMdata.scale;
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
    DOF_bendl = FEMdata.DOF_bendl;
    DOF_bendt = FEMdata.DOF_bendt;
    NGP = FEMdata.NGP;
    dirac_delta = FEMdata.dirac_delta;
    % Concentraded loads - numerical values
    Px = FEMdata.loads.Px;
    Py = FEMdata.loads.Py;
    Pz = FEMdata.loads.Pz;
    Mx = FEMdata.loads.Mx;
    My = FEMdata.loads.My;
    Mz = FEMdata.loads.Mz;
    Pa = FEMdata.loads.Pa;
    Pl = FEMdata.loads.Pl;
    Pt = FEMdata.loads.Pt;
    Tq = FEMdata.loads.Tq;
    Ml = FEMdata.loads.Ml;
    Mt = FEMdata.loads.Mt;
    Bm = FEMdata.loads.Bm;
    % Concentraded loads - as functions of x
    Px_of_x = FEMdata.loads.Px_of_x;
    Py_of_x = FEMdata.loads.Py_of_x;
    Pz_of_x = FEMdata.loads.Pz_of_x;
    Mx_of_x = FEMdata.loads.Mx_of_x;
    My_of_x = FEMdata.loads.My_of_x;
    Mz_of_x = FEMdata.loads.Mz_of_x;
    Pa_of_x = FEMdata.loads.Pa_of_x;
    Pl_of_x = FEMdata.loads.Pl_of_x;
    Pt_of_x = FEMdata.loads.Pt_of_x;
    Tq_of_x = FEMdata.loads.Tq_of_x;
    Ml_of_x = FEMdata.loads.Ml_of_x;
    Mt_of_x = FEMdata.loads.Mt_of_x;
    Bm_of_x = FEMdata.loads.Bm_of_x; 
    % Concentrated loads positions
    apx = FEMdata.loads.apx;
    apy = FEMdata.loads.apy;
    apz = FEMdata.loads.apz;
    amx = FEMdata.loads.amx;
    amy = FEMdata.loads.amy;
    amz = FEMdata.loads.amz;
    apa = FEMdata.loads.apa;
    apl = FEMdata.loads.apl;
    apt = FEMdata.loads.apt;
    atq = FEMdata.loads.atq;
    aml = FEMdata.loads.aml;
    amt = FEMdata.loads.amt;
    abm = FEMdata.loads.abm;
    % Distributed loads
    fa_of_x = FEMdata.loads.fa_of_x;
    tq_of_x = FEMdata.loads.tq_of_x;
    ql_of_x = FEMdata.loads.ql_of_x;
    qt_of_x = FEMdata.loads.qt_of_x;
    fx_of_x = FEMdata.loads.fx_of_x;
    mx_of_x = FEMdata.loads.mx_of_x;
    qy_of_x = FEMdata.loads.qy_of_x;
    qz_of_x = FEMdata.loads.qz_of_x;
    bm_of_x = FEMdata.loads.bm_of_x;   
    % Concentraded sources - numerical values
    CSx = FEMdata.sources.CSx;
    CSy = FEMdata.sources.CSy;
    CSz = FEMdata.sources.CSz;
    CSu = FEMdata.sources.CSu;
    CSv = FEMdata.sources.CSv;
    CSw = FEMdata.sources.CSw;
    CSt = FEMdata.sources.CSt;
    % Concentraded sources - as functions of x
    CSx_of_x = FEMdata.sources.CSx_of_x;
    CSy_of_x = FEMdata.sources.CSy_of_x;
    CSz_of_x = FEMdata.sources.CSz_of_x;
    CSu_of_x = FEMdata.sources.CSu_of_x;
    CSv_of_x = FEMdata.sources.CSv_of_x;
    CSw_of_x = FEMdata.sources.CSw_of_x;
    CSt_of_x = FEMdata.sources.CSt_of_x;
    % Concentraded sources positions
    akx = FEMdata.sources.akx;
    aky = FEMdata.sources.aky;
    akz = FEMdata.sources.akz;
    aku = FEMdata.sources.aku;
    akv = FEMdata.sources.akv;
    akw = FEMdata.sources.akw;
    akt = FEMdata.sources.akt;
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
    H_dist = FEMdata.funs.H_dist;
    H_psi = FEMdata.funs.H_psi;
    H_phi = FEMdata.funs.H_phi;
    H_zeta = FEMdata.funs.H_zeta;
    f_dist_l = FEMdata.funs.f_dist_l;
    f_dist_g = FEMdata.funs.f_dist_g;
    D = FEMdata.funs.D;
    D_cf = FEMdata.funs.D_cf;
    % Outputs:
    varargout = {unit_sys,beam_theory,warp_DOF,N_beams,L,b_alpha,b_beta,b_gamma,constitutive_model,Ne_b,element_order,elem_connect,BC_nodes,scale,Ne,N_nodes,Ndof,edof,enn,ndof,elem_nodes,nodes_coords,n_div,EDOFs,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,DOF_bendl,DOF_bendt,NGP,dirac_delta,Px,Py,Pz,Mx,My,Mz,Pa,Pl,Pt,Tq,Ml,Mt,Bm,Px_of_x,Py_of_x,Pz_of_x,Mx_of_x,My_of_x,Mz_of_x,Pa_of_x,Pl_of_x,Pt_of_x,Tq_of_x,Ml_of_x,Mt_of_x,Bm_of_x,apx,apy,apz,amx,amy,amz,apa,apl,apt,atq,aml,amt,abm,fa_of_x,tq_of_x,ql_of_x,qt_of_x,fx_of_x,mx_of_x,qy_of_x,qz_of_x,bm_of_x,CSx,CSy,CSz,CSu,CSv,CSw,CSt,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,akx,aky,akz,aku,akv,akw,akt,cfl_of_x,cft_of_x,psi,dpsi,d2psi,d3psi,phi,dphi,d2phi,zeta,dzeta,B,B_cf,H_dist,H_psi,H_phi,H_zeta,f_dist_l,f_dist_g,D,D_cf};
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
              % Interpolation, constitutive, force and strain-displacement function matrices
              B = FEMdata.funs.B;
              B_cf = FEMdata.funs.B_cf;
              H_dist = FEMdata.funs.H_dist;
              f_dist_l = FEMdata.funs.f_dist_l;
              f_dist_g = FEMdata.funs.f_dist_g;
              D = FEMdata.funs.D;
              D_cf = FEMdata.funs.D_cf;
              % Outputs:
              varargout = {Ne,Ndof,edof,enn,EDOFs,NGP,B,B_cf,H_dist,f_dist_l,f_dist_g,D,D_cf};        
        case 'concentraded_loads'
            % General data
            warp_DOF = FEMdata.warp_DOF;
            edof = FEMdata.edof;
            enn = FEMdata.enn;
            DOF_u = FEMdata.DOF_u;
            DOF_phix = FEMdata.DOF_phix;
            DOF_bendl = FEMdata.DOF_bendl;
            DOF_bendt = FEMdata.DOF_bendt;
            % Interpolation, constitutive, force and strain-displacement function matrices
            H_psi = FEMdata.funs.H_psi;
            H_phi = FEMdata.funs.H_phi;
            H_zeta = FEMdata.funs.H_zeta;
            % Outputs:
            varargout = {warp_DOF,edof,enn,DOF_u,DOF_phix,DOF_bendl,DOF_bendt,H_psi,H_phi,H_zeta};
        case 'bcs'
              % General data
              warp_DOF = FEMdata.warp_DOF;
              BC_nodes = FEMdata.BC_nodes;
              N_nodes = FEMdata.N_nodes;
              ndof = FEMdata.ndof;
              elem_nodes = FEMdata.elem_nodes;
              EDOFs = FEMdata.EDOFs;
              DOF_u = FEMdata.DOF_u;
              DOF_phix = FEMdata.DOF_phix;
              DOF_bendl = FEMdata.DOF_bendl;
              DOF_bendt = FEMdata.DOF_bendt;
              % Outputs:
              varargout = {warp_DOF,BC_nodes,N_nodes,ndof,elem_nodes,EDOFs,DOF_u,DOF_phix,DOF_bendl,DOF_bendt};
          case 'outputs'
              % General data
              beam_theory = FEMdata.beam_theory;
              warp_DOF = FEMdata.warp_DOF;
              constitutive_model = FEMdata.constitutive_model;
              scale = FEMdata.scale;
              Ne = FEMdata.Ne;
              enn = FEMdata.enn;
              n_div = FEMdata.n_div;
              DOF_u = FEMdata.DOF_u;
              DOF_v = FEMdata.DOF_v;
              DOF_w = FEMdata.DOF_w;
              DOF_phix = FEMdata.DOF_phix;
              DOF_phiy = FEMdata.DOF_phiy;
              DOF_phiz = FEMdata.DOF_phiz;
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
              % Outputs:
              varargout = {beam_theory,warp_DOF,constitutive_model,scale,Ne,enn,n_div,DOF_u,DOF_v,DOF_w,DOF_phix,DOF_phiy,DOF_phiz,psi,dpsi,d2psi,d3psi,phi,dphi,d2phi,zeta,dzeta};         
        case 'plots'
            % General data
            unit_sys = FEMdata.unit_sys;
            warp_DOF = FEMdata.warp_DOF;
            N_beams = FEMdata.N_beams;
            L = FEMdata.L;
            b_alpha = FEMdata.b_alpha;
            b_beta = FEMdata.b_beta;
            b_gamma = FEMdata.b_gamma;
            BC_nodes = FEMdata.BC_nodes;
            scale = FEMdata.scale;
            Ne = FEMdata.Ne;
            % Concentraded loads - numerical values
            Px = FEMdata.loads.Px;
            Py = FEMdata.loads.Py;
            Pz = FEMdata.loads.Pz;
            Mx = FEMdata.loads.Mx;
            My = FEMdata.loads.My;
            Mz = FEMdata.loads.Mz;
            Pa = FEMdata.loads.Pa;
            Pl = FEMdata.loads.Pl;
            Pt = FEMdata.loads.Pt;
            Tq = FEMdata.loads.Tq;
            Ml = FEMdata.loads.Ml;
            Mt = FEMdata.loads.Mt;
            Bm = FEMdata.loads.Bm;
            % Concentrated loads positions
            apx = FEMdata.loads.apx;
            apy = FEMdata.loads.apy;
            apz = FEMdata.loads.apz;
            amx = FEMdata.loads.amx;
            amy = FEMdata.loads.amy;
            amz = FEMdata.loads.amz;
            apa = FEMdata.loads.apa;
            apl = FEMdata.loads.apl;
            apt = FEMdata.loads.apt;
            atq = FEMdata.loads.atq;
            aml = FEMdata.loads.aml;
            amt = FEMdata.loads.amt;
            abm = FEMdata.loads.abm;
            % Distributed loads
            fa_of_x = FEMdata.loads.fa_of_x;
            tq_of_x = FEMdata.loads.tq_of_x;
            ql_of_x = FEMdata.loads.ql_of_x;
            qt_of_x = FEMdata.loads.qt_of_x;
            fx_of_x = FEMdata.loads.fx_of_x;
            mx_of_x = FEMdata.loads.mx_of_x;
            qy_of_x = FEMdata.loads.qy_of_x;
            qz_of_x = FEMdata.loads.qz_of_x;
            bm_of_x = FEMdata.loads.bm_of_x;
            % Concentraded sources positions
            akx = FEMdata.sources.akx;
            aky = FEMdata.sources.aky;
            akz = FEMdata.sources.akz;
            aku = FEMdata.sources.aku;
            akv = FEMdata.sources.akv;
            akw = FEMdata.sources.akw;
            akt = FEMdata.sources.akt;
            % Distributed sources
            cfl_of_x = FEMdata.sources.cfl_of_x;
            cft_of_x = FEMdata.sources.cft_of_x;
            % Outputs:
            varargout = {unit_sys,warp_DOF,N_beams,L,b_alpha,b_beta,b_gamma,BC_nodes,scale,Ne,Px,Py,Pz,Mx,My,Mz,Pa,Pl,Pt,Tq,Ml,Mt,Bm,apx,apy,apz,amx,amy,amz,apa,apl,apt,atq,aml,amt,abm,fa_of_x,tq_of_x,ql_of_x,qt_of_x,fx_of_x,mx_of_x,qy_of_x,qz_of_x,bm_of_x,akx,aky,akz,aku,akv,akw,akt,cfl_of_x,cft_of_x};
        otherwise
              error(['Unexpected option: ' varargin{1}])
    end
end

