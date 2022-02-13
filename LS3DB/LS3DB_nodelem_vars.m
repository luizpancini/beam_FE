function [FEMdata,edata] = LS3DB_nodelem_vars(FEMdata,edata,E,G,A,J,Gamma,Iyy,Izz,Ksy,Ksz,fx_of_x,mx_of_x,qy_of_x,qz_of_x,fa_of_x,ql_of_x,qt_of_x,tq_of_x,bm_of_x,Px_of_x,Py_of_x,Pz_of_x,Pa_of_x,Pl_of_x,Pt_of_x,Mx_of_x,My_of_x,Mz_of_x,Ml_of_x,Mt_of_x,Tq_of_x,Bm_of_x,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,CSbl_of_x,CSbt_of_x,cfl_of_x,cft_of_x)

% Unpack FEM data
[warp_DOF,L,b_alpha,b_beta,b_gamma,Ne_b,element_order,elem_connect,elem_nodes,apx,apy,apz,amx,amy,amz,apa,apl,apt,atq,aml,amt,abm,akx,aky,akz,aku,akv,akw,akt,akbl,akbt] = unpack_FEMdata(FEMdata,'nodelem_vars');

%% Degrees of freedom data
% Element number of nodes
if element_order == "linear"
    enn = 2;         
elseif element_order == "quadratic"
    enn = 3;         
else
    error('Specify element_order as linear or quadratic');
end
% Node's number of degrees of freedom
if warp_DOF  % Allowance for the warping DOF
    ndof = 7;        
else
    ndof = 6;
end
% Element number of degrees of freedom for element matrices
edof = enn*ndof;
% Number of Gauss points for integration  
NGP = 2*enn;  

% Element's DOFs
EDOFs = cell(enn,1); % Each cell entry contains corresponding local DOFs of element's nodes
DOF_u = zeros(1,enn); DOF_v = DOF_u; DOF_w = DOF_u; DOF_phix = DOF_u; DOF_phiy = DOF_u; DOF_phiz = DOF_u; DOF_dphix = DOF_u; 
for node=1:enn
    EDOFs{node} = ndof*node-(ndof-1):ndof*node;
    DOF_u(node) = EDOFs{node}(1);
    DOF_v(node) = EDOFs{node}(2);
    DOF_w(node) = EDOFs{node}(3);
    DOF_phix(node) = EDOFs{node}(4); 
    DOF_phiy(node) = EDOFs{node}(5);
    DOF_phiz(node) = EDOFs{node}(6);
    if warp_DOF, DOF_dphix(node) = EDOFs{node}(7); else, DOF_dphix(node) = NaN; end
end

%% Element and nodal variables
% Elements' nodes
if elem_connect == "sequenced"              
    elem_nodes = zeros(sum(Ne_b),enn); % #rows = #elements, #columns = #nodes/element
    for i=1:sum(Ne_b)
        if enn == 2
            elem_nodes(i,:) = [i, i+1];
        elseif enn == 3
            j = 1+2*(i-1);
            elem_nodes(i,:) = [j, j+1, j+2];
        end
    end
else % Unsequenced, so check for errors
    if size(elem_nodes,2) ~= enn
        error("Specify exactly " + num2str(enn) + " nodes per " + element_order + " element in the matrix elem_nodes");   
    end
end
% Total number of elements, nodes and DOFs
N_nodes = max(max(elem_nodes));     % Total number of nodes
Ne = sum(Ne_b);                     % Total number of elements
Ndof = ndof*N_nodes;                % Total number of degrees-of-freedom   
% Initialize element properties and data
n_div = 11;                         % Default number of divisions for interpolation in each element
nodes_coords = zeros(N_nodes,3);    % Initialize nodal coordinates
L0 = cell(Ne,1); x1 = cell(Ne,1); beam = cell(Ne,1); e_on_beam = cell(Ne,1); e_node_range = cell(Ne,1); e_dof_range = cell(Ne,1); Jac = cell(Ne,1); R0 = cell(Ne,1); T0 = cell(Ne,1); x_vec_nodes = cell(Ne,1); x_vec_interp = cell(Ne,1); x_und_nodes = cell(Ne,1); y_und_nodes = cell(Ne,1); z_und_nodes = cell(Ne,1); x_und_cont = cell(Ne,1); y_und_cont = cell(Ne,1); z_und_cont = cell(Ne,1); x0 = cell(Ne,1); y0 = cell(Ne,1); z0 = cell(Ne,1);
E_of_x = cell(Ne,1); G_of_x = cell(Ne,1); A_of_x = cell(Ne,1); J_of_x = cell(Ne,1); Iyy_of_x = cell(Ne,1); Izz_of_x = cell(Ne,1); Gamma_of_x = cell(Ne,1); Ksy_of_x = cell(Ne,1); Ksz_of_x = cell(Ne,1); Gamx_of_x = cell(Ne,1); Gamy_of_x = cell(Ne,1); Gamz_of_x = cell(Ne,1);
E_of_xi = cell(Ne,1); G_of_xi = cell(Ne,1); A_of_xi = cell(Ne,1); J_of_xi = cell(Ne,1); Iyy_of_xi = cell(Ne,1); Izz_of_xi = cell(Ne,1); Gamma_of_xi = cell(Ne,1); Ksy_of_xi = cell(Ne,1); Ksz_of_xi = cell(Ne,1); Gamx_of_xi = cell(Ne,1); Gamy_of_xi = cell(Ne,1); Gamz_of_xi = cell(Ne,1);
fx_of_xi = cell(Ne,1); qy_of_xi = cell(Ne,1); qz_of_xi = cell(Ne,1); mx_of_xi = cell(Ne,1); Px_of_xi = cell(Ne,1); Py_of_xi = cell(Ne,1); Pz_of_xi = cell(Ne,1); My_of_xi = cell(Ne,1); Mz_of_xi = cell(Ne,1); Mx_of_xi = cell(Ne,1); CSx_of_xi = cell(Ne,1); CSy_of_xi = cell(Ne,1); CSz_of_xi = cell(Ne,1);
fa_of_xi = cell(Ne,1); ql_of_xi = cell(Ne,1); qt_of_xi = cell(Ne,1); tq_of_xi = cell(Ne,1); bm_of_xi = cell(Ne,1); Pa_of_xi = cell(Ne,1); Pl_of_xi = cell(Ne,1); Pt_of_xi = cell(Ne,1); Ml_of_xi = cell(Ne,1); Mt_of_xi = cell(Ne,1); Tq_of_xi = cell(Ne,1); Bm_of_xi = cell(Ne,1); CSu_of_xi = cell(Ne,1); CSv_of_xi = cell(Ne,1); CSw_of_xi = cell(Ne,1); CSt_of_xi = cell(Ne,1); CSbl_of_xi = cell(Ne,1); CSbt_of_xi = cell(Ne,1); cfl_of_xi = cell(Ne,1); cft_of_xi = cell(Ne,1); cf_on_element = cell(Ne,1);
xi_apx = cell(Ne,1); xi_apy = cell(Ne,1); xi_apz = cell(Ne,1); xi_amx = cell(Ne,1); xi_amy = cell(Ne,1); xi_amz = cell(Ne,1); xi_apa = cell(Ne,1); xi_apl = cell(Ne,1); xi_apt = cell(Ne,1); xi_atq = cell(Ne,1); xi_aml = cell(Ne,1); xi_amt = cell(Ne,1); xi_abm = cell(Ne,1); xi_akx = cell(Ne,1); xi_aky = cell(Ne,1); xi_akz = cell(Ne,1); xi_aku = cell(Ne,1); xi_akv = cell(Ne,1); xi_akw = cell(Ne,1); xi_akt = cell(Ne,1); xi_akbl = cell(Ne,1);  xi_akbt = cell(Ne,1);xi_ama = cell(Ne,1); 
indexat = @(expr, index) expr(index); % To get specific element index
% Loop over elements
for e=1:Ne
    % beam: beam to which element e belongs
    beam{e} = 1; 
    while sum(Ne_b(1:beam{e})) < e
        beam{e} = beam{e}+1;
    end
    % e_on_beam: which element of that beam e is
    if beam{e} == 1
        e_on_beam{e} = e;
    else
        e_on_beam{e} = e-sum(Ne_b(1:beam{e}-1));
    end
    % Element node range, undeformed length, starting position and Jacobian
    e_node_range{e} = elem_nodes(e,:);       
    L0{e} = L(beam{e})/Ne_b(beam{e}); 
    x1{e} = L(beam{e})/Ne_b(beam{e})*(e_on_beam{e}-1); 
    Jac{e} = L0{e}/2;
    % Element global DOFs range
    e_dof_range{e} = zeros(edof,1);
    for n=1:enn
        e_dof_range{e}(ndof*n-(ndof-1):ndof*n) = ndof*e_node_range{e}(n)-(ndof-1):ndof*e_node_range{e}(n);
    end 
    % Element rotation and compound rotation tensors 
    sa = sin(b_alpha(beam{e})); ca = cos(b_alpha(beam{e})); cb = cos(b_beta(beam{e})); sb = sin(b_beta(beam{e})); cg = cos(b_gamma(beam{e})); sg = sin(b_gamma(beam{e}));
    R0{e} = [cb*ca, ca*sg*sb - cg*sa, sg*sa + cg*ca*sb
             cb*sa, cg*ca + sg*sb*sa, cg*sb*sa - ca*sg
               -sb,            cb*sg,            cg*cb]; % Bryant angles rotation matrix
    T0{e} = zeros(edof);                                 % Compound rotation matrix                
    for n=1:enn
        if warp_DOF
            T0{e}(EDOFs{n},EDOFs{n}) = [R0{e} zeros(3,4); zeros(3) R0{e} zeros(3,1); zeros(1,6) 1];
        else
            T0{e}(EDOFs{n},EDOFs{n}) = [R0{e} zeros(3); zeros(3) R0{e}];
        end
    end  
    % Element undeformed nodal coordinates
    x_vec_nodes{e} = zeros(enn,1); x_vec_nodes{e}(1) = x1{e}; x_vec_nodes{e}(end) = x1{e}+L0{e};      % Vector of nodal coordinates in local frame
    x_vec_interp{e} = linspace(x1{e},x1{e}+L0{e},n_div)';                                             % Vector of interpolated coordinates in local frame
    nodes_coords(e_node_range{e}(end),:) = nodes_coords(e_node_range{e}(1),:) + L0{e}*R0{e}(:,1)';    % Matrix of all nodes' global coordinates
    for n=2:enn-1
        nodes_coords(elem_nodes(e,n),:) = (nodes_coords(e_node_range{e}(end),:) + nodes_coords(e_node_range{e}(1),:))/(enn-1);  % Equally spaced inside nodes
        x_vec_nodes{e}(n) = x1{e}+(n-1)*L0{e}/(enn-1);                                                                          % Equally spaced inside nodes
    end
    if e_on_beam{e} == 1
        x0{e} = nodes_coords(e_node_range{e}(1),1);            % Starting global x coordinate of element's beam
        y0{e} = nodes_coords(e_node_range{e}(1),2);            % Starting global y coordinate of element's beam
        z0{e} = nodes_coords(e_node_range{e}(1),3);            % Starting global z coordinate of element's beam
    else
        x0{e} = x0{e-1};            
        y0{e} = y0{e-1};            
        z0{e} = z0{e-1};            
    end
    x_und_nodes{e} = x0{e} + x_vec_nodes{e}*R0{e}(1,1);    % Global x coordinates of element nodes
    y_und_nodes{e} = y0{e} + x_vec_nodes{e}*R0{e}(2,1);    % Global y coordinates of element nodes
    z_und_nodes{e} = z0{e} + x_vec_nodes{e}*R0{e}(3,1);    % Global z coordinates of element nodes
    x_und_cont{e} = x0{e} + x_vec_interp{e}*R0{e}(1,1);
    y_und_cont{e} = y0{e} + x_vec_interp{e}*R0{e}(2,1);
    z_und_cont{e} = z0{e} + x_vec_interp{e}*R0{e}(3,1);
    % Element constitutive properties as functions of the local coordinate x
    E_of_x{e} = @(x) indexat(E(x),beam{e});
    G_of_x{e} = @(x) indexat(G(x),beam{e});
    A_of_x{e} = @(x) indexat(A(x),beam{e});
    J_of_x{e} = @(x) indexat(J(x),beam{e});
    Gamma_of_x{e} = @(x) indexat(Gamma(x),beam{e});
    Iyy_of_x{e} = @(x) indexat(Iyy(x),beam{e});
    Izz_of_x{e} = @(x) indexat(Izz(x),beam{e});
    Ksy_of_x{e} = @(x) indexat(Ksy(x),beam{e});
    Ksz_of_x{e} = @(x) indexat(Ksz(x),beam{e});
    Gamx_of_x{e} = @(x) 0;
    Gamy_of_x{e} = @(x) E_of_x{e}(x).*Izz_of_x{e}(x)./(Ksy_of_x{e}(x).*G_of_x{e}(x).*A_of_x{e}(x)*L0{e}^2);
    Gamz_of_x{e} = @(x) E_of_x{e}(x).*Iyy_of_x{e}(x)./(Ksz_of_x{e}(x).*G_of_x{e}(x).*A_of_x{e}(x)*L0{e}^2);
    % Element constitutive properties as functions of the parent coordinate xi
    x = @(xi) L0{e}/2*(1+xi)+x1{e}; 
    E_of_xi{e} = @(xi) indexat(E(x(xi)),beam{e});
    G_of_xi{e} = @(xi) indexat(G(x(xi)),beam{e});
    A_of_xi{e} = @(xi) indexat(A(x(xi)),beam{e});
    J_of_xi{e} = @(xi) indexat(J(x(xi)),beam{e});
    Gamma_of_xi{e} = @(xi) indexat(Gamma(x(xi)),beam{e});
    Iyy_of_xi{e} = @(xi) indexat(Iyy(x(xi)),beam{e});
    Izz_of_xi{e} = @(xi) indexat(Izz(x(xi)),beam{e});
    Ksy_of_xi{e} = @(xi) indexat(Ksy(x(xi)),beam{e});
    Ksz_of_xi{e} = @(xi) indexat(Ksz(x(xi)),beam{e});
    Gamx_of_xi{e} = @(xi) 0; 
    Gamy_of_xi{e} = @(xi) E_of_xi{e}(xi).*Izz_of_xi{e}(xi)./(Ksy_of_xi{e}(xi).*G_of_xi{e}(xi).*A_of_xi{e}(xi)*L0{e}^2);
    Gamz_of_xi{e} = @(xi) E_of_xi{e}(xi).*Iyy_of_xi{e}(xi)./(Ksz_of_xi{e}(xi).*G_of_xi{e}(xi).*A_of_xi{e}(xi)*L0{e}^2);
    % Element distributed loads as functions of the parent coordinate xi
    fx_of_xi{e} = @(xi) fx_of_x{beam{e}}(x(xi));
    qy_of_xi{e} = @(xi) qy_of_x{beam{e}}(x(xi));
    qz_of_xi{e} = @(xi) qz_of_x{beam{e}}(x(xi));
    mx_of_xi{e} = @(xi) mx_of_x{beam{e}}(x(xi));
    fa_of_xi{e} = @(xi) fa_of_x{beam{e}}(x(xi));
    ql_of_xi{e} = @(xi) ql_of_x{beam{e}}(x(xi));
    qt_of_xi{e} = @(xi) qt_of_x{beam{e}}(x(xi));
    tq_of_xi{e} = @(xi) tq_of_x{beam{e}}(x(xi));
    bm_of_xi{e} = @(xi) bm_of_x{beam{e}}(x(xi));
    % Element distributed sources as functions of the parent coordinate xi
    cfl_of_xi{e} = @(xi) cfl_of_x{beam{e}}(x(xi));
    cft_of_xi{e} = @(xi) cft_of_x{beam{e}}(x(xi));
    if func2str(cfl_of_x{beam{e}}) ~= "@(x)0" || func2str(cft_of_x{beam{e}}) ~= "@(x)0", cf_on_element{e} = 1; else, cf_on_element{e} = 0; end
    % Element concentraded loads as functions of the parent coordinate xi
    Px_of_xi{e} = @(xi) Px_of_x{beam{e}}(x(xi));
    Py_of_xi{e} = @(xi) Py_of_x{beam{e}}(x(xi));
    Pz_of_xi{e} = @(xi) Pz_of_x{beam{e}}(x(xi));
    Mx_of_xi{e} = @(xi) Mx_of_x{beam{e}}(x(xi));
    My_of_xi{e} = @(xi) My_of_x{beam{e}}(x(xi));
    Mz_of_xi{e} = @(xi) Mz_of_x{beam{e}}(x(xi));
    Pa_of_xi{e} = @(xi) Pa_of_x{beam{e}}(x(xi));
    Pl_of_xi{e} = @(xi) Pl_of_x{beam{e}}(x(xi)); 
    Pt_of_xi{e} = @(xi) Pt_of_x{beam{e}}(x(xi));
    Ml_of_xi{e} = @(xi) Ml_of_x{beam{e}}(x(xi));
    Mt_of_xi{e} = @(xi) Mt_of_x{beam{e}}(x(xi));
    Tq_of_xi{e} = @(xi) Tq_of_x{beam{e}}(x(xi));
    Bm_of_xi{e} = @(xi) Bm_of_x{beam{e}}(x(xi));
    % Element concentraded sources as functions of the parent coordinate xi
    CSx_of_xi{e} = @(xi) CSx_of_x{beam{e}}(x(xi));
    CSy_of_xi{e} = @(xi) CSy_of_x{beam{e}}(x(xi));
    CSz_of_xi{e} = @(xi) CSz_of_x{beam{e}}(x(xi));
    CSu_of_xi{e} = @(xi) CSu_of_x{beam{e}}(x(xi));
    CSv_of_xi{e} = @(xi) CSv_of_x{beam{e}}(x(xi));
    CSw_of_xi{e} = @(xi) CSw_of_x{beam{e}}(x(xi));
    CSt_of_xi{e} = @(xi) CSt_of_x{beam{e}}(x(xi));
    CSbl_of_xi{e} = @(xi) CSbl_of_x{beam{e}}(x(xi));
    CSbt_of_xi{e} = @(xi) CSbt_of_x{beam{e}}(x(xi));
    % Element concentraded load positions as functions of the parent coordinate xi and the respective valid indices
    [xi_apx{e},apx] = get_loads_sources_positions(apx,Jac{e},beam{e},x1{e});
    [xi_apy{e},apy] = get_loads_sources_positions(apy,Jac{e},beam{e},x1{e});
    [xi_apz{e},apz] = get_loads_sources_positions(apz,Jac{e},beam{e},x1{e});
    [xi_amx{e},amx] = get_loads_sources_positions(amx,Jac{e},beam{e},x1{e});
    [xi_amy{e},amy] = get_loads_sources_positions(amy,Jac{e},beam{e},x1{e});
    [xi_amz{e},amz] = get_loads_sources_positions(amz,Jac{e},beam{e},x1{e});
    [xi_apa{e},apa] = get_loads_sources_positions(apa,Jac{e},beam{e},x1{e});
    [xi_apl{e},apl] = get_loads_sources_positions(apl,Jac{e},beam{e},x1{e});
    [xi_apt{e},apt] = get_loads_sources_positions(apt,Jac{e},beam{e},x1{e});
    [xi_atq{e},atq] = get_loads_sources_positions(atq,Jac{e},beam{e},x1{e});
    [xi_aml{e},aml] = get_loads_sources_positions(aml,Jac{e},beam{e},x1{e});
    [xi_amt{e},amt] = get_loads_sources_positions(amt,Jac{e},beam{e},x1{e});
    [xi_abm{e},abm] = get_loads_sources_positions(abm,Jac{e},beam{e},x1{e});
    % Element concentraded sources positions as functions of the parent coordinate xi and the respective valid indices
    [xi_akx{e},akx] = get_loads_sources_positions(akx,Jac{e},beam{e},x1{e});
    [xi_aky{e},aky] = get_loads_sources_positions(aky,Jac{e},beam{e},x1{e});
    [xi_akz{e},akz] = get_loads_sources_positions(akz,Jac{e},beam{e},x1{e});
    [xi_aku{e},aku] = get_loads_sources_positions(aku,Jac{e},beam{e},x1{e});
    [xi_akv{e},akv] = get_loads_sources_positions(akv,Jac{e},beam{e},x1{e});
    [xi_akw{e},akw] = get_loads_sources_positions(akw,Jac{e},beam{e},x1{e});
    [xi_akt{e},akt] = get_loads_sources_positions(akt,Jac{e},beam{e},x1{e});
    [xi_akbl{e},akbl] = get_loads_sources_positions(akbl,Jac{e},beam{e},x1{e});
    [xi_akbt{e},akbt] = get_loads_sources_positions(akbt,Jac{e},beam{e},x1{e});
end
% Truncate small components of nodal coordinates
tol = 1e-6; nodes_coords(abs(nodes_coords)<tol) = 0;

%% Add variables to respective structures
% FEM data
FEMdata.enn = enn;
FEMdata.ndof = ndof;
FEMdata.edof = edof;
FEMdata.EDOFs = EDOFs;
FEMdata.DOF_u = DOF_u;
FEMdata.DOF_v = DOF_v;
FEMdata.DOF_w = DOF_w;
FEMdata.DOF_phix = DOF_phix;
FEMdata.DOF_phiy = DOF_phiy;
FEMdata.DOF_phiz = DOF_phiz;
FEMdata.DOF_dphix = DOF_dphix;
FEMdata.NGP = NGP;
FEMdata.N_nodes = N_nodes;
FEMdata.Ne = Ne;
FEMdata.Ndof = Ndof;
FEMdata.elem_nodes = elem_nodes;
FEMdata.nodes_coords = nodes_coords;
FEMdata.n_div = n_div;

% Elements' data
edata.e_dof_range = e_dof_range;
edata.e_node_range = e_node_range;
edata.L0 = L0;
edata.x1 = x1;
edata.beam = beam;
edata.e_on_beam = e_on_beam;
edata.R0 = R0;
edata.T0 = T0;
edata.Jac = Jac;
edata.x_vec_nodes = x_vec_nodes;
edata.x_vec_interp = x_vec_interp;
edata.x0 = x0;
edata.y0 = y0;
edata.z0 = z0;
edata.x_und_nodes = x_und_nodes;
edata.y_und_nodes = y_und_nodes;
edata.z_und_nodes = z_und_nodes;
edata.x_und_cont = x_und_cont;
edata.y_und_cont = y_und_cont;
edata.z_und_cont = z_und_cont;
edata.E_of_x = E_of_x;
edata.A_of_x = A_of_x;
edata.G_of_x = G_of_x;
edata.J_of_x = J_of_x;
edata.Gamma_of_x = Gamma_of_x;
edata.Iyy_of_x = Iyy_of_x;
edata.Izz_of_x = Izz_of_x;
edata.Ksy_of_x = Ksy_of_x;
edata.Ksz_of_x = Ksz_of_x;
edata.Gamx_of_x = Gamx_of_x;
edata.Gamy_of_x = Gamy_of_x;
edata.Gamz_of_x = Gamz_of_x;
edata.E_of_xi = E_of_xi;
edata.A_of_xi = A_of_xi;
edata.G_of_xi = G_of_xi;
edata.J_of_xi = J_of_xi;
edata.Gamma_of_xi = Gamma_of_xi;
edata.Iyy_of_xi = Iyy_of_xi;
edata.Izz_of_xi = Izz_of_xi;
edata.Ksy_of_xi = Ksy_of_xi;
edata.Ksz_of_xi = Ksz_of_xi;
edata.Gamx_of_xi = Gamx_of_xi;
edata.Gamy_of_xi = Gamy_of_xi;
edata.Gamz_of_xi = Gamz_of_xi;
edata.fx_of_xi = fx_of_xi;
edata.qy_of_xi = qy_of_xi;
edata.qz_of_xi = qz_of_xi;
edata.mx_of_xi = mx_of_xi;
edata.fa_of_xi = fa_of_xi;
edata.ql_of_xi = ql_of_xi;
edata.qt_of_xi = qt_of_xi;
edata.tq_of_xi = tq_of_xi;
edata.bm_of_xi = bm_of_xi;
edata.Px_of_xi = Px_of_xi;
edata.Py_of_xi = Py_of_xi;
edata.Pz_of_xi = Pz_of_xi;
edata.Mx_of_xi = Mx_of_xi;
edata.My_of_xi = My_of_xi;
edata.Mz_of_xi = Mz_of_xi;
edata.Pa_of_xi = Pa_of_xi;
edata.Pl_of_xi = Pl_of_xi;
edata.Pt_of_xi = Pt_of_xi;
edata.Tq_of_xi = Tq_of_xi;
edata.Ml_of_xi = Ml_of_xi;
edata.Mt_of_xi = Mt_of_xi;
edata.Bm_of_xi = Bm_of_xi;
edata.CSx_of_xi = CSx_of_xi;
edata.CSy_of_xi = CSy_of_xi;
edata.CSz_of_xi = CSz_of_xi;
edata.CSu_of_xi = CSu_of_xi;
edata.CSv_of_xi = CSv_of_xi;
edata.CSw_of_xi = CSw_of_xi;
edata.CSt_of_xi = CSt_of_xi;
edata.CSbl_of_xi = CSbl_of_xi;
edata.CSbt_of_xi = CSbt_of_xi;
edata.cfl_of_xi = cfl_of_xi;
edata.cft_of_xi = cft_of_xi;
edata.cf_on_element = cf_on_element;
edata.xi_apx = xi_apx;
edata.xi_apy = xi_apy;
edata.xi_apz = xi_apz;
edata.xi_amx = xi_amx;
edata.xi_amy = xi_amy;
edata.xi_amz = xi_amz;
edata.xi_apa = xi_apa;
edata.xi_apl = xi_apl;
edata.xi_apt = xi_apt;
edata.xi_atq = xi_atq;
edata.xi_aml = xi_aml;
edata.xi_amt = xi_amt;
edata.xi_abm = xi_abm;
edata.xi_akx = xi_akx;
edata.xi_aky = xi_aky;
edata.xi_akz = xi_akz;
edata.xi_aku = xi_aku;
edata.xi_akv = xi_akv;
edata.xi_akw = xi_akw;
edata.xi_akt = xi_akt;
edata.xi_akbl = xi_akbl;
edata.xi_akbt = xi_akbt;

%% Nested functions
    function [xi_a_valid,a] = get_loads_sources_positions(a,Jac,beam,x1)
        xi_a = round(Jac^-1*(a{beam}-x1)-1,4);  % Positions of source in the xi domain    
        valid_ind = (xi_a >= -1 & xi_a <= 1);   % Valid indices (-1 <= xi <= 1) 
        xi_a_valid = xi_a(valid_ind);           % Maintain only the positions with valid indices
        a{beam}(valid_ind) = [];                % Cancel the valid positions so they are not applied to adjecent element as well
    end

end