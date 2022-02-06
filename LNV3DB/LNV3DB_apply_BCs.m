function FEMdata = LNV3DB_apply_BCs(FEMdata)

% Unpack FEM data
[BC_nodes,N_nodes,ndof,elem_nodes,EDOFs,DOF_u,DOF_v,DOF_w,DOF_phix] = unpack_FEMdata(FEMdata,'bcs');

%% Apply BCs
remove_list = [];  
% Loop over BCs
for i=1:length(BC_nodes)
    % Get current node and its BCs
    node = str2double(BC_nodes{i}(1));
    BC_type = regexpi(BC_nodes{i}(2),'[a-z]','match','once'); if BC_type == "S", BC_type = "SS"; end
    if node > N_nodes
        warning("BC applied to inexistent node - model may be underconstrained");
        continue
    end
    % Get current node's DOFs
    e = find(any(elem_nodes==node,2),1);                                    % First element that contains that node
    local_node = find(any(elem_nodes(e,:)==node,1),1);                      % Corresponding local node in that element
    e_nodes = elem_nodes(e,:);                                              % Element's node range
    node_DOFs = ndof*e_nodes(local_node)-(ndof-1):ndof*e_nodes(local_node); % Node's global DOFs
    % Collect BC'ed DOFs        
    if BC_type == "C"                                                       % Clamp: restrict all DOFs
        DOFs_to_BC = EDOFs{1};
        BCed_dofs = node_DOFs(DOFs_to_BC);                                         
    elseif BC_type == "SS"                                                  % Simple support: restrict all displacements and rotation
        DOFs_to_BC = [DOF_u(1) DOF_v(1) DOF_w(1) DOF_phix(1)];
        BCed_dofs = node_DOFs(DOFs_to_BC); 
    elseif BC_type == "R"                                                   % Restrict DOF # to specified values
        DOFs_to_BC = str2double(regexpi(BC_nodes{i}(2),'[0-9]','match'))';
        BCed_dofs =  node_DOFs(DOFs_to_BC);
    else
        error("Unknown BC type: " + BC_type);
    end
    remove_list = [remove_list; BCed_dofs'];
end
FEMdata.remove_list = unique(remove_list);          % List of all BC'ed global DOFs
d_modes = length(remove_list);                      % Number of inaccessible nodes
FEMdata.acc_modes = length(FEMdata.K) - d_modes;    % Number of accessible nodes

% Reduced matrices for computation of eigenvalues - remove BC'ed DOFs
FEMdata.Kr = FEMdata.K; FEMdata.Mr = FEMdata.M;                 
FEMdata.Kr(remove_list,:) = []; FEMdata.Kr(:,remove_list) = [];
FEMdata.Mr(remove_list,:) = []; FEMdata.Mr(:,remove_list) = [];

