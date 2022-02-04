function FEMdata = LS3DB_apply_BCs(FEMdata)

% Unpack FEM data
[BC_nodes,N_nodes,ndof,elem_nodes,EDOFs,DOF_u,DOF_phix,DOF_bendl,DOF_bendt] = unpack_FEMdata(FEMdata,'bcs');

%% Apply BCs
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
        DOFs_values = zeros(length(BCed_dofs),1);
    elseif BC_type == "SS"                                                  % Simple support: restrict all displacements and twist rotation 
        DOFs_to_BC = [DOF_u(1) DOF_bendl(1) DOF_bendt(1) DOF_phix(1)];
        BCed_dofs = node_DOFs(DOFs_to_BC); 
        DOFs_values = zeros(length(BCed_dofs),1);
    elseif BC_type == "R"                                                   % Restrict DOF # to specified values
        DOFs_to_BC = str2double(regexpi(BC_nodes{i}(2),'[0-9]','match'))';
        BCed_dofs =  node_DOFs(DOFs_to_BC);
        DOFs_values = (str2double(BC_nodes{i}(3:end)))';
        if isempty(DOFs_values)                                                 
            DOFs_values = zeros(length(DOFs_to_BC),1);                      % If DOF values are not explicitly specified, assume they are zero
        end
        if length(DOFs_values) ~= length(DOFs_to_BC)
            error("DOFs to be restricted and respective values do not have the same size");
        end
    else
        error("Unknown BC type: " + BC_type);
    end
    % Adjust stiffness matrix by changing lines corresponding to BCed DOFs
    zero_lines = zeros(length(BCed_dofs),length(FEMdata.F));
    d = diag(FEMdata.K(BCed_dofs,BCed_dofs));
    FEMdata.K(BCed_dofs,:) = zero_lines; 
    original_diag = d.*eye(length(BCed_dofs));
    FEMdata.K(BCed_dofs,BCed_dofs) = original_diag;
    % Force vector = specified displacement * diagonal stiffnesses:
    FEMdata.F(BCed_dofs) = DOFs_values.*d;
end

