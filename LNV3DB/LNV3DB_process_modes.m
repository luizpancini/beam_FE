function SOLdata = LNV3DB_process_modes(FEMdata,SOLdata,lambda,Wr)

% Unpack FEM data and solution data
[unit_sys,N_modes,Ndof,acc_modes,remove_list] = unpack_FEMdata(FEMdata,'process_modes');

% Check if number of modes sought is achievable
if N_modes > acc_modes
    warning('Number of acessible nodes is smaller than requested');
    FEMdata.N_modes = acc_modes;     % Update FEMdata
end
N_modes = min([N_modes, acc_modes]); % Solution number of modes

% List of nodes of reduced matrices
in_list = 1:Ndof;
in_list(remove_list) = [];

% Get natural frequencies and adjusted eigenvectors
omega = sqrt(real(abs(lambda)));                % Natural frequencies [rad/s]
W = zeros(Ndof,length(in_list));                % Initialize eigenvectors
W(in_list,:) = Wr;                              % Reestablish full eigenvectors
omega = flip(omega); W = flip(W,2);             % Order from lowest to highest frequency
W = W(:,1:N_modes); omega = omega(1:N_modes);   % Get only the number of requested (or maximum possible) nodes 
if unit_sys.freq == "Hz"
    omega = omega/(2*pi);                       % Transformation from rad/s to Hz
end

% Normalize eigenvectors
for m=1:N_modes
    W_max = max(abs(W(:,m)));
    W(:,m) = W(:,m)/W_max;
end

% Add to solution structure
SOLdata.N_modes = N_modes;
SOLdata.omega = omega;
SOLdata.W = W;