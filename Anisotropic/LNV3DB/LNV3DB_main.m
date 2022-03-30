%% Pre-processing
tic
LNV3DB_preprocess;

%% Concentraded sources as functions of local coordinate
[FEMdata,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,CSbl_of_x,CSbt_of_x,Ma_of_x,MaRix_of_x,MaRiy_of_x,MaRiz_of_x] = LNV3DB_loads_sources_funs(FEMdata);

%% Get element and node variables
[FEMdata,edata] = LNV3DB_nodelem_vars(FEMdata,edata,E,G,rho,A,J,Gamma,Is,Iyy,Izz,Ksy,Ksz,e1,e2,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,CSbl_of_x,CSbt_of_x,Ma_of_x,MaRix_of_x,MaRiy_of_x,MaRiz_of_x,cfl_of_x,cft_of_x);

%% Interpolation functions
FEMdata = LNV3DB_interp_funs(FEMdata);

%% Constitutive matrix 
FEMdata = LNV3DB_constitutive_mat(FEMdata);

%% Assemble global matrices
FEMdata = LNV3DB_global_matrices(FEMdata,edata);

%% Boundary conditions 
FEMdata = LNV3DB_apply_BCs(FEMdata);
time_pre = toc;

%% Solve the eigenproblem
tic
% Get eigenvectors and eigenvalues
[Wr,lambda] = eigenshuffle(FEMdata.Mr\FEMdata.Kr); SOLdata.Wr = Wr; SOLdata.lambda = lambda;
% Process modes
SOLdata = LNV3DB_process_modes(FEMdata,SOLdata,lambda,Wr);
time_solve = toc;

%% Calculate outputs
tic
SOLdata = LNV3DB_get_outputs(edata,FEMdata,SOLdata);
 
%% Plot the results
LNV3DB_plotter(edata,FEMdata,SOLdata);
time_post = toc;

% Total computer time
time_total = time_pre + time_solve + time_post;