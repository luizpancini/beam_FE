%% Pre-processing
tic
LS3DB_preprocess;

%% Concentraded loads and sources as functions of local coordinate
[FEMdata,Px_of_x,Py_of_x,Pz_of_x,Pa_of_x,Pl_of_x,Pt_of_x,Mx_of_x,My_of_x,Mz_of_x,Ml_of_x,Mt_of_x,Tq_of_x,Bm_of_x,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x] = LS3DB_loads_sources_funs(FEMdata);

%% Get element and node variables
[FEMdata,edata] = LS3DB_nodelem_vars(FEMdata,edata,E,G,A,J,Gamma,Iyy,Izz,Ksy,Ksz,fx_of_x,mx_of_x,qy_of_x,qz_of_x,fa_of_x,ql_of_x,qt_of_x,tq_of_x,bm_of_x,Px_of_x,Py_of_x,Pz_of_x,Pa_of_x,Pl_of_x,Pt_of_x,Mx_of_x,My_of_x,Mz_of_x,Ml_of_x,Mt_of_x,Tq_of_x,Bm_of_x,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,cfl_of_x,cft_of_x);

%% Interpolation functions
FEMdata = LS3DB_interp_funs(FEMdata);

%% Constitutive matrix 
FEMdata = LS3DB_constitutive_mat(FEMdata);

%% Assemble global matrices
clearvars -except FEMdata edata
FEMdata = LS3DB_global_matrices(FEMdata,edata);

%% Boundary conditions 
FEMdata = LS3DB_apply_BCs(FEMdata);
time_pre = toc;

%% Solve the linear system 
tic
warning('off')
SOLdata.u_global = FEMdata.K\FEMdata.F;
warning('on')
time_solve = toc;

%% Calculate outputs
tic
SOLdata = LS3DB_get_outputs(edata,FEMdata,SOLdata);
 
%% Plot the results
LS3DB_plotter(edata,FEMdata,SOLdata);
time_post = toc;

% Total computer time
time_total = time_pre + time_solve + time_post;