function [FEMdata,Px_of_x,Py_of_x,Pz_of_x,Pa_of_x,Pl_of_x,Pt_of_x,Mx_of_x,My_of_x,Mz_of_x,Ml_of_x,Mt_of_x,Tq_of_x,Bm_of_x,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x] = LS3DB_loads_sources_funs(FEMdata)      

%% Numerical approximation of concentraded loads and sources - in the local x domain
% delta for approximation of Dirac's delta function
dirac_delta = 1e-4; FEMdata.dirac_delta = dirac_delta;
% Initialize outputs
def_fun_x = cell(FEMdata.N_beams,1); def_fun_x(:) = {@(x) 0};
Px_of_x = def_fun_x; Py_of_x = def_fun_x; Pz_of_x = def_fun_x; Pa_of_x = def_fun_x; Pl_of_x = def_fun_x; Pt_of_x = def_fun_x; Mx_of_x = def_fun_x; My_of_x = def_fun_x; Mz_of_x = def_fun_x; Tq_of_x = def_fun_x; Ml_of_x = def_fun_x; Mt_of_x = def_fun_x; Bm_of_x = def_fun_x; CSx_of_x = def_fun_x; CSy_of_x = def_fun_x; CSz_of_x = def_fun_x; CSu_of_x = def_fun_x; CSv_of_x = def_fun_x; CSw_of_x = def_fun_x; CSt_of_x = def_fun_x;  
% Loop over beams
for b=1:FEMdata.N_beams
    % Global frame loads/sources
    Px_of_x{b} = add_loads_and_sources(FEMdata.loads.Px{b},FEMdata.loads.apx{b},FEMdata.L(b)); % Concentraded x-direction loads applied at x=apx 
    Py_of_x{b} = add_loads_and_sources(FEMdata.loads.Py{b},FEMdata.loads.apy{b},FEMdata.L(b)); % Concentraded y-direction loads applied at x=apl 
    Pz_of_x{b} = add_loads_and_sources(FEMdata.loads.Pz{b},FEMdata.loads.apz{b},FEMdata.L(b)); % Concentraded z-direction loads applied at x=apt 
    Mx_of_x{b} = add_loads_and_sources(FEMdata.loads.Mx{b},FEMdata.loads.amx{b},FEMdata.L(b)); % Concentraded x-direction moments applied at x=amx 
    My_of_x{b} = add_loads_and_sources(FEMdata.loads.My{b},FEMdata.loads.amy{b},FEMdata.L(b)); % Concentraded y-direction moments applied at x=amy 
    Mz_of_x{b} = add_loads_and_sources(FEMdata.loads.Mz{b},FEMdata.loads.amz{b},FEMdata.L(b)); % Concentraded z-direction moments applied at x=amz 
    CSx_of_x{b} = add_loads_and_sources(FEMdata.sources.CSx{b},FEMdata.sources.akx{b},FEMdata.L(b)); % Concentrated x-direction springs at x=akx 
    CSy_of_x{b} = add_loads_and_sources(FEMdata.sources.CSy{b},FEMdata.sources.aky{b},FEMdata.L(b)); % Concentrated y-direction springs at x=aky
    CSz_of_x{b} = add_loads_and_sources(FEMdata.sources.CSz{b},FEMdata.sources.akz{b},FEMdata.L(b)); % Concentrated z-direction springs at x=akz
    % Local frame loads/sources
    Pa_of_x{b} = add_loads_and_sources(FEMdata.loads.Pa{b},FEMdata.loads.apa{b},FEMdata.L(b)); % Concentrated axial loads applied at x=apa
    Pl_of_x{b} = add_loads_and_sources(FEMdata.loads.Pl{b},FEMdata.loads.apl{b},FEMdata.L(b)); % Concentrated lateral loads applied at x=apl
    Pt_of_x{b} = add_loads_and_sources(FEMdata.loads.Pt{b},FEMdata.loads.apt{b},FEMdata.L(b)); % Concentrated transverse loads applied at x=apt
    Tq_of_x{b} = add_loads_and_sources(FEMdata.loads.Tq{b},FEMdata.loads.atq{b},FEMdata.L(b)); % Concentrated torques applied at x=atq
    Ml_of_x{b} = add_loads_and_sources(FEMdata.loads.Ml{b},FEMdata.loads.aml{b},FEMdata.L(b)); % Concentrated lateral moments applied at x=aml
    Mt_of_x{b} = add_loads_and_sources(FEMdata.loads.Mt{b},FEMdata.loads.amt{b},FEMdata.L(b)); % Concentrated transverse moments applied at x=amt
    Bm_of_x{b} = add_loads_and_sources(FEMdata.loads.Bm{b},FEMdata.loads.abm{b},FEMdata.L(b)); % Concentrated bimoments applied at x=abm
    CSu_of_x{b} = add_loads_and_sources(FEMdata.sources.CSu{b},FEMdata.sources.aku{b},FEMdata.L(b)); % Concentrated axial springs at x=aku
    CSv_of_x{b} = add_loads_and_sources(FEMdata.sources.CSv{b},FEMdata.sources.akv{b},FEMdata.L(b)); % Concentrated lateral springs at x=akv
    CSw_of_x{b} = add_loads_and_sources(FEMdata.sources.CSw{b},FEMdata.sources.akw{b},FEMdata.L(b)); % Concentrated transverse springs at x=akw 
    CSt_of_x{b} = add_loads_and_sources(FEMdata.sources.CSt{b},FEMdata.sources.akt{b},FEMdata.L(b)); % Concentrated torsional springs at x=akt  
end

% Add to FEM data
FEMdata.loads.Px_of_x = Px_of_x;
FEMdata.loads.Py_of_x = Py_of_x;
FEMdata.loads.Pz_of_x = Pz_of_x;
FEMdata.loads.Mx_of_x = Mx_of_x;
FEMdata.loads.My_of_x = My_of_x;
FEMdata.loads.Mz_of_x = Mz_of_x;
FEMdata.loads.Pa_of_x = Pa_of_x;
FEMdata.loads.Pl_of_x = Pl_of_x;
FEMdata.loads.Pt_of_x = Pt_of_x;
FEMdata.loads.Tq_of_x = Tq_of_x;
FEMdata.loads.Ml_of_x = Ml_of_x;
FEMdata.loads.Mt_of_x = Mt_of_x;
FEMdata.loads.Bm_of_x = Bm_of_x;
FEMdata.sources.CSx_of_x = CSx_of_x;
FEMdata.sources.CSy_of_x = CSy_of_x;
FEMdata.sources.CSz_of_x = CSz_of_x;
FEMdata.sources.CSu_of_x = CSu_of_x;
FEMdata.sources.CSv_of_x = CSv_of_x;
FEMdata.sources.CSw_of_x = CSw_of_x;
FEMdata.sources.CSt_of_x = CSt_of_x;

%% Nested functions

    function load_fun = add_loads_and_sources(F,a,L)
        load_fun = @(x) 0;
        for j=1:length(a)
            load = F; loc = a;
            % Adjust load/source location if on beam extremity
            if (loc+dirac_delta) > L
                loc = loc-dirac_delta;
            elseif (loc-dirac_delta) < 0
                loc = loc+dirac_delta;
            end
            % Add contribution
            load_fun = @(x) load_fun(x) + load .*(loc-dirac_delta <= x & x <= loc+dirac_delta);  % Dirac approximation of concentraded load/source applied at x=loc
        end
    end   

end