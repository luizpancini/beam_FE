function [FEMdata,CSx_of_x,CSy_of_x,CSz_of_x,CSu_of_x,CSv_of_x,CSw_of_x,CSt_of_x,CSbl_of_x,CSbt_of_x,Ma_of_x,MaRix_of_x,MaRiy_of_x,MaRiz_of_x] = LNV3DB_loads_sources_funs(FEMdata)      

%% Numerical approximation of concentraded loads and sources - in the local x domain
% delta for approximation of Dirac's delta function
dirac_delta = 1e-4; FEMdata.dirac_delta = dirac_delta;
% Initialize outputs
def_fun_x = cell(FEMdata.N_beams,1); def_fun_x(:) = {@(x) 0};
CSx_of_x = def_fun_x; CSy_of_x = def_fun_x; CSz_of_x = def_fun_x; CSu_of_x = def_fun_x; CSv_of_x = def_fun_x; CSw_of_x = def_fun_x; CSt_of_x = def_fun_x; CSbl_of_x = def_fun_x; CSbt_of_x = def_fun_x; Ma_of_x = def_fun_x; MaRix_of_x = def_fun_x; MaRiy_of_x = def_fun_x; MaRiz_of_x = def_fun_x;  
% Loop over beams
for b=1:FEMdata.N_beams
    % Global frame sources 
    CSx_of_x{b} = add_sources(FEMdata.sources.CSx{b},FEMdata.sources.akx{b},FEMdata.L(b));      % Concentrated x-direction springs at x=akx 
    CSy_of_x{b} = add_sources(FEMdata.sources.CSy{b},FEMdata.sources.aky{b},FEMdata.L(b));      % Concentrated y-direction springs at x=aky
    CSz_of_x{b} = add_sources(FEMdata.sources.CSz{b},FEMdata.sources.akz{b},FEMdata.L(b));      % Concentrated z-direction springs at x=akz
    % Local frame sources
    CSu_of_x{b} = add_sources(FEMdata.sources.CSu{b},FEMdata.sources.aku{b},FEMdata.L(b));      % Concentrated axial springs at x=aku
    CSv_of_x{b} = add_sources(FEMdata.sources.CSv{b},FEMdata.sources.akv{b},FEMdata.L(b));      % Concentrated lateral springs at x=akv
    CSw_of_x{b} = add_sources(FEMdata.sources.CSw{b},FEMdata.sources.akw{b},FEMdata.L(b));      % Concentrated transverse springs at x=akw 
    CSt_of_x{b} = add_sources(FEMdata.sources.CSt{b},FEMdata.sources.akt{b},FEMdata.L(b));      % Concentrated torsional springs at x=akt 
    CSbl_of_x{b} = add_sources(FEMdata.sources.CSbl{b},FEMdata.sources.akbl{b},FEMdata.L(b));   % Concentrated lateral rotation springs at x=akbl
    CSbt_of_x{b} = add_sources(FEMdata.sources.CSbt{b},FEMdata.sources.akbt{b},FEMdata.L(b));   % Concentrated transverse rotation springs at x=akbt
    Ma_of_x{b} = add_sources(FEMdata.sources.Ma{b},FEMdata.sources.ama{b},FEMdata.L(b));        % Concentrated masses' translational inertias at x=ama
    MaRix_of_x{b} = add_sources(FEMdata.sources.MaRix{b},FEMdata.sources.ama{b},FEMdata.L(b));  % Concentrated masses' rotary inertias about the local x-axis at x=ama
    MaRiy_of_x{b} = add_sources(FEMdata.sources.MaRiy{b},FEMdata.sources.ama{b},FEMdata.L(b));  % Concentrated masses' rotary inertias about the local y-axis at x=ama
    MaRiz_of_x{b} = add_sources(FEMdata.sources.MaRiz{b},FEMdata.sources.ama{b},FEMdata.L(b));  % Concentrated masses' rotary inertias about the local z-axis at x=ama
end

% Add to FEM data
FEMdata.sources.CSx_of_x = CSx_of_x;
FEMdata.sources.CSy_of_x = CSy_of_x;
FEMdata.sources.CSz_of_x = CSz_of_x;
FEMdata.sources.CSu_of_x = CSu_of_x;
FEMdata.sources.CSv_of_x = CSv_of_x;
FEMdata.sources.CSw_of_x = CSw_of_x;
FEMdata.sources.CSt_of_x = CSt_of_x;
FEMdata.sources.CSbl_of_x = CSbl_of_x;
FEMdata.sources.Ma_of_x = Ma_of_x;
FEMdata.sources.MaRix_of_x = MaRix_of_x;
FEMdata.sources.MaRiy_of_x = MaRiy_of_x;
FEMdata.sources.MaRiz_of_x = MaRiz_of_x;

%% Nested functions
    function source_fun = add_sources(F,a,L)
        source_fun = @(x) 0;
        for j=1:length(a)
            source = F; loc = a;
            % Adjust load/source location if on beam extremity
            if (loc+dirac_delta) > L
                loc = loc-dirac_delta;
            elseif (loc-dirac_delta) < 0
                loc = loc+dirac_delta;
            end
            % Add contribution
            source_fun = @(x) source_fun(x) + source .*(loc-dirac_delta <= x & x <= loc+dirac_delta);  % Dirac approximation of concentraded source applied at x=loc
        end
    end   

end