% Pre-processing

%% Create structures with all the FEM, element and solution data
FEMdata = struct();
edata = struct();
SOLdata = struct();

% FEM data
FEMdata.unit_sys.units = unit_sys;
FEMdata.beam_theory = beam_theory;        
FEMdata.N_beams = N_beams;                                
FEMdata.L = L;                     
FEMdata.b_alpha = b_alpha;                       
FEMdata.b_beta = b_beta;                       
FEMdata.b_gamma = b_gamma; 
FEMdata.constitutive_model = constitutive_model;              
FEMdata.Ne_b = Ne_b;          
FEMdata.element_order = element_order;          
FEMdata.elem_connect = elem_connect;        
FEMdata.BC_nodes = BC_nodes; 
FEMdata.RI = RI;
FEMdata.N_modes = N_modes;           

%% If loads/sources not input, set default (empty/zero) values
def_cell = cell(N_beams,1); def_cell(:) = {[]};
def_fun = cell(N_beams,1); def_fun(:) = {@(x) 0};
% Concentraded sources
if ~exist('CSx','var'), CSx = def_cell; else, CSx(end+1:N_beams) = {[]}; CSx = CSx'; end; FEMdata.sources.CSx = CSx;
if ~exist('CSy','var'), CSy = def_cell; else, CSy(end+1:N_beams) = {[]}; CSy = CSy'; end; FEMdata.sources.CSy = CSy;
if ~exist('CSz','var'), CSz = def_cell; else, CSz(end+1:N_beams) = {[]}; CSz = CSz'; end; FEMdata.sources.CSz = CSz;
if ~exist('CSu','var'), CSu = def_cell; else, CSu(end+1:N_beams) = {[]}; CSu = CSu'; end; FEMdata.sources.CSu = CSu;
if ~exist('CSv','var'), CSv = def_cell; else, CSv(end+1:N_beams) = {[]}; CSv = CSv'; end; FEMdata.sources.CSv = CSv;
if ~exist('CSw','var'), CSw = def_cell; else, CSw(end+1:N_beams) = {[]}; CSw = CSw'; end; FEMdata.sources.CSw = CSw;
if ~exist('CSt','var'), CSt = def_cell; else, CSt(end+1:N_beams) = {[]}; CSt = CSt'; end; FEMdata.sources.CSt = CSt;
if ~exist('CSbl','var'), CSbl = def_cell; else, CSbl(end+1:N_beams) = {[]}; CSbl = CSbl'; end; FEMdata.sources.CSbl = CSbl;
if ~exist('CSbt','var'), CSbt = def_cell; else, CSbt(end+1:N_beams) = {[]}; CSbt = CSbt'; end; FEMdata.sources.CSbt = CSbt;
if ~exist('Ma','var'), Ma = def_cell; else, Ma(end+1:N_beams) = {[]}; Ma = Ma'; end; FEMdata.sources.Ma = Ma;
if ~exist('MaRix','var'), MaRix = def_cell; else, MaRix(end+1:N_beams) = {[]}; MaRix = MaRix'; end; FEMdata.sources.MaRix = MaRix;
if ~exist('MaRiy','var'), MaRiy = def_cell; else, MaRiy(end+1:N_beams) = {[]}; MaRiy = MaRiy'; end; FEMdata.sources.MaRiy = MaRiy;
if ~exist('MaRiz','var'), MaRiz = def_cell; else, MaRiz(end+1:N_beams) = {[]}; MaRiz = MaRiz'; end; FEMdata.sources.MaRiz = MaRiz;
% Concentraded sources positions 
if ~exist('akx','var'), akx = def_cell; else, akx(end+1:N_beams) = {[]}; akx = akx'; end; FEMdata.sources.akx = akx;
if ~exist('aky','var'), aky = def_cell; else, aky(end+1:N_beams) = {[]}; aky = aky'; end; FEMdata.sources.aky = aky;
if ~exist('akz','var'), akz = def_cell; else, akz(end+1:N_beams) = {[]}; akz = akz'; end; FEMdata.sources.akz = akz;
if ~exist('aku','var'), aku = def_cell; else, aku(end+1:N_beams) = {[]}; aku = aku'; end; FEMdata.sources.aku = aku;
if ~exist('akv','var'), akv = def_cell; else, akv(end+1:N_beams) = {[]}; akv = akv'; end; FEMdata.sources.akv = akv;
if ~exist('akw','var'), akw = def_cell; else, akw(end+1:N_beams) = {[]}; akw = akw'; end; FEMdata.sources.akw = akw;
if ~exist('akt','var'), akt = def_cell; else, akt(end+1:N_beams) = {[]}; akt = akt'; end; FEMdata.sources.akt = akt;
if ~exist('akbl','var'), akbl = def_cell; else, akbl(end+1:N_beams) = {[]}; akbl = akbl'; end; FEMdata.sources.akbl = akbl;
if ~exist('akbt','var'), akbt = def_cell; else, akbt(end+1:N_beams) = {[]}; akbt = akbt'; end; FEMdata.sources.akbt = akbt;
if ~exist('ama','var'), ama = def_cell; else, ama(end+1:N_beams) = {[]}; ama = ama'; end; FEMdata.sources.ama = ama;
% Distributed sources
if ~exist('cfl_of_x','var'), cfl_of_x = def_fun; else, cfl_of_x(end+1:N_beams) = {@(x) 0}; cfl_of_x(cellfun(@isempty,cfl_of_x)) = {@(x) 0}; cfl_of_x = cfl_of_x'; end; FEMdata.sources.cfl_of_x = cfl_of_x;
if ~exist('cft_of_x','var'), cft_of_x = def_fun; else, cft_of_x(end+1:N_beams) = {@(x) 0}; cft_of_x(cellfun(@isempty,cft_of_x)) = {@(x) 0}; cft_of_x = cft_of_x'; end; FEMdata.sources.cft_of_x = cft_of_x;
    
% Set as column vectors
for b=1:FEMdata.N_beams  
    if ~iscolumn(FEMdata.sources.akx{b}), FEMdata.sources.akx{b} = FEMdata.sources.akx{b}'; end
    if ~iscolumn(FEMdata.sources.aky{b}), FEMdata.sources.aky{b} = FEMdata.sources.aky{b}'; end
    if ~iscolumn(FEMdata.sources.akz{b}), FEMdata.sources.akz{b} = FEMdata.sources.akz{b}'; end
    if ~iscolumn(FEMdata.sources.aku{b}), FEMdata.sources.aku{b} = FEMdata.sources.aku{b}'; end
    if ~iscolumn(FEMdata.sources.akv{b}), FEMdata.sources.akv{b} = FEMdata.sources.akv{b}'; end
    if ~iscolumn(FEMdata.sources.akw{b}), FEMdata.sources.akw{b} = FEMdata.sources.akw{b}'; end
    if ~iscolumn(FEMdata.sources.akt{b}), FEMdata.sources.akt{b} = FEMdata.sources.akt{b}'; end
    if ~iscolumn(FEMdata.sources.akbl{b}), FEMdata.sources.akbl{b} = FEMdata.sources.akbl{b}'; end
    if ~iscolumn(FEMdata.sources.akbt{b}), FEMdata.sources.akbt{b} = FEMdata.sources.akbt{b}'; end
    if ~iscolumn(FEMdata.sources.ama{b}), FEMdata.sources.ama{b} = FEMdata.sources.ama{b}'; end
end

%% Check inputs
% No shear flexibility in the Euler-Bernoulli theory
if beam_theory == "EB", Ksy = @(x) 1e9*ones(N_beams,1); Ksz = @(x) 1e9*ones(N_beams,1); end 
% Check element connectivity
if FEMdata.elem_connect == "unsequenced" 
    if ~exist('elem_nodes','var')
        error('Specify the elem_nodes matrix, with each row containing the nodes of that corresponding element');
    else
        FEMdata.elem_nodes = elem_nodes;
    end
elseif FEMdata.elem_connect ~= "sequenced" && FEMdata.elem_connect ~= "unsequenced"
    error('Specify elem_connect as sequenced or unsequenced');
end
% Check element number vector
if size(FEMdata.Ne_b) ~= FEMdata.N_beams
    error('The size of the Ne_b vector must be equal to N_beams');
end

% Unit system
FEMdata.unit_sys.angle = "rad";
if contains(unit_sys,"m")
    FEMdata.unit_sys.length = "m";    
elseif contains(unit_sys,"in")
    FEMdata.unit_sys.length = "in";
elseif contains(unit_sys,"ft")
    FEMdata.unit_sys.length = "ft";
end
if contains(unit_sys,"rads")
    FEMdata.unit_sys.freq = "rad/s";
elseif contains(unit_sys,"Hz")
    FEMdata.unit_sys.freq = "Hz";
else
    error('Specify frequency unit as Hz or rads');
end

% Check validity
check_sources_inputs(CSx,akx,N_beams,L);
check_sources_inputs(CSy,aky,N_beams,L);
check_sources_inputs(CSz,akz,N_beams,L);
check_sources_inputs(CSu,aku,N_beams,L);
check_sources_inputs(CSv,akv,N_beams,L);
check_sources_inputs(CSw,akw,N_beams,L);
check_sources_inputs(CSt,akt,N_beams,L);
check_sources_inputs(CSbl,akbl,N_beams,L);
check_sources_inputs(CSbt,akbt,N_beams,L);
check_sources_inputs(Ma,ama,N_beams,L);
check_sources_inputs(MaRix,ama,N_beams,L);
check_sources_inputs(MaRiy,ama,N_beams,L);
check_sources_inputs(MaRiz,ama,N_beams,L);

%% Functions
function check_sources_inputs(P,a,N_beams,L)
    if length(P) > length(a)
        error("Load without location specified");
    elseif length(a) > length(P)
        error("Location without load specified");
    end
    for b=1:N_beams
        if length(P{b}) ~= length(a{b})
            error("Load-location pair ill-defined");
        end
        if any(a{b}>L(b))
            error("Load applied outside of beam");
        elseif any(a{b}<0)
            error("Load applied at negative location");
        end
    end
end