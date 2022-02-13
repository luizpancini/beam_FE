% Pre-processing

%% Create structures with all the FEM, element and solution data
FEMdata = struct();
edata = struct();
SOLdata = struct();

% FEM data
FEMdata.unit_sys.units = unit_sys;
FEMdata.beam_theory = beam_theory;
FEMdata.warp_DOF = warp_DOF;         
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
FEMdata.scale = scale;           
if exist('elem_nodes','var'), FEMdata.elem_nodes = elem_nodes; end

%% If loads/sources not input, set default (empty/zero) values
def_cell = cell(N_beams,1); def_cell(:) = {[]};
def_fun = cell(N_beams,1); def_fun(:) = {@(x) 0};
% Concentraded loads
if ~exist('Px','var'), Px = def_cell; else, Px(end+1:N_beams) = {[]}; Px = Px'; end; FEMdata.loads.Px = Px;
if ~exist('Py','var'), Py = def_cell; else, Py(end+1:N_beams) = {[]}; Py = Py'; end; FEMdata.loads.Py = Py;
if ~exist('Pz','var'), Pz = def_cell; else, Pz(end+1:N_beams) = {[]}; Pz = Pz'; end; FEMdata.loads.Pz = Pz;
if ~exist('Mx','var'), Mx = def_cell; else, Mx(end+1:N_beams) = {[]}; Mx = Mx'; end; FEMdata.loads.Mx = Mx;
if ~exist('My','var'), My = def_cell; else, My(end+1:N_beams) = {[]}; My = My'; end; FEMdata.loads.My = My;
if ~exist('Mz','var'), Mz = def_cell; else, Mz(end+1:N_beams) = {[]}; Mz = Mz'; end; FEMdata.loads.Mz = Mz;
if ~exist('Pa','var'), Pa = def_cell; else, Pa(end+1:N_beams) = {[]}; Pa = Pa'; end; FEMdata.loads.Pa = Pa;
if ~exist('Pl','var'), Pl = def_cell; else, Pl(end+1:N_beams) = {[]}; Pl = Pl'; end; FEMdata.loads.Pl = Pl;
if ~exist('Pt','var'), Pt = def_cell; else, Pt(end+1:N_beams) = {[]}; Pt = Pt'; end; FEMdata.loads.Pt = Pt;
if ~exist('Tq','var'), Tq = def_cell; else, Tq(end+1:N_beams) = {[]}; Tq = Tq'; end; FEMdata.loads.Tq = Tq;
if ~exist('Ml','var'), Ml = def_cell; else, Ml(end+1:N_beams) = {[]}; Ml = Ml'; end; FEMdata.loads.Ml = Ml;
if ~exist('Mt','var'), Mt = def_cell; else, Mt(end+1:N_beams) = {[]}; Mt = Mt'; end; FEMdata.loads.Mt = Mt;
if ~exist('Bm','var'), Bm = def_cell; else, Bm(end+1:N_beams) = {[]}; Bm = Bm'; end; FEMdata.loads.Bm = Bm;
% Concentraded loads positions
if ~exist('apx','var'), apx = def_cell; else, apx(end+1:N_beams) = {[]}; apx = apx'; end; FEMdata.loads.apx = apx;
if ~exist('apy','var'), apy = def_cell; else, apy(end+1:N_beams) = {[]}; apy = apy'; end; FEMdata.loads.apy = apy;
if ~exist('apz','var'), apz = def_cell; else, apz(end+1:N_beams) = {[]}; apz = apz'; end; FEMdata.loads.apz = apz;
if ~exist('amx','var'), amx = def_cell; else, amx(end+1:N_beams) = {[]}; amx = amx'; end; FEMdata.loads.amx = amx;
if ~exist('amy','var'), amy = def_cell; else, amy(end+1:N_beams) = {[]}; amy = amy'; end; FEMdata.loads.amy = amy;
if ~exist('amz','var'), amz = def_cell; else, amz(end+1:N_beams) = {[]}; amz = amz'; end; FEMdata.loads.amz = amz;
if ~exist('apa','var'), apa = def_cell; else, apa(end+1:N_beams) = {[]}; apa = apa'; end; FEMdata.loads.apa = apa;
if ~exist('apl','var'), apl = def_cell; else, apl(end+1:N_beams) = {[]}; apl = apl'; end; FEMdata.loads.apl = apl;
if ~exist('apt','var'), apt = def_cell; else, apt(end+1:N_beams) = {[]}; apt = apt'; end; FEMdata.loads.apt = apt;
if ~exist('atq','var'), atq = def_cell; else, atq(end+1:N_beams) = {[]}; atq = atq'; end; FEMdata.loads.atq = atq;
if ~exist('aml','var'), aml = def_cell; else, aml(end+1:N_beams) = {[]}; aml = aml'; end; FEMdata.loads.aml = aml;
if ~exist('amt','var'), amt = def_cell; else, amt(end+1:N_beams) = {[]}; amt = amt'; end; FEMdata.loads.amt = amt;
if ~exist('abm','var'), abm = def_cell; else, abm(end+1:N_beams) = {[]}; abm = abm'; end; FEMdata.loads.abm = abm;
% Distributed loads
if ~exist('fa_of_x','var'), fa_of_x = def_fun; else, fa_of_x(end+1:N_beams) = {@(x) 0}; fa_of_x(cellfun(@isempty,fa_of_x)) = {@(x) 0}; fa_of_x = fa_of_x'; end; FEMdata.loads.fa_of_x = fa_of_x;
if ~exist('tq_of_x','var'), tq_of_x = def_fun; else, tq_of_x(end+1:N_beams) = {@(x) 0}; tq_of_x(cellfun(@isempty,tq_of_x)) = {@(x) 0}; tq_of_x = tq_of_x'; end; FEMdata.loads.tq_of_x = tq_of_x;
if ~exist('ql_of_x','var'), ql_of_x = def_fun; else, ql_of_x(end+1:N_beams) = {@(x) 0}; ql_of_x(cellfun(@isempty,ql_of_x)) = {@(x) 0}; ql_of_x = ql_of_x'; end; FEMdata.loads.ql_of_x = ql_of_x;
if ~exist('qt_of_x','var'), qt_of_x = def_fun; else, qt_of_x(end+1:N_beams) = {@(x) 0}; qt_of_x(cellfun(@isempty,qt_of_x)) = {@(x) 0}; qt_of_x = qt_of_x'; end; FEMdata.loads.qt_of_x = qt_of_x;
if ~exist('fx_of_x','var'), fx_of_x = def_fun; else, fx_of_x(end+1:N_beams) = {@(x) 0}; fx_of_x(cellfun(@isempty,fx_of_x)) = {@(x) 0}; fx_of_x = fx_of_x'; end; FEMdata.loads.fx_of_x = fx_of_x;
if ~exist('mx_of_x','var'), mx_of_x = def_fun; else, mx_of_x(end+1:N_beams) = {@(x) 0}; mx_of_x(cellfun(@isempty,mx_of_x)) = {@(x) 0}; mx_of_x = mx_of_x'; end; FEMdata.loads.mx_of_x = mx_of_x;
if ~exist('qy_of_x','var'), qy_of_x = def_fun; else, qy_of_x(end+1:N_beams) = {@(x) 0}; qy_of_x(cellfun(@isempty,qy_of_x)) = {@(x) 0}; qy_of_x = qy_of_x'; end; FEMdata.loads.qy_of_x = qy_of_x;
if ~exist('qz_of_x','var'), qz_of_x = def_fun; else, qz_of_x(end+1:N_beams) = {@(x) 0}; qz_of_x(cellfun(@isempty,qz_of_x)) = {@(x) 0}; qz_of_x = qz_of_x'; end; FEMdata.loads.qz_of_x = qz_of_x;
if ~exist('bm_of_x','var'), bm_of_x = def_fun; else, bm_of_x(end+1:N_beams) = {@(x) 0}; bm_of_x(cellfun(@isempty,bm_of_x)) = {@(x) 0}; bm_of_x = bm_of_x'; end; FEMdata.loads.bm_of_x = bm_of_x;
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
% Distributed sources
if ~exist('cfl_of_x','var'), cfl_of_x = def_fun; else, cfl_of_x(end+1:N_beams) = {@(x) 0}; cfl_of_x(cellfun(@isempty,cfl_of_x)) = {@(x) 0}; cfl_of_x = cfl_of_x'; end; FEMdata.sources.cfl_of_x = cfl_of_x;
if ~exist('cft_of_x','var'), cft_of_x = def_fun; else, cft_of_x(end+1:N_beams) = {@(x) 0}; cft_of_x(cellfun(@isempty,cft_of_x)) = {@(x) 0}; cft_of_x = cft_of_x'; end; FEMdata.sources.cft_of_x = cft_of_x;
    
% Set as column vectors
for b=1:FEMdata.N_beams
    if ~iscolumn(FEMdata.loads.Px{b}), FEMdata.loads.Px{b} = FEMdata.loads.Px{b}'; end
    if ~iscolumn(FEMdata.loads.Py{b}), FEMdata.loads.Py{b} = FEMdata.loads.Py{b}'; end
    if ~iscolumn(FEMdata.loads.Pz{b}), FEMdata.loads.Pz{b} = FEMdata.loads.Pz{b}'; end
    if ~iscolumn(FEMdata.loads.Mx{b}), FEMdata.loads.Mx{b} = FEMdata.loads.Mx{b}'; end
    if ~iscolumn(FEMdata.loads.My{b}), FEMdata.loads.My{b} = FEMdata.loads.My{b}'; end
    if ~iscolumn(FEMdata.loads.Mz{b}), FEMdata.loads.Mz{b} = FEMdata.loads.Mz{b}'; end
    if ~iscolumn(FEMdata.loads.Pa{b}), FEMdata.loads.Pa{b} = FEMdata.loads.Pa{b}'; end
    if ~iscolumn(FEMdata.loads.Pl{b}), FEMdata.loads.Pl{b} = FEMdata.loads.Pl{b}'; end
    if ~iscolumn(FEMdata.loads.Pt{b}), FEMdata.loads.Pt{b} = FEMdata.loads.Pt{b}'; end
    if ~iscolumn(FEMdata.loads.Tq{b}), FEMdata.loads.Tq{b} = FEMdata.loads.Tq{b}'; end
    if ~iscolumn(FEMdata.loads.Ml{b}), FEMdata.loads.Ml{b} = FEMdata.loads.Ml{b}'; end
    if ~iscolumn(FEMdata.loads.Mt{b}), FEMdata.loads.Mt{b} = FEMdata.loads.Mt{b}'; end
    if ~iscolumn(FEMdata.loads.Bm{b}), FEMdata.loads.Bm{b} = FEMdata.loads.Bm{b}'; end
    if ~iscolumn(FEMdata.loads.apx{b}), FEMdata.loads.apx{b} = FEMdata.loads.apx{b}'; end
    if ~iscolumn(FEMdata.loads.apy{b}), FEMdata.loads.apy{b} = FEMdata.loads.apy{b}'; end
    if ~iscolumn(FEMdata.loads.apz{b}), FEMdata.loads.apz{b} = FEMdata.loads.apz{b}'; end
    if ~iscolumn(FEMdata.loads.amx{b}), FEMdata.loads.amx{b} = FEMdata.loads.amx{b}'; end
    if ~iscolumn(FEMdata.loads.amy{b}), FEMdata.loads.amy{b} = FEMdata.loads.amy{b}'; end
    if ~iscolumn(FEMdata.loads.amz{b}), FEMdata.loads.amz{b} = FEMdata.loads.amz{b}'; end
    if ~iscolumn(FEMdata.loads.apa{b}), FEMdata.loads.apa{b} = FEMdata.loads.apa{b}'; end
    if ~iscolumn(FEMdata.loads.apl{b}), FEMdata.loads.apl{b} = FEMdata.loads.apl{b}'; end
    if ~iscolumn(FEMdata.loads.apt{b}), FEMdata.loads.apt{b} = FEMdata.loads.apt{b}'; end
    if ~iscolumn(FEMdata.loads.atq{b}), FEMdata.loads.atq{b} = FEMdata.loads.atq{b}'; end
    if ~iscolumn(FEMdata.loads.aml{b}), FEMdata.loads.aml{b} = FEMdata.loads.aml{b}'; end
    if ~iscolumn(FEMdata.loads.amt{b}), FEMdata.loads.amt{b} = FEMdata.loads.amt{b}'; end
    if ~iscolumn(FEMdata.loads.abm{b}), FEMdata.loads.abm{b} = FEMdata.loads.abm{b}'; end    
    if ~iscolumn(FEMdata.sources.akx{b}), FEMdata.sources.akx{b} = FEMdata.sources.akx{b}'; end
    if ~iscolumn(FEMdata.sources.aky{b}), FEMdata.sources.aky{b} = FEMdata.sources.aky{b}'; end
    if ~iscolumn(FEMdata.sources.akz{b}), FEMdata.sources.akz{b} = FEMdata.sources.akz{b}'; end
    if ~iscolumn(FEMdata.sources.aku{b}), FEMdata.sources.aku{b} = FEMdata.sources.aku{b}'; end
    if ~iscolumn(FEMdata.sources.akv{b}), FEMdata.sources.akv{b} = FEMdata.sources.akv{b}'; end
    if ~iscolumn(FEMdata.sources.akw{b}), FEMdata.sources.akw{b} = FEMdata.sources.akw{b}'; end
    if ~iscolumn(FEMdata.sources.akt{b}), FEMdata.sources.akt{b} = FEMdata.sources.akt{b}'; end
    if ~iscolumn(FEMdata.sources.akbl{b}), FEMdata.sources.akbl{b} = FEMdata.sources.akbl{b}'; end
    if ~iscolumn(FEMdata.sources.akbt{b}), FEMdata.sources.akbt{b} = FEMdata.sources.akbt{b}'; end
end

%% Check inputs
% No shear flexibility in the Euler-Bernoulli theory
if beam_theory == "EB", Ksy = @(x) 1e9*ones(N_beams,1); Ksz = @(x) 1e9*ones(N_beams,1); end 
% Check element connectivity
if FEMdata.elem_connect == "unsequenced" 
    if ~exist('elem_nodes','var')
        error('Specify the elem_nodes matrix, with each row containing the nodes of that corresponding element');
    else
        edata.elem_nodes = elem_nodes;
    end
elseif FEMdata.elem_connect ~= "sequenced" && FEMdata.elem_connect ~= "unsequenced"
    error('Specify elem_connect as sequenced or unsequenced');
elseif ~exist('elem_nodes','var')
    FEMdata.elem_nodes = [];
end
% Check element number vector
if size(FEMdata.Ne_b) ~= FEMdata.N_beams
    error('The size of the Ne_b vector must be equal to N_beams');
end
% Check GJ/EGam ratio
if warp_DOF 
    if element_order == "linear", max_ratio = 1; tip = 1; elseif element_order == "quadratic", max_ratio = 1e1; tip = 0; end
    ratios = G(0).*J(0)./(E(0).*Gamma(0))./Ne_b;
    high_ratios_ind = round(ratios,1) > max_ratio;
    if any(high_ratios_ind) % If the ratio GJ/EGam is very high
        if tip, tip = ", or increasing the element order"; else, tip = ""; end
        warning("Ratio GJ/E*Gam very high for beam(s) " + num2str(high_ratios_ind) + " - results for internal torque may be innaccurate. Consider increasing the number of elements, setting warp_DOF to zero" + tip);
    end
end    

% Unit system
if unit_sys == "SI-N.m"
    FEMdata.unit_sys.length = "m";
    FEMdata.unit_sys.angle = "rad";
    FEMdata.unit_sys.force = "N";
elseif unit_sys == "SI-kN.m"
    FEMdata.unit_sys.length = "m";
    FEMdata.unit_sys.angle = "rad";
    FEMdata.unit_sys.force = "kN";
elseif unit_sys == "E-lb.in"
    FEMdata.unit_sys.length = "in";
    FEMdata.unit_sys.angle = "rad";
    FEMdata.unit_sys.force = "lb";
elseif unit_sys == "E-lb.ft"
    FEMdata.unit_sys.length = "ft";
    FEMdata.unit_sys.angle = "rad";
    FEMdata.unit_sys.force = "lb";
end

% Check validity
check_loads_sources_inputs(Px,apx,N_beams,L);
check_loads_sources_inputs(Py,apy,N_beams,L);
check_loads_sources_inputs(Pz,apz,N_beams,L);
check_loads_sources_inputs(Mx,amx,N_beams,L);
check_loads_sources_inputs(My,amy,N_beams,L);
check_loads_sources_inputs(Mz,amz,N_beams,L);
check_loads_sources_inputs(Pa,apa,N_beams,L);
check_loads_sources_inputs(Pl,apl,N_beams,L);
check_loads_sources_inputs(Pt,apt,N_beams,L);
check_loads_sources_inputs(Tq,atq,N_beams,L);
check_loads_sources_inputs(Ml,aml,N_beams,L);
check_loads_sources_inputs(Mt,amt,N_beams,L);
check_loads_sources_inputs(Bm,abm,N_beams,L);
check_loads_sources_inputs(CSx,akx,N_beams,L);
check_loads_sources_inputs(CSy,aky,N_beams,L);
check_loads_sources_inputs(CSz,akz,N_beams,L);
check_loads_sources_inputs(CSu,aku,N_beams,L);
check_loads_sources_inputs(CSv,akv,N_beams,L);
check_loads_sources_inputs(CSw,akw,N_beams,L);
check_loads_sources_inputs(CSt,akt,N_beams,L);
check_loads_sources_inputs(CSbl,akbl,N_beams,L);
check_loads_sources_inputs(CSbt,akbt,N_beams,L);

%% Functions
function check_loads_sources_inputs(P,a,N_beams,L)
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