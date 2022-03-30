function FEMdata = LS3DB_constitutive_mat(FEMdata)

%% Constitutive matrices
D = cell(FEMdata.N_beams,1); D_cf = cell(FEMdata.N_beams,1);
for b=1:FEMdata.N_beams
    if FEMdata.props.constitutive_model{b} == "iso"
        if FEMdata.warp_DOF
            d = @(xi,props) [props.E_of_xi(xi).*props.A_of_xi(xi); props.Ksy_of_xi(xi).*props.G_of_xi(xi).*props.A_of_xi(xi); props.Ksz_of_xi(xi).*props.G_of_xi(xi).*props.A_of_xi(xi); props.G_of_xi(xi).*props.J_of_xi(xi); props.E_of_xi(xi).*props.Iyy_of_xi(xi); props.E_of_xi(xi).*props.Izz_of_xi(xi); props.E_of_xi(xi).*props.Gamma_of_xi(xi)];
        else
            d = @(xi,props) [props.E_of_xi(xi).*props.A_of_xi(xi); props.Ksy_of_xi(xi).*props.G_of_xi(xi).*props.A_of_xi(xi); props.Ksz_of_xi(xi).*props.G_of_xi(xi).*props.A_of_xi(xi); props.G_of_xi(xi).*props.J_of_xi(xi); props.E_of_xi(xi).*props.Iyy_of_xi(xi); props.E_of_xi(xi).*props.Izz_of_xi(xi)];
        end
        d_cf = @(xi,cfl,cft) [cfl(xi); cft(xi)];
        D{b} = @(xi,props) diag(d(xi,props)); 
        D_cf{b} = @(xi,cfl,cft) diag(d_cf(xi,cfl,cft));
    elseif FEMdata.props.constitutive_model{b} == "aniso"
        D{b} = @(xi,props) props.C_aniso_of_xi(xi);
        d_cf = @(xi,cfl,cft) [cfl(xi); cft(xi)];
        D_cf{b} = @(xi,cfl,cft) diag(d_cf(xi,cfl,cft));
    else
        error("Constitutive model not available");
    end
end

%% Add variables to FEMdata
FEMdata.funs.D = D;
FEMdata.funs.D_cf = D_cf;

