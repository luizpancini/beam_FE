function FEMdata = LS3DB_constitutive_mat(FEMdata)

%% Constitutive matrices
D = cell(FEMdata.N_beams,1); D_cf = cell(FEMdata.N_beams,1);
for i=1:FEMdata.N_beams
    if lower(FEMdata.constitutive_model) == "iso"
        if FEMdata.warp_DOF
            d = @(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) [E(xi).*A(xi); Ksy(xi).*G(xi).*A(xi); Ksz(xi).*G(xi).*A(xi); G(xi).*J(xi); E(xi).*Iyy(xi); E(xi).*Izz(xi); E(xi).*Gamma(xi)];
        else
            d = @(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) [E(xi).*A(xi); Ksy(xi).*G(xi).*A(xi); Ksz(xi).*G(xi).*A(xi); G(xi).*J(xi); E(xi).*Iyy(xi); E(xi).*Izz(xi)];
        end
        d_cf = @(xi,cfl,cft) [cfl(xi); cft(xi)];
        D{i} = @(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) diag(d(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz)); 
        D_cf{i} = @(xi,cfl,cft) diag(d_cf(xi,cfl,cft));
    else
        error("Constitutive model not available");
    end
end

%% Add variables to FEMdata
FEMdata.funs.D = D;
FEMdata.funs.D_cf = D_cf;

