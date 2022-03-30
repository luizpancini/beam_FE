function FEMdata = LNV3DB_constitutive_mat(FEMdata)

% Initialize outputs
D = cell(FEMdata.N_beams,1); D_cf = cell(FEMdata.N_beams,1); D_Mt = cell(FEMdata.N_beams,1); D_Mr = cell(FEMdata.N_beams,1);

% Loop over beams
for i=1:FEMdata.N_beams
    % Stress-strain constitutive matrix
    if lower(FEMdata.constitutive_model) == "iso"
        d = @(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) [E(xi).*A(xi); Ksy(xi).*G(xi).*A(xi); Ksz(xi).*G(xi).*A(xi); G(xi).*J(xi); E(xi).*Iyy(xi); E(xi).*Izz(xi); E(xi).*Gamma(xi)];
        D{i} = @(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz) diag(d(xi,E,A,G,J,Gamma,Iyy,Izz,Ksy,Ksz));                                                                                                                               
    else
        error("Constitutive model not available");
    end
    % Elastic foundation stress-strain constitutive matrix
    d_cf = @(xi,cfl,cft) [cfl(xi); cft(xi)];
    D_cf{i} = @(xi,cfl,cft) diag(d_cf(xi,cfl,cft));
    % Translational inertia "constitutive" matrix
    D_Mt{i} = @(xi,rho,A,Is,e1,e2) rho(xi).*[A(xi) 0              0                           0
                                             0     A(xi)          0               -e1(xi).*A(xi)
                                             0     0              A(xi)            e2(xi).*A(xi)
                                             0     -e1(xi).*A(xi) e2(xi).*A(xi)           Is(xi)];
    % Rotary inertia "constitutive" matrix
    d_Mr = @(xi,rho,Iyy,Izz,Gamma) rho(xi).*[Izz(xi); Iyy(xi); Gamma(xi)];
    D_Mr{i} = @(xi,rho,Iyy,Izz,Gamma) diag(d_Mr(xi,rho,Iyy,Izz,Gamma)); 
end

% Add variables to FEMdata
FEMdata.funs.D = D;
FEMdata.funs.D_cf = D_cf;
FEMdata.funs.D_Mt = D_Mt;
FEMdata.funs.D_Mr = D_Mr;
