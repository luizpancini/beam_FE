function LS3DB_plotter(edata,FEMdata,SOLdata)

% Unpack FEM data
[unit_sys,warp_DOF,N_beams,L,b_alpha,b_beta,b_gamma,BC_nodes,scale,Ne,Px,Py,Pz,Mx,My,Mz,Pa,Pl,Pt,Tq,Ml,Mt,Bm,apx,apy,apz,amx,amy,amz,apa,apl,apt,atq,aml,amt,abm,fa_of_x,tq_of_x,ql_of_x,qt_of_x,fx_of_x,mx_of_x,qy_of_x,qz_of_x,bm_of_x,akx,aky,akz,aku,akv,akw,akt,cfl,cft] = unpack_FEMdata(FEMdata,'plots');

%% Plot style - default settings
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
plot_opt.ax_size = 20;
plot_opt.lw = 1;
plot_opt.ms = 5;
plot_opt.c = jet(N_beams);
xLabel = "$x$ [" + unit_sys.length + "]";

%% Plot deformed structure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]); 
figure1.Position = [488  67  680  695]; 
figure1.Name = "Deformed shape, displacements and internal forces"; 
tabgp = uitabgroup(figure1);
tab0 = uitab('Parent',tabgp,'Title','Deformation','BackgroundColor',[1 1 1]);
axes0 = axes('Parent',tab0,'FontSize',plot_opt.ax_size,'FontName','times new roman');
hold(axes0,'on'); 
axes0.View = [60 25]; 
axis equal;
xlabel("$x$ [" + unit_sys.length + "]",'FontWeight','normal','FontSize',plot_opt.ax_size);
ylabel("$y$ [" + unit_sys.length + "]",'FontWeight','normal','FontSize',plot_opt.ax_size);
zlabel("$z$ [" + unit_sys.length + "]",'FontWeight','normal','FontSize',plot_opt.ax_size);
% Loop over elements
for e=1:Ne
    % Unpack element and solution data
    [e_node_range,beam,R0,x_vec_interp,x0,y0,z0,x_und_nodes,y_und_nodes,z_und_nodes,x_und_cont,y_und_cont,z_und_cont] = unpack_edata(edata,e,'plot_structure');
    [x_def_nodes,y_def_nodes,z_def_nodes,x_def_interp,y_def_interp,z_def_interp] = unpack_SOLdata(SOLdata,e,'plot_structure');
    % Plot structure  
    plot3(x_und_nodes,y_und_nodes,z_und_nodes,'ko--',x_def_nodes,y_def_nodes,z_def_nodes,'bo',x_def_interp,y_def_interp,z_def_interp,'b-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',axes0);
    % Plot BCs
    BC_nodes = LS3DB_draw_BCs_springs_and_masses(BC_nodes,e_node_range,L(beam),R0,b_alpha(beam),b_beta(beam),b_gamma(beam),cfl{beam},cft{beam},akx{beam},aky{beam},akz{beam},aku{beam},akv{beam},akw{beam},akt{beam},[],x0,y0,z0,x_vec_interp,x_und_nodes,y_und_nodes,z_und_nodes);
    % Plot loads
    LS3DB_draw_loads(L(beam),R0,fx_of_x{beam},qy_of_x{beam},qz_of_x{beam},mx_of_x{beam},fa_of_x{beam},ql_of_x{beam},qt_of_x{beam},tq_of_x{beam},bm_of_x{beam},Px{beam},Py{beam},Pz{beam},Mx{beam},My{beam},Mz{beam},Pa{beam},Pl{beam},Pt{beam},Ml{beam},Mt{beam},Tq{beam},Bm{beam},apx{beam},apy{beam},apz{beam},amx{beam},amy{beam},amz{beam},apa{beam},apl{beam},apt{beam},aml{beam},amt{beam},atq{beam},abm{beam},x_vec_interp,x_und_cont,y_und_cont,z_und_cont,x0,y0,z0);
end
p0 = plot3(nan,nan,nan,'ko--',nan,nan,nan,'bo-');
lgd0 = legend(p0,'Undeformed','Deformed');
set(lgd0,'Location','best','FontSize',plot_opt.ax_size*0.6);
annotation(tab0,'textbox',[.01 .93 .15 .06],'String',"Scale = " + scale + "x",'FitBoxToText','on');
 
%% Plot outputs: generalized displacements and internal forces
% Generalized displacements
plot_generalized_output('u','plot_u',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$u$ [" + unit_sys.length + "]",plot_opt); % Axial displacement
plot_generalized_output('v','plot_v',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$v$ [" + unit_sys.length + "]",plot_opt); % Lateral displacement
plot_generalized_output('w','plot_w',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$w$ [" + unit_sys.length + "]",plot_opt); % Transverse displacement
plot_generalized_output('phix','plot_phix',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$\phi_x$ [" + unit_sys.angle + "]",plot_opt); % Twist angle
plot_generalized_output('phiy','plot_phiy',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$\phi_y$ [" + unit_sys.angle + "]",plot_opt); % Transverse rotation angle
plot_generalized_output('phiz','plot_phiz',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$\phi_z$ [" + unit_sys.angle + "]",plot_opt); % Lateral rotation angle
if warp_DOF, plot_generalized_output('dphix','plot_dphix',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$\phi_x^{\prime}$ [" + unit_sys.angle + "/" + unit_sys.length + "]",plot_opt); end % Twist angle derivative
% Generalized forces
plot_generalized_output('Nx','plot_Nx',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$N_x$ [" + unit_sys.force + "]",plot_opt); % Normal force
plot_generalized_output('Vy','plot_Vy',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$V_y$ [" + unit_sys.force + "]",plot_opt); % Lateral shear force
plot_generalized_output('Vz','plot_Vz',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$V_z$ [" + unit_sys.force + "]",plot_opt); % Transverse shear force
plot_generalized_output('Tq','plot_Tq',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$T$ [" + unit_sys.force + "." + unit_sys.length + "]",plot_opt); % Twisting torque
plot_generalized_output('My','plot_My',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$M_y$ [" + unit_sys.force + "." + unit_sys.length + "]",plot_opt); % Transverse bending moment
plot_generalized_output('Mz','plot_Mz',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$M_z$ [" + unit_sys.force + "." + unit_sys.length + "]",plot_opt); % Lateral bending moment
if warp_DOF, plot_generalized_output('Bm','plot_Bm',tabgp,Ne,N_beams,edata,SOLdata,xLabel,"$B_{\omega}$ [" + unit_sys.force + "." + unit_sys.length + "]",plot_opt); end % Bimoment

%% Nested functions

    function plot_generalized_output(Title,SOLdata_opt,tabgp,Ne,N_beams,edata,SOLdata,xLabel,yLabel,opt)      
        tab = uitab('Parent',tabgp,'Title',Title,'BackgroundColor',[1 1 1]);
        ax = axes('Parent',tab,'FontSize',opt.ax_size,'FontName','times new roman');
        hold(ax,'on');
        for i=1:Ne
            [b,x_vec_nodes,x_vec_int] = unpack_edata(edata,i,'plot_gen_outs');
            [out_local_nodes,out_local_interp] = unpack_SOLdata(SOLdata,i,SOLdata_opt);
            plot(x_vec_nodes,out_local_nodes,'o',x_vec_int,out_local_interp,'--','LineWidth',opt.lw,'MarkerSize',opt.ms,'Color',opt.c(b,:),'Parent',ax);
        end
        xlabel(xLabel,'FontWeight','normal','FontSize',opt.ax_size);
        ylabel(yLabel,'FontWeight','normal','FontSize',opt.ax_size);
        p = plot(nan,nan,'ko',nan,nan,'k--');
        lgd = legend(p,'Nodal values','Interpolation');
        set(lgd,'Location','best','FontSize',opt.ax_size*0.6);
        if N_beams > 1
            colormap(opt.c); cb = colorbar; caxis([1 N_beams]); cb.Ticks = 1:1:N_beams; cb.Label.String = 'Beam';
        end       
    end
end