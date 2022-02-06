function LNV3DB_plotter(edata,FEMdata,SOLdata)

% Unpack FEM data
[unit_sys,L,b_alpha,b_beta,b_gamma,BC_nodes,Ne,N_modes,akx,aky,akz,aku,akv,akw,akt,akbl,akbt,ama,cfl,cft] = unpack_FEMdata(FEMdata,'plots');

%% Plot style - default settings
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
plot_opt.ax_size = 20;
plot_opt.lw = 1;
plot_opt.ms = 5;
plot_opt.c = jet(N_modes);
xLabel = "$x$ [" + unit_sys.length + "]";

%% Plot mode shapes
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]); 
figure1.Position = [488  67  680  695]; 
figure1.Name = "Mode shapes"; 
tabgp = uitabgroup(figure1);
tab0 = uitab('Parent',tabgp,'Title','Mode shapes','BackgroundColor',[1 1 1]);
axes0 = axes('Parent',tab0,'FontSize',plot_opt.ax_size,'FontName','times new roman');
hold(axes0,'on'); 
axes0.View = [60 25]; 
axis equal;
xlabel("$x$ [" + unit_sys.length + "]",'FontWeight','normal','FontSize',plot_opt.ax_size);
ylabel("$y$ [" + unit_sys.length + "]",'FontWeight','normal','FontSize',plot_opt.ax_size);
zlabel("$z$ [" + unit_sys.length + "]",'FontWeight','normal','FontSize',plot_opt.ax_size);
p = gobjects(N_modes+1,1); p_str = cell(N_modes+1,1); p_str{1} = 'Undeformed';
% Loop over modes
for m=1:N_modes
    % Loop over elements
    for e=1:Ne
        % Unpack element and solution data
        [e_node_range,beam,R0,x_vec_interp,x0,y0,z0,x_und_nodes,y_und_nodes,z_und_nodes] = unpack_edata(edata,e,'plot_structure');
        [omega,x_def_nodes,y_def_nodes,z_def_nodes,x_def_interp,y_def_interp,z_def_interp] = unpack_SOLdata(SOLdata,e,'plot_structure',m);
        % Plot BCs
        BC_nodes = LNV3DB_draw_BCs_springs_and_masses(BC_nodes,e_node_range,L(beam),R0,b_alpha(beam),b_beta(beam),b_gamma(beam),cfl{beam},cft{beam},akx{beam},aky{beam},akz{beam},aku{beam},akv{beam},akw{beam},akt{beam},akbl{beam},akbt{beam},ama{beam},x0,y0,z0,x_vec_interp,x_und_nodes,y_und_nodes,z_und_nodes);
        % Plot structure
        if m == 1
            p(1) = plot3(x_und_nodes,y_und_nodes,z_und_nodes,'ko--','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',axes0);
        end
        plot3(x_def_nodes,y_def_nodes,z_def_nodes,'o','Color',plot_opt.c(m,:),'LineStyle','none','MarkerSize',plot_opt.ms,'Parent',axes0);
        p(m+1) = plot3(x_def_interp,y_def_interp,z_def_interp,'-','Color',plot_opt.c(m,:),'LineWidth',plot_opt.lw,'Parent',axes0);
        p_str{m+1} = "Mode " + m + " ($\omega$ = " + round(omega(m),5,'significant') + " " + unit_sys.freq + ")";
    end
end
lgd0 = legend(p,p_str);
set(lgd0,'Location','best','FontSize',plot_opt.ax_size*0.6);

%% Plot outputs: generalized displacements' mode shapes
p(1) = []; p_str(1) = [];
plot_generalized_output('u','plot_u',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$u$ [" + unit_sys.length + "]",plot_opt); % Axial displacement
plot_generalized_output('v','plot_v',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$v$ [" + unit_sys.length + "]",plot_opt); % Lateral displacement
plot_generalized_output('w','plot_w',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$w$ [" + unit_sys.length + "]",plot_opt); % Transverse displacement
plot_generalized_output('phix','plot_phix',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$\phi_x$ [" + unit_sys.angle + "]",plot_opt); % Twist angle
plot_generalized_output('phiy','plot_phiy',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$\phi_y$ [" + unit_sys.angle + "]",plot_opt); % Transverse rotation angle
plot_generalized_output('phiz','plot_phiz',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$\phi_z$ [" + unit_sys.angle + "]",plot_opt); % Lateral rotation angle
plot_generalized_output('dphix','plot_dphix',tabgp,Ne,N_modes,edata,SOLdata,xLabel,"$\phi_x^{\prime}$ [" + unit_sys.angle + "/" + unit_sys.length + "]",plot_opt); % Twist angle derivative

%% Nested functions
    function plot_generalized_output(Title,SOLdata_opt,tabgp,Ne,N_modes,edata,SOLdata,xLabel,yLabel,opt)      
        tab = uitab('Parent',tabgp,'Title',Title,'BackgroundColor',[1 1 1]);
        ax = axes('Parent',tab,'FontSize',opt.ax_size,'FontName','times new roman');
        hold(ax,'on');
        for n=1:N_modes
            for i=1:Ne
                [x_vec_nodes,x_vec_int] = unpack_edata(edata,i,'plot_gen_outs');
                [out_local_nodes,out_local_interp] = unpack_SOLdata(SOLdata,i,SOLdata_opt,n);
                plot(x_vec_nodes,out_local_nodes,'o','LineWidth',opt.lw,'MarkerSize',opt.ms,'Color',opt.c(n,:),'Parent',ax);
                p(n) = plot(x_vec_int,out_local_interp,'--','LineWidth',opt.lw,'MarkerSize',opt.ms,'Color',opt.c(n,:),'Parent',ax);
            end
        end
        xlabel(xLabel,'FontWeight','normal','FontSize',opt.ax_size);
        ylabel(yLabel,'FontWeight','normal','FontSize',opt.ax_size);
        lgd = legend(p,p_str);
        set(lgd,'Location','best','FontSize',opt.ax_size*0.6);    
    end
end