function BC_nodes = LS3DB_draw_BCs_springs_and_masses(BC_nodes,e_node_range,L,R0,alpha,beta,gamma,cfl_of_x,cft_of_x,akx,aky,akz,aku,akv,akw,akt,ama,x0b,y0b,z0b,x_interp,x_und_nodes,y_und_nodes,z_und_nodes)

%% Plot BCs
remove_nodes_list = [];
for n=1:length(BC_nodes)
    % BCed nodes and corresponding local nodes for current element (if any)
    node = str2double(BC_nodes{n}(1));
    BC_type = regexpi(BC_nodes{n}(2),'[a-z]','match','once'); if BC_type == "S", BC_type = "SS"; end
    [~,BCed_local_nodes] = find(e_node_range == node);
    if ~isempty(BCed_local_nodes)
        remove_nodes_list = [remove_nodes_list; n];
        % Current local node and its coordinates
        local_node = BCed_local_nodes;
        X = x_und_nodes(local_node);
        Y = y_und_nodes(local_node);
        Z = z_und_nodes(local_node); 
        P = [X; Y; Z]; % Starting point coordinates
        % Plot according to BC type
        if BC_type == "C"   
            plot_clamp(e_node_range,local_node,alpha,beta,gamma,P,L);
        elseif BC_type == "SS" 
            plot_support(e_node_range,local_node,alpha,beta,gamma,P,L);
        else % BC_type = "R": Specified DOF value
            plot_BCed_DOF(BC_nodes{n},X,Y,Z,L);
        end
    end
end
BC_nodes(remove_nodes_list) = []; % Erase that BC so it is not plotted twice

%% Plot springs
% Local elastic foundations
plot_elastic_found(cfl_of_x,x_interp,R0,x0b,y0b,z0b,alpha,beta,gamma,'y',L);              % Lateral elastic foundations
plot_elastic_found(cft_of_x,x_interp,R0,x0b,y0b,z0b,alpha,beta,gamma,'z',L);              % Transverse elastic foundations
% Global frame concentraded springs
plot_spring(akx,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'x',L,'global');                 % x-direction springs
plot_spring(aky,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'y',L,'global','fun',cfl_of_x);  % y-direction springs
plot_spring(akz,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'z',L,'global','fun',cft_of_x);  % z-direction springs
% Local frame concentraded springs
plot_spring(aku,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'x',L,'local');                  % Axial springs
plot_spring(akv,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'y',L,'local','fun',cfl_of_x);   % Lateral springs
plot_spring(akw,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'z',L,'local','fun',cft_of_x);   % Transverse springs
plot_spring(akt,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,'t',L,'local');                  % Torsional springs

%% Plot concentraded masses
plot_mass(ama,x_interp,R0,x0b,y0b,z0b,L);

%% Nested functions

    %% Get rotation matrix
    function R = get_R(alpha,beta,gamma,da,db,dc)
        % Orientation angles
        ca = cos(alpha+da); sa = sin(alpha+da); cb = cos(beta+db); sb = sin(beta+db); cg = cos(gamma+dc); sg = sin(gamma+dc);
        % Rotation matrix
        R = [cb*ca, ca*sg*sb - cg*sa, sg*sa + cg*ca*sb
             cb*sa, cg*ca + sg*sb*sa, cg*sb*sa - ca*sg
               -sb,            cb*sg,            cg*cb]; 
    end

    %% Plot clamp
    function plot_clamp(e_node_range,local_node,alpha,beta,gamma,P0,L)
        % Get angle's increment according to local node position
        if local_node == 1, node_pos = "first"; elseif local_node == length(e_node_range), node_pos = "last"; else, node_pos = "middle"; end
        if node_pos == "first", da = 0; db = pi; dc = 0; elseif node_pos == "last", da = 0; db = 0; dc = 0; else, warning('A middle node has been clamped'); da = 0; db = 0; dc = 0; end 
        % Rotation matrix
        R = get_R(alpha,beta,gamma,da,db,dc);
        % Plot
        plot_clamped_base(P0,R,L);
    end

    %% Plot clamped base
    function plot_clamped_base(P0,R,L,varargin)
        % Handle inputs
        c = 'g'; lw = 1; size = L/20; axis = "x"; size_fac = 1;
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'axis'
                    axis = varargin{2};
                case 'size'
                    size = varargin{2};    
                case 'size_fac'
                    size_fac = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        % Unrotated/untranslated clamped base coordinates
        clx = 1/4*size*[1    0    0    1    0    0   1   0   0   1   0   0   1   0   0   1   0   0   1   0   0    1    0    0    1    0    0    1    0    0    1    0    0    1    0    0    1    0    0    1    0    0    1    0    0    1    0    0    0    0    0];
        cly =     size*[1/2  1/2  1/2  1/2  1/2  1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/4 1/4 1/4 0   0   0   -1/4 -1/4 -1/4 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/4 -1/4 -1/4 0    0    0    1/4  1/4  1/4  1/2  -1/2 -1/2 1/2];
        clz =     size*[-1/2 -1/2 -1/4 -1/4 -1/4 0   0   0   1/4 1/4 1/4 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2  1/2  1/2  1/2  1/2  1/2  1/4  1/4  1/4  0    0    0    -1/4 -1/4 -1/4 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 -1/2 1/2  -1/2 1/2];
        if axis == "x"
            cl_coords = [clx; cly; clz];
        elseif axis == "y"
            cl_coords = [cly; clx; clz];
        elseif axis == "z" || axis == "t"
            cl_coords = [clz; cly; clx];
        end
        cl_coords = size_fac*cl_coords;
        % Translate/rotate clamp
        clamp = P0 + R*cl_coords; 
        % Plot
        plot3(clamp(1,:),clamp(2,:),clamp(3,:),'Color',c,'LineWidth',lw);        
    end
    
    %% Plot simple support
    function plot_support(e_node_range,local_node,alpha,beta,gamma,P0,L)
        % Get angle's increment according to local node position
        if local_node == 1, node_pos = "first"; elseif local_node == length(e_node_range), node_pos = "last"; else, node_pos = "middle"; end
        if node_pos == "first", da = 0; db = pi; dc = 0; elseif node_pos == "last", da = 0; db = 0; dc = 0; else, da = pi/2; db = 0; dc = 0; end
        % Rotation matrix
        R = get_R(alpha,beta,gamma,da,db,dc);
        % Plot the support base
        plot_support_base(P0,R,L);
    end
    
    %% Plot simple support base
    function plot_support_base(P0,R,L,varargin)
        % Handle inputs
        c = 'g'; lw = 1; size = L/20; axis = "x"; 
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'axis'
                    axis = varargin{2};
                case 'size'
                    size = varargin{2};    
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        % Unrotated/untranslated support coordinates
        ssx = size*[1    0 1   0 1    0 1];
        ssy = size*[1/2  0 1/2 0 -1/2 0 -1/2];
        ssz = size*[-1/2 0 1/2 0 1/2  0 -1/2];
        if axis == "x"
            ss_coords = [ssx; ssy; ssz];
        elseif axis == "y"
            ss_coords = [ssy; ssx; ssz];
        elseif axis == "z"
            ss_coords = [ssz; ssx; ssy];
        end
        % Translate/rotate support
        support = P0 + R*ss_coords; 
        % Plot support
        plot3(support(1,:),support(2,:),support(3,:),'Color',c,'LineWidth',lw);
        % Clamped base starting point
        P0 = P0 + size*R(:,1);
        % Plot clamped base
        plot_clamped_base(P0,R,L)
    end
    
    %% Plot BC'ed DOF
    function plot_BCed_DOF(BC_nodes,X,Y,Z,L)
        DOFs_to_BC = str2double(regexpi(BC_nodes(2),'[0-9]','match'))';
        DOFs_values = str2double(BC_nodes(3:end));
        if isempty(DOFs_values) || (max(abs(DOFs_values))) == 0 % Plot only if DOFs have been restricted to zeros
            for j=1:length(DOFs_to_BC)
                DOF = DOFs_to_BC(j);
                switch DOF
                    case 1
                        plot_DOF_cone("x",X,Y,Z,L);
                    case 2
                        plot_DOF_cone("y",X,Y,Z,L);
                    case 3
                        plot_DOF_cone("z",X,Y,Z,L);
                    case 4
                        plot_DOF_cone("x",X,Y,Z,L);
                        plot_DOF_cone("x",X+2*pi*L/10,Y,Z,L);
                    case 5
                        plot_DOF_cone("y",X,Y,Z,L);
                        plot_DOF_cone("y",X,Y+2*pi*L/10,Z,L);
                    case 6
                        plot_DOF_cone("z",X,Y,Z,L);
                        plot_DOF_cone("z",X,Y,Z+2*pi*L/10,L);
                    case 7
                        plot_DOF_cone("x",X,Y,Z,L);
                        plot_DOF_cone("x",X+2*pi*L/10,Y,Z,L);
                        plot_DOF_cone("x",X+4*pi*L/10,Y,Z,L);
                end
            end
        end
    end

    %% Plot BC'ed DOF "cone"
    function plot_DOF_cone(axis,X,Y,Z,L,varargin)
        % Handle inputs
        c = 'g'; n_div = 21; size = L/200;
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'n_div'
                    n_div = varargin{2}; 
                case 'size'
                    size = varargin{2};     
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        % DOF restriction cone plot variables
        r = linspace(0,2*pi,n_div);
        th = linspace(0,2*pi,n_div);
        [H,T] = meshgrid(r,th);
        % Cone coordinates
        x_cone = size/2*H.*cos(T);
        y_cone = size/2*H.*sin(T);
        z_cone = size*H;      
        % Change according to axis
        switch lower(axis)
            case "x"
                xx_cone = z_cone;
                yy_cone = y_cone;
                zz_cone = x_cone;
            case "y"
                xx_cone = y_cone;
                yy_cone = z_cone;
                zz_cone = x_cone;
            case "z"
                xx_cone = x_cone;
                yy_cone = y_cone;
                zz_cone = z_cone;
        end        
        % Plot
        h = surf(X+xx_cone,Y+yy_cone,Z+zz_cone);
        set(h,'FaceColor',c,'EdgeColor','none');       
    end

    %% Plot elastic foundation
    function plot_elastic_found(fun,x_interp,R0,x0b,y0b,z0b,alpha,beta,gamma,axis,L)
        if func2str(fun) == "@(x)0", return; end
        ak = x_interp;
        plot_spring(ak,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,axis,L/2,'local')
    end

    %% Plot spring
    function plot_spring(ak,x_interp,x0b,y0b,z0b,R0,alpha,beta,gamma,axis,L,frame,varargin)
        % Handle inputs
        fun = @(x) 0; lsp = L/10;
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'lsp'
                    lsp = varargin{2};
                case 'fun'
                    fun = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        % Change spring orientation if there is elastic foundation
        da = 0; db = 0; dc = 0;
        if func2str(fun) ~= "@(x)0", da = 0; db = 0; dc = pi; end 
        % Indices of springs in current element
        valid_ind = ak<=max(x_interp)+5e-15; 
        % Spring's starting point
        ak_valid = ak(valid_ind); if isempty(ak_valid), return, end
        P0s = [x0b; y0b; z0b] + R0(:,1)*ak_valid';
        % For y,z-axis and torsional springs, the rotation matrix is readily calculated
        if axis ~= "x" 
            if frame == "global"
                R = get_R(0,0,0,da,db,dc);
            elseif frame == "local"
                R = get_R(alpha,beta,gamma,da,db,dc);
            end
        end
        if axis == "y" || axis == "t"
            dir = R(:,2);
        elseif axis == "z"
            dir = R(:,3);
        end
        % Plot each spring in turn
        for i=1:length(ak_valid)          
            % For x-axis springs, the rotation matrix will be dependable on the closest node        
            if axis == "x"
                if ak_valid(i) <= (x_interp(1)+x_interp(end))/2, closest_node = "first"; else, closest_node = "last"; end 
                if frame == "global"
                    if closest_node == "first", da = pi; end
                    R = get_R(0,0,0,da,db,dc);
                elseif frame == "local"
                    if closest_node == "first", db = pi; end
                    R = get_R(alpha,beta,gamma,da,db,dc);
                end
                dir = R(:,1);
            end
            % Spring
            P0 = P0s(:,i);
            plot_spring_coords(P0,R,axis,lsp)
            % Clamped base
            if axis == "t", lsp = lsp/2; end
            P0 = P0 + lsp*dir;
            plot_clamped_base(P0,R,L,'axis',axis,'size_fac',1/2);
        end
    end

    %% Plot spring coordinates
    function plot_spring_coords(P0,R,axis,lsp,varargin)
        % Handle inputs
        c = 'g'; lw = 1; n_div = 61; r = lsp/2; 
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'n_div'
                    n_div = varargin{2};    
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        % Unrotated/untranslated spring coordinates
        if axis ~= "t"
            spx =  lsp*[0, 0.29, 0.36, 0.43, 0.50, 0.57, 0.64, 0.71, 1];
            spy = zeros(1,length(spx));
            spz = 0.2*lsp*[0, 0, 1, -1, 1, -1, 1, 0, 0];        
            if axis == "x"
                s = [spx; spy; spz];
            elseif axis == "y"
                s = [spy; spx; spz];
            elseif axis == "z"
                s = [spz; spy; spx];
            end
        else
            spiral_fun = @(tt) [zeros(1,length(tt)); r*cos(3*tt).*tt/(2*pi); r*sin(3*tt).*tt/(2*pi)];
            s = spiral_fun(linspace(0,2*pi,n_div));
        end
        % Translate/rotate spring
        spring = P0 + R*s; 
        % Plot
        plot3(spring(1,:),spring(2,:),spring(3,:),'Color',c,'LineWidth',lw);
    end

    %% Plot concentrated masses 
    function plot_mass(ama,x_interp,R0,x0b,y0b,z0b,L,varargin)
        % Handle inputs
        c = 'g'; lw = 1; ms = L/20;
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'ms'
                    ms = varargin{2};    
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        % Indices of masses in current element
        valid_ind = ama<=max(x_interp)+5e-15; 
        % Masses' starting points
        ama_valid = ama(valid_ind); if isempty(ama_valid); return; end
        P0s = [x0b; y0b; z0b] + R0(:,1)*ama_valid';
        % Mass coordinates - always plot in the vertical (z-axis) direction
        mx =  ms*[0, 0,  1/2, 1/2,  -1/2, -1/2, 1/2, -1/2, 1/2,  -1/2, -1/2, 1/2, 1/2,  -1/2, -1/2, 1/2, 1/2, 1/2,   1/2, -1/2, -1/2];
        my =  ms*[0, 0,  1/2, -1/2, -1/2, 1/2,  1/2, -1/2, -1/2, 1/2,  1/2,  1/2, -1/2, -1/2, 1/2,  1/2, 1/2, -1/2, -1/2, -1/2, -1/2];
        mz =  ms*[0, -1, -1,  -1,   -1,   -1,   -1,  -1,   -1,   -1,   -2,   -2,  -2,   -2,   -2,   -2,  -1,  -1,   -2,   -2,   -1];
        m = [mx; my; mz];
        % Plot each mass in turn
        for i=1:length(ama_valid)
            P0 = P0s(:,i);
            mass = P0 + m;
            plot3(mass(1,:),mass(2,:),mass(3,:),'Color',c,'LineWidth',lw);
        end
    end

end