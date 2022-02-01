function varargout = LS3DB_draw_loads(L,R0,fx_of_x,qy_of_x,qz_of_x,mx_of_x,fa_of_x,ql_of_x,qt_of_x,tq_of_x,bm_of_x,Px,Py,Pz,Mx,My,Mz,Pa,Pl,Pt,Ml,Mt,Tq,Bm,apx,apy,apz,amx,amy,amz,apa,apl,apt,aml,amt,atq,abm,x_interp,x_und,y_und,z_und,x0b,y0b,z0b,varargin)

%% Handle inputs
% Default values
norm_fac = 0.2; n_div = 21; fxh = []; qyh = []; qzh = []; fah = []; qlh = []; qth = []; Pxh = []; Pyh = []; Pzh = []; Pah = []; Plh = []; Pth = []; Mxh = []; Myh = []; Mzh = []; Tqh = []; Bmh = []; Mlh = []; Mth = []; tqh = cell(length(x_interp),1); tqh(:) = {[]}; mxh = tqh; bmh = tqh;
% Load optional arguments
while ~isempty(varargin)
    switch lower(varargin{1})
          case 'norm_fac'
              norm_fac = varargin{2}; % Relative size of arrow to length of beam (L)
          case 'ndiv'
              n_div = varargin{2}; % Size of vector for estimation of maximum distributed loads over current beam
          case 'fxh'
              fxh = varargin{2}; 
          case 'qyh'
              qyh = varargin{2};
          case 'qzh'
              qzh = varargin{2};
          case 'dmxh'
              mxh = varargin{2};
          case 'pxh'
              Pxh = varargin{2};
          case 'pyh'
              Pyh = varargin{2};
          case 'pzh'
              Pzh = varargin{2};    
          case 'mxh'
              Mxh = varargin{2};
          case 'myh'
              Myh = varargin{2};
          case 'mzh'
              Mzh = varargin{2};    
          case 'fah'
              fah = varargin{2};
          case 'qlh'
              qlh = varargin{2};
          case 'qth'
              qth = varargin{2};
          case 'tqh'
              tqh = varargin{2};
          case 'pah'
              Pah = varargin{2};
          case 'plh'
              Plh = varargin{2};
          case 'pth'
              Pth = varargin{2};
          case 'mlh'
              Mlh = varargin{2};
          case 'mth'
              Mth = varargin{2};
          case 'th'
              Tqh = varargin{2};
          otherwise
              error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end
x_dense = linspace(0,L,n_div); % Denser vector for estimation of maximum of distributed loads over current beam

%% Plots
% Global frame distributed loads
draw_dist_forces(fx_of_x,x_dense,x_interp,x_und,y_und,z_und,[1; 0; 0],L,norm_fac/2,'lw',1.5);    % Distributed x-direction forces
draw_dist_forces(qy_of_x,x_dense,x_interp,x_und,y_und,z_und,[0; 1; 0],L,norm_fac/2);             % Distributed y-direction forces
draw_dist_forces(qz_of_x,x_dense,x_interp,x_und,y_und,z_und,[0; 0; 1],L,norm_fac/2);             % Distributed z-direction forces
draw_dist_moments(mx_of_x,x_dense,x_interp,x_und,y_und,z_und,eye(3),L,'h',mxh);                  % Distributed x-direction moments
% Local frame distributed loads
draw_dist_forces(fa_of_x,x_dense,x_interp,x_und,y_und,z_und,R0(:,1),L,norm_fac/2,'lw',1.5);      % Distributed axial forces
draw_dist_forces(ql_of_x,x_dense,x_interp,x_und,y_und,z_und,R0(:,2),L,norm_fac/2);               % Distributed lateral forces
draw_dist_forces(qt_of_x,x_dense,x_interp,x_und,y_und,z_und,R0(:,3),L,norm_fac/2);               % Distributed transverse forces
draw_dist_moments(tq_of_x,x_dense,x_interp,x_und,y_und,z_und,R0,L,'h',tqh);                      % Distributed twisting moments  
% Global frame concentraded loads
draw_con_force(Px,apx,L,x0b,y0b,z0b,x_interp,[1; 0; 0],R0,norm_fac);         % Concentraded x-direction forces
draw_con_force(Py,apy,L,x0b,y0b,z0b,x_interp,[0; 1; 0],R0,norm_fac);         % Concentraded y-direction forces
draw_con_force(Pz,apz,L,x0b,y0b,z0b,x_interp,[0; 0; 1],R0,norm_fac);         % Concentraded z-direction forces
draw_con_moment(Mx,amx,L,x0b,y0b,z0b,x_interp,eye(3),R0,'CircleAxis','x');   % Concentraded x-direction bending moments
draw_con_moment(My,amy,L,x0b,y0b,z0b,x_interp,eye(3),R0,'CircleAxis','y');   % Concentraded y-direction bending moments
draw_con_moment(Mz,amz,L,x0b,y0b,z0b,x_interp,eye(3),R0,'CircleAxis','z');   % Concentraded z-direction bending moments
% Local frame concentraded loads
draw_con_force(Pa,apa,L,x0b,y0b,z0b,x_interp,R0(:,1),R0,norm_fac);           % Concentraded axial forces
draw_con_force(Pl,apl,L,x0b,y0b,z0b,x_interp,R0(:,2),R0,norm_fac);           % Concentraded lateral forces
draw_con_force(Pt,apt,L,x0b,y0b,z0b,x_interp,R0(:,3),R0,norm_fac);           % Concentraded transverse forces
draw_con_moment(Tq,atq,L,x0b,y0b,z0b,x_interp,R0,R0,'CircleAxis','x');       % Concentraded torques
draw_con_moment(Ml,aml,L,x0b,y0b,z0b,x_interp,R0,R0,'CircleAxis','y');       % Concentraded lateral bending moments
draw_con_moment(Mt,amt,L,x0b,y0b,z0b,x_interp,R0,R0,'CircleAxis','z');       % Concentraded transverse bending moments

% Prepare output handles
varargout = {fxh; qyh; qzh; mxh; Pxh; Pyh; Pzh; Mxh; Myh; Mzh; fah; qlh; qth; tqh; bmh; Pah; Plh; Pth; Tqh; Mlh; Mth; Bmh}; 

%% Nested functions

    %% Draw distributed forces
    function h = draw_dist_forces(fun,x_dense,x_interp,x_und,y_und,z_und,dir,L,norm_fac,varargin)      
        % Handle inputs
        c = 'm'; lw = 1; h = []; 
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'h'
                    h = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
             
        if func2str(fun) ~= "@(x)0"              
            % Vectors for quiver
            f_vec_norm = fun(x_interp')/max(abs(fun(x_dense)));  % Normalized vector of distributed loads
            P2 = [x_und, y_und, z_und];
            P1 = P2 - (norm_fac*L*dir*f_vec_norm)';
            DP = P2-P1;          
            % Plot
            if isempty(h)
                h = quiver3(P1(:,1),P1(:,2),P1(:,3),DP(:,1),DP(:,2),DP(:,3),0,'Color',c,'LineWidth',lw,'MaxHeadSize',0.25);
            else
                h.XData = P1(:,1); h.YData = P1(:,2); h.ZData = P1(:,3); h.UData = DP(:,1); h.VData = DP(:,2); h.WData = DP(:,3);
            end
        end
    end

    %% Draw distributed moments
    function h = draw_dist_moments(fun,x_dense,x_interp,x_und,y_und,z_und,R,L,varargin)       
        % Handle inputs
        c = 'm'; h = []; 
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'h'
                    h = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
             
        if func2str(fun) ~= "@(x)0"
            % Get normalized vector of distributed loads
            m_vec_norm = fun(x_interp')/max(abs(fun(x_dense))); % Normalized vector of distributed loads
            if length(m_vec_norm) == 1, m_vec_norm = m_vec_norm*ones(length(x_und),1); end % Adjust for cases where function is input as a constant
            % Loop over interpolation vector
            for i=1:length(x_interp)
                v_M = [x_und(i); y_und(i); z_und(i)];
                h{i} = plot_circ_vec(m_vec_norm(i),v_M,R,L,h{i},'c',c);   
            end
        end
    end

    %% Draw concentraded forces
    function h = draw_con_force(F,a,L,x0b,y0b,z0b,x_interp,dir,R,norm_fac,varargin)
        % Handle inputs
        c = 'r'; lw = 1; h = []; 
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'h'
                    h = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        
        % Get valid loads and position vectors
        valid_ind = a<=max(x_interp)+5e-15;     % Indices of loads in current element
        a_vec = a(valid_ind);                   % Load position vector
        norm_F_vec = F(valid_ind)/max(abs(F));  % Normalized load vector   
        if isempty(a_vec), return; end
        
        % Vectors for quiver
        P2 = [x0b y0b z0b] + a_vec*R(:,1)';     % Endpoints of force vectors
        P1 = P2 - norm_fac*L*norm_F_vec*(dir)'; % Starting points of force vectors
        DP = P2-P1;                             % Force vectors
        
        % Plot
        if isempty(h)
            h = quiver3(P1(:,1),P1(:,2),P1(:,3),DP(:,1),DP(:,2),DP(:,3),0,'Color',c,'LineWidth',lw,'MaxHeadSize',0.25);
        else
            h.XData = P1(:,1); h.YData = P1(:,2); h.ZData = P1(:,3); h.UData = DP(:,1); h.VData = DP(:,2); h.WData = DP(:,3);
        end
    end

    %% Draw concentraded moments
    function h = draw_con_moment(M,a,L,x0b,y0b,z0b,x_interp,R,R0,varargin)        
        % Handle inputs
        c_axis = "x"; h = [];
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'circleaxis'
                    c_axis = varargin{2};
                case 'h'
                    h = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        
        % Get valid loads and position vectors
        valid_ind = a<=max(x_interp)+5e-15;     % Indices of loads in current element
        a_vec = a(valid_ind);                   % Load position vector
        if isempty(a_vec), return; end
        norm_M_vec = M(valid_ind)/max(abs(M));  % Normalized load vector
        v_M = [x0b; y0b; z0b] + R0(:,1)*a_vec';  % Global coordinates of loads positions
        if isempty(h), h = cell(length(a_vec),1); end
        
        % Loop over all loads
        for i=1:length(a_vec)
            h{i} = plot_circ_vec(norm_M_vec(i),v_M(:,i),R,L,h{i},'circleaxis',c_axis);          
        end
    end
    
    %% Plot a circular vector for moments
    function h = plot_circ_vec(M,v_M,R,L,h,varargin)       
        % Handle inputs
        c = 'r'; lw = 1; r = L*abs(M)/20; angle = 9*pi/5; ah = 0.3*r; c_axis = "x"; c_div = 30; 
        % Loading optional arguments
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'c'
                    c = varargin{2};
                case 'lw'
                    lw = varargin{2};
                case 'r'
                    r = varargin{2};
                case 'angle'
                    angle = varargin{2};
                case 'arrowhead'
                    ah = varargin{2};
                case 'circleaxis'
                    c_axis = varargin{2};
                case 'circledivisions'
                    c_div = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
        
        % Circle definition
        if c_axis == "x"
            circ_fun = @(tt) [zeros(1,length(tt)); r*cos(tt); r*sin(tt)];
        elseif c_axis == "y"
            circ_fun = @(tt,x,y,z) [r*cos(tt); zeros(1,length(tt)); r*sin(tt)];
        elseif c_axis == "z"
            circ_fun = @(tt,x,y,z) [r*cos(tt); r*sin(tt); zeros(1,length(tt))];
        else
            error('Choose circle axis as x, y or z');
        end
        circumference = R*circ_fun(linspace(0,angle,c_div)) + v_M;
        
        % Arrow parameters
        dir = sign(M);    % Positive anticlockwise
        if c_axis == "x"
            if dir == 1, tt = angle; else, tt = 0; end
            arrow = [0       0  0;
                -ah/2    0  ah/2;
                -dir*ah  0  -dir*ah];
            Rtt = [1       0        0
                0  cos(tt) -sin(tt);
                0 sin(tt)  cos(tt)];
        elseif c_axis == "y"
            if dir == 1, tt = 0; else, tt = angle; end
            arrow = [ah/2 0 -ah/2;
                0 0 0;
                dir*ah  0 dir*ah];
            Rtt = [cos(tt) 0 -sin(tt);
                0       1        0;
                sin(tt) 0 cos(tt)];
        elseif c_axis == "z"
            if dir == 1, tt = angle; else, tt = 0; end
            arrow = [-ah/2 0 ah/2;
                -dir*ah  0 -dir*ah;
                0 0 0];
            Rtt = [cos(tt) -sin(tt) 0;
                sin(tt)  cos(tt) 0;
                0        0       1];
        end
        arrowhead = R*(circ_fun(tt) + Rtt*arrow) + v_M;
        
        % Plot
        if isempty(h)
            h1 = plot3(circumference(1,:),circumference(2,:),circumference(3,:),'Color',c,'LineWidth',lw);
            h2 = plot3(arrowhead(1,:),arrowhead(2,:),arrowhead(3,:),'Color',c,'LineWidth',lw);
        else
            h1.XData = circumference(1,:); h1.YData = circumference(2,:); h1.ZData = circumference(3,:);
            h2.XData = arrowhead(1,:); h2.YData = arrowhead(2,:); h2.ZData = arrowhead(3,:);
        end
        h = [h1; h2];
    end
end