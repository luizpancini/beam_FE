function I = gauss_legendre(func,a,b,pts,varargin)
% gauss_legendre: Gauss-Legendre integration quadrature
% I = gauss_legendre(func,a,b,pts,varargin)
% inputs:
% func = name of function to be integrated
% a, b = integration limits
% pts = number of points for integration (default = 7)
% pl,p2,... = additional parameters used by func
% output:
% I = integral estimate

if nargin<3
    error('at least 3 input arguments required');
end
if nargin<4 || isempty(pts)
    pts = 7;
end
if abs(round(pts)) ~= pts 
    error('number of points must be an integer greater than zero');
end
if abs(round(pts)) == pts && pts > 7
    pts = 7;
end

a1 = (b+a)/2;
a2 = (b-a)/2;
if pts == 1
    c0 = 2; x0 = 0;
    I = a2*c0*func(a1+a2*x0,varargin{:});
elseif pts == 2
    c0 = 1; x0 = -1/sqrt(3);
    c1 = 1; x1 = 1/sqrt(3);
    I = a2*(c0*func(a1+a2*x0,varargin{:}) + c1*func(a1+a2*x1,varargin{:}));
elseif pts == 3
    c0 = 5/9; x0 = -sqrt(3/5);
    c1 = 8/9; x1 = 0;
    c2 = 5/9; x2 = sqrt(3/5);
    I = a2*(c0*func(a1+a2*x0,varargin{:}) + c1*func(a1+a2*x1,varargin{:}) + ...
            c2*func(a1+a2*x2,varargin{:}));
elseif pts == 4
    c0 = (18-sqrt(30))/36; x0 = -sqrt(525+70*sqrt(30))/35;
    c1 = (18+sqrt(30))/36; x1 = -sqrt(525-70*sqrt(30))/35;
    c2 = (18+sqrt(30))/36; x2 = sqrt(525-70*sqrt(30))/35;
    c3 = (18-sqrt(30))/36; x3 = sqrt(525+70*sqrt(30))/35;
    I = a2*(c0*func(a1+a2*x0,varargin{:}) + c1*func(a1+a2*x1,varargin{:}) + ...
            c2*func(a1+a2*x2,varargin{:}) + c3*func(a1+a2*x3,varargin{:}));
elseif pts == 5
    c0 = (322-13*sqrt(70))/900; x0 = -sqrt(245+14*sqrt(70))/21;
    c1 = (322+13*sqrt(70))/900; x1 = -sqrt(245-14*sqrt(70))/21;
    c2 = 128/225; x2 = 0;
    c3 = (322+13*sqrt(70))/900; x3 = sqrt(245-14*sqrt(70))/21;
    c4 = (322-13*sqrt(70))/900; x4 = sqrt(245+14*sqrt(70))/21;
    I = a2*(c0*func(a1+a2*x0,varargin{:}) + c1*func(a1+a2*x1,varargin{:}) + ...
            c2*func(a1+a2*x2,varargin{:}) + c3*func(a1+a2*x3,varargin{:}) + ...
            c4*func(a1+a2*x4,varargin{:}));
elseif pts == 6
    c0 = 0.171324492379170; x0 = -0.932469514203152; 
    c1 = 0.360761573048139; x1 = -0.661209386466265;
    c2 = 0.467913934572691; x2 = -0.238619186083197;
    c3 = 0.467913934572691; x3 =  0.238619186083197;
    c4 = 0.360761573048131; x4 =  0.661209386466265;
    c5 = 0.171324492379170; x5 =  0.932469514203152;
    I = a2*(c0*func(a1+a2*x0,varargin{:}) + c1*func(a1+a2*x1,varargin{:}) + ...
            c2*func(a1+a2*x2,varargin{:}) + c3*func(a1+a2*x3,varargin{:}) + ...
            c4*func(a1+a2*x4,varargin{:}) + c5*func(a1+a2*x5,varargin{:}));
elseif pts == 7
    c0 = 0.1294849661688697; x0 = -0.9491079123427585; 
    c1 = 0.2797053914892766; x1 = -0.7415311855993945;
    c2 = 0.3818300505051189; x2 = -0.4058451513773972;
    c3 = 0.4179591836734694; x3 = 0;
    c4 = 0.3818300505051189; x4 = 0.4058451513773972;
    c5 = 0.2797053914892766; x5 = 0.7415311855993945;
    c6 = 0.1294849661688697; x6 = 0.9491079123427585;
    I = a2*(c0*func(a1+a2*x0,varargin{:}) + c1*func(a1+a2*x1,varargin{:}) + ...
            c2*func(a1+a2*x2,varargin{:}) + c3*func(a1+a2*x3,varargin{:}) + ...
            c4*func(a1+a2*x4,varargin{:}) + c5*func(a1+a2*x5,varargin{:}) +...
            c6*func(a1+a2*x6,varargin{:}));
end

