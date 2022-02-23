%% script to draw a very useless but very beautiful rose with matlab.
% https://fr.mathworks.com/matlabcentral/communitycontests/contests/4/entries/701?source=15572&s_eid=psm_15572

n = 800;
p = pi;
[R, T] = ndgrid(linspace(0,1,n),linspace(-2,20*p,n));
x = 1-((0.5)*((5/4)*(1-mod(3.6*T,2*p)/p).^2-0.25).^2);
U = 2*exp(-T/(8*p));
L = sin(U);
J = cos(U);
y = 1.99*(R.^2).*(1.2*R-1).^2.*L;
K = x.*(R.*L+y.*J);
X = K.*sin(T);
Y = K.*cos(T);
Z = x.*(R.*J-y.*L);
surf(X,Y,Z,'LineStyle','none');
grid, axis off;
colormap(jet);