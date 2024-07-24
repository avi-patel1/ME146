%% Grid 

L = 1;
N = 100;
dx = 1/N;
x = dx/2:dx:1-dx/2;

% prescribe initial concentration and the velocity
c= zeros(N,1);

u = ones(N+1,1);
%u =[sin(0);sin(2*pi*x)'];

% Boundary Conditions
cleft = 1;
cright = 1;

% Peclet Number
Pe = 100;

%% time integration step:

dt_adv = 0.5 * dx/max(abs(u));
dt_diff = 0.1*Pe*dx^2;
dt = min(dt_adv, dt_diff);

t=0;
tmax=3;

% time loop
while t < tmax
   
    % advective flux
    % upwind
    Fadv = u.*(u>0).*[cleft;c]+u.*(u<0).*[c;cright];
    
    % diffusive flux
    Fdiff = (-1/Pe)*([c;c(end)] -[cleft;c]) / dx;

    c = c - dt/dx*(Fadv(2:N+1) - Fadv(1:N) + Fdiff(2:N+1) - Fdiff(1:N));
    
    
    figure(1); axis([0,1,0,1])
    scatter(x,c, 'k');
    title("1D Transport: t = " + t);
    xlabel('x position');
    ylabel('concentration');
    t = t + dt;
end