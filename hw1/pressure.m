% Avi Patel 
% 1D pressure
% 1/29/2024

L =1 ; % domain length
N = 100; % discretize the domain
dx = L/N;
x = dx/2:dx:L-dx/2;
dt = 10 ^ -5;
theta = 1; % explicit if 0, implicit if 1

% physical parameters 
% Mobilities an Transmissibility
lam = ones(N,1); % exp(10*x);
p0 = zeros(N,1);

I = eye(N);
D = (dx/dt) * I;

T = (1/dx) * spdiags([-lam 2*lam -lam],-1:1, N,N );
A = D + theta * T;
b = (D - (1 - theta)*T) * p0;

A(1,1) = (dx/dt) + theta * (3/dx);
A(1,2) = theta * (-1 / dx);
A(N,N) = (dx/dt) + theta * (1/dx);
A(N,N-1) = theta * (-1/dx);

b(1) = (dx/dt)*p0(1) - (1-theta)*((1/dx)*p0(1) + (2/dx)*p0(1) - (2/dx)*p0(1)) + (2/dx) * 1;
b(N) = (dx/dt)*p0(N) - (1-theta) * ((-1/dx)*p0(N-1) + (-1/dx) * p0(N) );


% time loop
t = 0;
tmax = 5;
count = 0;
tic;

while t < tmax

    p = A \ b;
    b = (D - (1-theta)*T)*p;
    % modify b based on B.C.
    b(1) = (dx/dt)*p(1) - (1-theta) * ((1/dx)*p(1) + (2/dx)*p(1) - (2/dx)*p(2)) + (2/dx) * 1;
    b(N) = (dx/dt)*p(N) - (1-theta) * ((-1/dx)*p(N-1) + (-1/dx) * p(N) ); 

    figure(1);
    plot(x,p,'k');
    title("1D Pressure: t=" + round(t,4));
    t = t + dt;
end
toc;