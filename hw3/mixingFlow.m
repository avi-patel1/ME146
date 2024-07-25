% Avi Patel
% Mixing-Driven reactions Flow in 2D 
% 2/26/2024


clear clr;
close all;
import random_perm.*;

%% Define the grid 
Lx = 2;
Ly = 1;

Nx = 256; 
Ny = 128; 

dx = Lx/Nx;
dy = Ly/Ny;

x = dx/2:dx:Lx-dx/2;
y = dy/2:dy:Ly-dy/2;

[xx,yy] = meshgrid(x,y);

% initialize concentration field 
cA = zeros(Ny, Nx);

cB = ones(Ny, Nx);
cC = zeros(Ny, Nx);

Pe = 2000;
Da = 20;

ubar = 1;
par = 0; 

% initialize permeability field
var_lnk = 0.1;  % set to 0.1 or 2
corr_lenx = 4 * dx;
corr_leny = 4 * dx;
rng("default");

% uncomment to use random permeability field 
k_xy = ones(Ny, Nx); % random_perm(var_lnk, corr_lenx, corr_leny, Nx, Ny, Lx, Ly);

% initialize time loop
t = 0;
t_final = 2;

dt_diff = 0.1*Pe*dx*dx;
dt = 0.01*dt_diff;
count = 0;
cC_visc_lst = {};
dt2_lst = {};

tic;
while(t < t_final)
    
    R = 0; % set -3, 0, 3
    mu_c = exp(-R.*cA);
    cA_left = 0.99+0.01.*rand(Ny,1);
    cB_left = 0;
    cC_left = 0;

    lambda = k_xy ./ mu_c; 
    
    pmat = zeros(Nx,Ny);
    flipud(pmat);
    
    % transmissibility
    Tx = zeros(Ny, Nx+1);
    Tx(:,2:Nx) = (2 ./ (1./lambda(:,1:Nx-1) + 1./lambda(:,2:Nx)) ) *dy/dx; % interior walls

    
    Ty = zeros(Ny+1, Nx);
    Ty(2:Ny, :) = ( 2./ (1./lambda(1:Ny-1,:) + 1./lambda(2:Ny,:)) ) * dx/dy; % interior walls
        
    % prescribed pressure b.c. on right side 
    Tx(:,Nx+1) = 2*lambda(:,Ny)*dy/dx;
    
    % Create the A matrix
    Tyvecm1 = reshape(Ty(1:Ny,:), Nx*Ny, 1);
    Tyvecp1 = reshape(Ty(2:Ny+1,:),Nx*Ny,1);
    TxvecmNy = reshape(Tx(:,1:Nx),Nx*Ny,1);
    TxvecpNy = reshape(Tx(:,2:Nx+1),Nx*Ny,1);
    
    T0vec = Tyvecm1 + Tyvecp1 + TxvecmNy + TxvecpNy;
    
    Amat = spdiags([-TxvecmNy,-Tyvecm1, T0vec, -Tyvecp1, -TxvecpNy],[Ny, 1,0,-1,-Ny], Nx*Ny, Nx*Ny);
    

    % create RHS 
    bmat = zeros(Ny,Nx);
    bmat(:,end) = Tx(:,end)*par;
    bvec = reshape(bmat,Nx*Ny,1);
    
    % prescribed flow b.c. RHS    
    bvec(1:Ny) = ubar * dy;
    

    %% solve the pressure
    
    pvec = Amat \ bvec; 
    pmat = reshape(pvec, Ny, Nx);
    
    % initialize velocity matrices
    u = zeros(Ny, Nx+1); 
    v = zeros(Ny+1, Nx);
    
    u(:, 2:Nx+1) = Tx(:, 2:Nx+1) .* (pmat(:, 1:Nx) - [pmat(:, 2:Nx),0*ones(Ny,1)]) / dx; % Tx.*  ([pl,pmat] - [pmat,pr])/dx;
    v(2:Ny,:) = Ty(2:Ny,:) .* (pmat(1:Ny-1,:) - pmat(2:Ny,:))/dy;

    Umax = max(abs(u(:)));
    Vmax = max(abs(v(:)));
    Umax = max([Umax, Vmax]);
    
    %% solve for the transport equation
       
    % advective flux
    Fxadv_A = zeros(Ny, Nx+1);
    Fxadv_B = zeros(Ny, Nx+1);
    Fxadv_C = zeros(Ny, Nx+1);

    Fxadv_A(:,2:Nx) = (u(:,2:Nx).*cA(:,1:Nx-1).*(u(:,2:Nx) > 0) + u(:,2:Nx).*cA(:,2:Nx).*(u(:,2:Nx) < 0)).*dy;
    Fxadv_B(:,2:Nx) = (u(:,2:Nx).*cB(:,1:Nx-1).*(u(:,2:Nx) > 0) + u(:,2:Nx).*cB(:,2:Nx).*(u(:,2:Nx) < 0)).*dy;
    Fxadv_C(:,2:Nx) = (u(:,2:Nx).*cC(:,1:Nx-1).*(u(:,2:Nx) > 0) + u(:,2:Nx).*cC(:,2:Nx).*(u(:,2:Nx) < 0)).*dy;

    Fxadv_A(:,1) = ubar*dy*cA_left; 
    Fxadv_A(:,Nx+1) = u(:,Nx+1).*dy.*cA(:,Nx);
    Fxadv_B(:,Nx+1) = u(:,Nx+1).*dy.*cB(:,Nx);
    
    Fxadv_C(:,1) = ubar*dy*cC_left;
    Fxadv_C(:,Nx+1) = u(:,Nx+1).*dy.*cC(:,Nx);

    Fyadv_A = zeros(Ny+1, Nx);
    Fyadv_B = zeros(Ny+1, Nx);
    Fyadv_C = zeros(Ny+1, Nx);

    Fyadv_A(2:Ny,:) = (v(2:Ny,:).*cA(1:Ny-1,:).*(v(2:Ny,:) > 0) + v(2:Ny,:).*cA(2:Ny,:).*(v(2:Ny,:) < 0)).*dx;
    Fyadv_B(2:Ny,:) = (v(2:Ny,:).*cB(1:Ny-1,:).*(v(2:Ny,:) > 0) + v(2:Ny,:).*cB(2:Ny,:).*(v(2:Ny,:) < 0)).*dx;
    Fyadv_C(2:Ny,:) = (v(2:Ny,:).*cC(1:Ny-1,:).*(v(2:Ny,:) > 0) + v(2:Ny,:).*cC(2:Ny,:).*(v(2:Ny,:) < 0)).*dx;

    % diffusive flux
    Fxdiff_A = zeros(Ny,Nx+1);
    Fxdiff_B = zeros(Ny,Nx+1);
    Fxdiff_C = zeros(Ny,Nx+1);

    Fxdiff_A(:,2:Nx) = 1/Pe*(cA(:,1:Nx-1) - cA(:,2:Nx)) .* dy/dx;
    Fxdiff_B(:,2:Nx) = 1/Pe*(cB(:,1:Nx-1) - cB(:,2:Nx)) .* dy/dx;
    Fxdiff_C(:,2:Nx) = 1/Pe*(cC(:,1:Nx-1) - cC(:,2:Nx)) .* dy/dx;

    Fx_A = (Fxadv_A + Fxdiff_A);
    Fx_B = (Fxadv_B + Fxdiff_B);
    Fx_C = (Fxadv_C + Fxdiff_C);

    Fydiff_A = zeros(Ny+1,Nx);
    Fydiff_B = zeros(Ny+1,Nx);
    Fydiff_C = zeros(Ny+1,Nx);

    Fydiff_A(2:Ny,:) = 1/Pe*(cA(1:Ny-1,:) - cA(2:Ny,:)) .* dx/dy; 
    Fydiff_B(2:Ny,:) = 1/Pe*(cB(1:Ny-1,:) - cB(2:Ny,:)) .* dx/dy; 
    Fydiff_C(2:Ny,:) = 1/Pe*(cC(1:Ny-1,:) - cC(2:Ny,:)) .* dx/dy; 

    Fy_A = (Fyadv_A + Fydiff_A);
    Fy_B = (Fyadv_B + Fydiff_B);
    Fy_C = (Fyadv_C + Fydiff_C);
       

    dt_adv= 0.5*dx/Umax;
    dt = min(dt_diff, dt_adv);

    % update the concentration
    cA = cA + (dt/(dx*dy))*(Fx_A(:,1:Nx)-Fx_A(:,2:Nx+1)+Fy_A(1:Ny,:)-Fy_A(2:Ny+1,:)) - dt*(Da*cA.*cB);
    cB = cB + (dt/(dx*dy))*(Fx_B(:,1:Nx)-Fx_B(:,2:Nx+1)+Fy_B(1:Ny,:)-Fy_B(2:Ny+1,:)) - dt*(Da*cA.*cB);
    cC = cC + (dt/(dx*dy))*(Fx_C(:,1:Nx)-Fx_C(:,2:Nx+1)+Fy_C(1:Ny,:)-Fy_C(2:Ny+1,:)) + dt*(Da*cA.*cB);
    
    cC_visc_lst{end + 1} = sum(cC,'all');
    dt2_lst{end+1} = t;
    t = t + dt;

    % uncomment to view concentration field
    if  mod(count,10) == 0
        %figure(1);
        subplot(3,1,1);
        imagesc(cA); colorbar; axis equal;
        title("Concentration Field of A: " + " , t = " + t)

        subplot(3,1,2);
        imagesc(cB); colorbar; axis equal;
        title("Concentration Field of B: " + " , t = " + t)

        subplot(3,1,3);
        imagesc(cC); colorbar; axis equal;
        title("Concentration Field of C: " + " , t = " + t)
    
        drawnow;
    end
    count = count + 1;
end
toc;

% res = zeros(length(pmat(:,1)),1);
% for i = 1:size(pmat,1)
%     for j = 1:size(pmat,2)
%         if i == j
%             res(i) = pmat(i,j);
%         end
%     end
% end

% Pressure Distribution plotting
% figure(1)
% plot(res);
% title("N = " + Nx)
% xlabel("x","FontSize",15)
% ylabel("p(x)", "FontSize",15)

% Breakthrough Curve Plotting
% hold off;
% figure(1)
% plot(cell2mat(c_out));
% title("R = " + R +", Pe = " + Pe)
% xlabel("t","FontSize",15)
% ylabel("c(t)","FontSize",15)
%end