% Avi Patel
% coupling transport with flow in 2D
% 2/21/2024


clear clr;
close all;
import random_perm.*;

%% Define the grid 
Lx = 1;
Ly = 1;

Nx = 100; 
Ny = Nx; 

dx = Lx/Nx;
dy = Ly/Ny;

x = dx/2:dx:Lx-dx/2;
y = dy/2:dy:Ly-dy/2;

[xx,yy] = meshgrid(x,y);

% initialize concentration field 
c = zeros(Ny, Nx);
Pe = 100;

% initialize permeability field
var_lnk = 2;  % set to 0.1 or 2
corr_lenx = 4 * dx;
corr_leny = 4 * dx;
rng("default");
k_xy = random_perm(var_lnk, corr_lenx, corr_leny, Nx, Ny, Lx, Ly);

% initialize time loop
t = 0;
t_final = 2;
count = 0;
c_out = {};
tic;
while(t < t_final)
    
    R = 0; % set -3, 0, 3
    mu_c = exp(-R.*c);

    lambda = k_xy ./ mu_c; 
    
    pmat = zeros(Nx,Ny);
    flipud(pmat);
    
    % transmissibility
    Tx = zeros(Ny, Nx+1);
    Tx(:,2:Nx) = (2 ./ (1./lambda(:,1:Nx-1) + 1./lambda(:,2:Nx)) ) *dy/dx; % interior walls
    
    % No flow b.c.
    Tx(:,1) = zeros(1,Ny);  
    Tx(:,Nx+1) = zeros(1,Ny);
    
    Ty = zeros(Ny+1, Nx);
    Ty(2:Ny, :) = ( 2./(1./lambda(1:Ny-1,:) + 1./lambda(2:Ny,:)) ) * dx/dy;
        
    % prescribed flow b.c.
    Tx(1,1) = 2*lambda(1,1) * dy/dx;
    Ty(1,1) = 2*lambda(1,1) * dx/dy;
    
    % prescribed pressure b.c. 
    Tx(Ny,Nx+1) = 2*lambda(Nx,Ny)*dy/dx;
    Ty(Ny+1,Nx) = 2*lambda(Nx,Ny)*dx/dy;
    
    % Create the A matrix
    Tyvecm1 = reshape(Ty(1:Ny,:), Nx*Ny, 1);
    Tyvecp1 = reshape(Ty(2:Ny+1,:),Nx*Ny,1);
    TxvecmNy = reshape(Tx(:,1:Nx),Nx*Ny,1);
    TxvecpNy = reshape(Tx(:,2:Nx+1),Nx*Ny,1);
    
    T0vec = Tyvecm1 + Tyvecp1 + TxvecmNy + TxvecpNy;
    
    Amat = spdiags([-TxvecmNy,-Tyvecm1, T0vec, -Tyvecp1, -TxvecpNy],[Ny, 1,0,-1,-Ny], Nx*Ny, Nx*Ny);
    
    Amat(1,1) = Tx(1,2) + Ty(2,1);

    % create RHS 
    bmat = zeros(Ny,Nx);
    bvec = reshape(bmat,Nx*Ny,1);
    
    % prescribed flow b.c. RHS
    bvec(1) =  1;
    

    %% solve the pressure equation
    
    pvec = Amat \ bvec; 
    pmat = reshape(pvec, Ny, Nx);
        
    % initialize velocity matrices
    u = zeros(Ny, Nx+1);
    v = zeros(Ny+1, Nx);
    
    u(:, 2:Nx) = Tx(:, 2:Nx) .* (pmat(:, 1:Nx-1) - pmat(:, 2:Nx)) / dx; % Tx.*  ([pl,pmat] - [pmat,pr])/dx;
    v(2:Ny,:) = Ty(2:Ny,:) .* (pmat(1:Ny-1,:) - pmat(2:Ny,:))/dy;

    u(1,1) = 0.5/dy;
    u(Ny,Nx+1) = 0.5/dy;
    v(1,1) = 0.5/dx;
    v(Ny+1, Nx) = 0.5/dx;

    Umax = max(u(:));
    Vmax = max(v(:));
    Umax = max([Umax, Vmax]);
    
    %% solve for the transport equation
       
    % advective flux
    Fxadv = zeros(Ny, Nx+1);
    Fxadv(:,2:Nx) = (u(:,2:Nx).*c(:,1:Nx-1).*(u(:,2:Nx) > 0) + u(:,2:Nx).*c(:,2:Nx).*(u(:,2:Nx) < 0)).*dy;

    Fxadv(1,1) = 0.5;
    Fxadv(Ny,Nx+1) = u(Ny,Nx+1)*dy*c(Ny,Nx);

    Fyadv = zeros(Ny+1, Nx);
    Fyadv(2:Ny,:) = (v(2:Ny,:).*c(1:Ny-1,:).*(v(2:Ny,:) > 0) + v(2:Ny,:).*c(2:Ny,:).*(v(2:Ny,:) < 0)).*dx; % Fxadv(:,2:Nx)

    Fyadv(1,1) = 0.5;
    Fyadv(Ny+1,Nx) = v(Ny+1,Nx)*dx*c(Ny,Nx);

    % diffusive flux
    Fxdiff = zeros(Ny,Nx+1);
    Fxdiff(:,2:Nx) = 1/Pe*(c(:,1:Nx-1) - c(:,2:Nx)) .* dy/dx;
    Fx = (Fxadv + Fxdiff);

    Fydiff = zeros(Ny+1,Nx);
    Fydiff(2:Ny,:) = 1/Pe*(c(1:Ny-1,:) - c(2:Ny,:)) .* dx/dy; 
    Fy = (Fyadv + Fydiff);
       
    % Dynamic time step
    dt_diff=0.15*Pe*dx*dx;
    dt_adv= 0.5*dx/Umax;
    dt=min(dt_adv,dt_diff);

    % update the concentration
    c = c + (dt/(dx*dy))*(Fx(:,1:Nx)-Fx(:,2:Nx+1)+Fy(1:Ny,:)-Fy(2:Ny+1,:));
    t = t + dt;

    % uncomment to view concentration field
    if  mod(count,50) == 0
        figure(1);
        imagesc(c); colorbar; axis equal;
        title("Concentration Field: R= " + R + " , t = " + t)
        drawnow;
    end
    count = count + 1;

end
toc;

% Uncomment to view pressure/concentration field along diagonal axis 
% res = zeros(length(pmat(:,1)),1);
% for i = 1:size(pmat,1)
%     for j = 1:size(pmat,2)
%         if i == j
%             res(i) = pmat(i,j);
%         end
%     end
% end

% res = zeros(length(c(:,1)),1);
% for i = 1:size(c,1)
%     for j = 1:size(c,2)
%         if i == j
%             res(i) = c(i,j);
%         end
%     end
% end

% % Concentration Distribution plotting
% figure(2)
% plot(res);
% title("N = " + Nx)
% xlabel("x","FontSize",15)
% ylabel("p(x)", "FontSize",15)

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