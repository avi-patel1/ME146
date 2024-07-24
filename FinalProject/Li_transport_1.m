%% Avi Patel 
%% 1D Tracer Transport for Li/Na: Plots 1D tracer transport at different time steps

%% Grid 

L = 1;
N = 100;
dz = 1/N;
z = dz/2:dz:1-dz/2;

% prescribe initial concentration and the velocity
c_Li = zeros(N,1);
c_Na = zeros(N,1);

% scaled initial Li concentrations from two different sources
cleft_Li_brine = 76000;
cleft_Li_sea =  10800;
csat = 226362.975;

% Boundary Conditions
% Atacama initial Li/Na concentration ratio
cleft_Li = 1/50.6666667; 
cleft_Na = 1;       

csat_brine = csat/cleft_Li_brine;
csat_sea = csat/cleft_Li_sea;

% Peclet Number
Pe_Li = 529; % 529, 2647, 13235
Pe_Na = 529; 

%% time integration step:
dt = dz^2;

t=0;
tmax=20;
count = 1;

timeVec =  [0,2,5,10,20]; % linspace(0,tmax, 500);

% time loop
while t < tmax

    z_bar = linspace(0,L,N+1);

    if t < csat_brine
        Fadv_Na = (1 - z_bar').*[cleft_Na;c_Na];
        Fadv_Na(1) = 1*cleft_Na; 
        Fadv_Na(N+1) = 0;
    
        Fdiff_Na = -1/Pe_Na * ([c_Na;c_Na(end)] -[cleft_Na;c_Na]) / dz;
        c_Na = c_Na - dt/dz*(Fadv_Na(2:N+1) - Fadv_Na(1:N) + Fdiff_Na(2:N+1) - Fdiff_Na(1:N));

    end
        
    % advective flux
    Fadv_Li = (1 - z_bar').*[cleft_Li;c_Li];
   
    % diffusive flux
    Fdiff_Li = -1/Pe_Li * ([c_Li;c_Li(end)] -[cleft_Li;c_Li]) / dz;

    F_Li = Fadv_Li + Fdiff_Li;
    
    c_Li = c_Li - dt/dz* (F_Li(2:N+1) - F_Li(1:N)); % (Fadv_Li(2:N+1) - Fadv_Li(1:N) + Fdiff_Li(2:N+1) - Fdiff_Li(1:N));
    
    t = t + dt;



    if abs(t - 2) < 10^-5
        subplot(1,2,1);
        plot((c_Li)./c_Na,z, color='blue',LineWidth=2);
        subplot(1,2,2);
        plot((c_Li)/cleft_Li,z, color='blue',LineWidth=2);
    end

    if abs(t - 5) < 10^-5
        subplot(1,2,1); hold on
        plot((c_Li)./c_Na,z,color='red',LineWidth=2);
        subplot(1,2,2); hold on
        plot((c_Li)/cleft_Li, z, color='red', LineWidth=2);
    end

    if abs(t - 10) < 10^-5
        subplot(1,2,1); hold on
        str = '#FFA500';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        plot((c_Li)./c_Na,z, color=color,LineWidth=2);
        subplot(1,2,2); hold on
        plot((c_Li)/cleft_Li, z, color=color, LineWidth=2);
    end

    if abs(t - 20) < 10^-5
        subplot(1,2,1); hold on
        str = '#800080';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        plot((c_Li)./c_Na,z, color=color, LineWidth=2);

        axis([0 1.6 0 1]); % axis([0 1.5 0 1]);
        pbaspect([0.5 1 1])
        legend('t=2', 't=5', 't=10', 't=20', Location='southeast')
        xlabel("$c_{\rm~Li}/c_{\rm~Na}$", 'Interpreter', 'latex',FontSize=15)
        ylabel("$\frac{z}{H}$", 'Interpreter', 'latex',FontSize=15)

        subplot(1,2,2); hold on
        plot((c_Li)/cleft_Li, z, color=color, LineWidth=2);
        pbaspect([0.5 1 1])
        xlabel("$c_{\rm~Li}/c^\circ_{\rm~Li}$", 'Interpreter', 'latex',FontSize=15)
        ylabel("$\frac{z}{H}$", 'Interpreter', 'latex',FontSize=15)
        
    end

end
