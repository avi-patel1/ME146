%% Avi Patel 
%% 1D Tracer Transport for Li/Na: Plots 1D tracer transport at different Pe values

%% Grid 

L = 1;
N = 100;
dz = 1/N;
z = dz/2:dz:1-dz/2;

% Peclet Number
Pe = {529, 2647, 13235};
colors = {'blue', 'red', 'green'};

%% time integration step:
dt = dz^2;

% time loop
for i = 1:length(Pe)
    disp(i)
    % prescribe initial concentration and the velocity
    c_Li = zeros(N,1);
    c_Na = zeros(N,1);
    
    % Boundary Conditions 
    % Atacama initial Li/Na concentration ratio
    cleft_Li =  1/50.6666667;  
    cleft_Na = 1;  

    % saturation parameter
    csat_Na = 263622.975/76000;

    t=0;
    tmax=20;
    count = 1;

    while t < tmax
    
        z_bar = linspace(0,L,N+1);
    
        if t < csat_Na
            Fadv_Na = (1 - z_bar').*[cleft_Na;c_Na]; % (1 - z_bar').*c_Na(1:N-1);
            Fadv_Na(1) = 1*cleft_Na; 
            Fadv_Na(N+1) = 0;
        
            Fdiff_Na = -1/Pe{i} * ([c_Na;c_Na(end)] -[cleft_Na;c_Na]) / dz;
            c_Na = c_Na - dt/dz*(Fadv_Na(2:N+1) - Fadv_Na(1:N) + Fdiff_Na(2:N+1) - Fdiff_Na(1:N));
        end
                    
        % advective flux
        Fadv_Li = (1 - z_bar').*[cleft_Li;c_Li];
        
        % diffusive flux
        Fdiff_Li = -1/Pe{i} * ([c_Li;c_Li(end)] -[cleft_Li;c_Li]) / dz;
    
        F_Li = Fadv_Li + Fdiff_Li;
        
        c_Li = c_Li - dt/dz*(F_Li(2:N+1) - F_Li(1:N)); % (Fadv_Li(2:N+1) - Fadv_Li(1:N) + Fdiff_Li(2:N+1) - Fdiff_Li(1:N));
       
        t = t + dt;
    end

    figure(1); hold on;
    plot((c_Li)./c_Na, z, color=colors{i},LineWidth=2)
    axis([0 0.62 0.95 1]);
    pbaspect([0.5 1 1])
    xlabel("$c_{\rm~Li}/c_{\rm~Na}$", 'Interpreter', 'latex',FontSize=15)
    ylabel("$\frac{z}{H}$", 'Interpreter', 'latex',FontSize=15)
    if i == 3
        legend('Pe=529', 'Pe=2647', 'Pe=13235', Location='southeast')
    end
        
end