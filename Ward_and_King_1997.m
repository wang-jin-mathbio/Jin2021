%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving Ward & King's model (1997) using the Newton-Raphson method      %
% From paper manuscript "Mathematical model of tumour spheroid            %
%     experiments with real-time cell cycle imaging"                      %
% By Wang Jin and Matthew Simpson                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
figure
dxi = 0.001; % Set uniform spacing in the transformed system
dt = 0.001; % Set uniform time step in the transformed system
t = 0.0; % Initialse time t
tf = 100.0; % Final simulation time
trecord = [25 50 75 100]; % Time points for outputing numerical results
tol = 1e-6; % Set convergence tolerance for the Newton-Raphson method
l = 1.0; % Current spheroid radius 
pl = 1.0; % Spheroid radius at previous time step
% Parameter values set to reproduce Fig 1-4 from Ward & King's paper (1997)
BA = 1.0; % Parameter B/A
sigma = 0.9; % Parameter sigma
delta = 0.5; % Parameter delta
beta = 0.005; % Parameter beta
CC = 0.1; % Parameter c_c
CD = 0.05; % Parameter c_d
m1 = 1; % Parameter m_1
m2 = 1; % Parameter m_2
%--------------------------------------------------------------------------
maxsteps = tf/dt; % Total number of time stpes
N = 1.0/dxi+1; % Total number of spatial nodes
Lrecord = zeros(maxsteps+1,1); % Array for storing l at each time step
Lrecord(1,1) = 1.0; % Initial condition for l
xi = zeros(1,N); % Initialise spatial coordinates for the transformed system
xx = zeros(1,N); % Initialise spatial coordinates for the original system
% Arrays for storing numerical results
nn = ones(1,N); % Living cell density n
pnn = ones(1,N); % Living cell density n at previous time steps
C = ones(1,N); % Oxygen density c
V = zeros(1,N); % Cell velocity v
%--------------------------------------------------------------------------
% Arrays for storing computing results directly from the Newton-Raphson method
deln = zeros(1,N); % delta n 
delO = zeros(1,N); % delta C
delV = zeros(1,N); % delta v
%--------------------------------------------------------------------------
% Construct spatial coordinates for the transformed system
for i = 1:N
    xi(1,i) = 0.0 + (i-1) * dxi; 
end
%--------------------------------------------------------------------------
% Initialise arrays for using the Newton-Raphson method
a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,N);
d = zeros(1,N);
%--------------------------------------------------------------------------
for i = 1:maxsteps
    t = t + dt; % Forward time
    kk = 0; % Clear previous iteration counter number
    delnn = ones(1,N); % Clear previous computational results
    while norm(delnn,Inf) > tol
        kk = kk + 1; % Update iteration counter
% Solve discretized equations for living cell density n written in Tridiagonal form
        % Left boundary
        dldt = V(1,N);
        km = mitosis(C(1,1),CC,m1);
        kd = death(C(1,1),BA,sigma,CD,m2);
        a(1,1) = 0.0;
        b(1,1) = -1.0/dt + (km-kd) - 2*nn(1,1)*(km-(1-delta)*kd);
        c(1,1) = 0.0;
        d(1,1) = (nn(1,1)-pnn(1,1))/dt - nn(1,1)*(km-kd) + nn(1,1)^2*(km-(1-delta)*kd);
        % Internal nodes
        for j = 2:N-1
            km = mitosis(C(1,j),CC,m1);
            kd = death(C(1,j),BA,sigma,CD,m2);
            if (xi(1,j)*dldt/l - V(1,j)/l) <= 0 
                a(1,j) = -(xi(1,j)*dldt/l - V(1,j)/l)/dxi;
                b(1,j) = -1.0/dt + (xi(1,j)*dldt/l-V(1,j)/l)/dxi + ...
                    (km-kd) - 2*nn(1,j)*(km-(1-delta)*kd);
                c(1,j) = 0.0;
                d(1,j) = (nn(1,j)-pnn(1,j))/dt - ...
                    (xi(1,j)*dldt/l-V(1,j)/l)*(nn(1,j) - ...
                        nn(1,j-1))/dxi - nn(1,j)*(km-kd) + ...
                            nn(1,j)^2*(km-(1-delta)*kd);
            end
            if (xi(1,j)*dldt/l - V(1,j)/l) > 0
                a(1,j) = 0.0;
                b(1,j) = -1.0/dt - (xi(1,j)*dldt/l-V(1,j)/l)/dxi + ...
                    (km-kd) - 2*nn(1,j)*(km-(1-delta)*kd);
                c(1,j) = (xi(1,j)*dldt/l-V(1,j)/l)/dxi;
                d(1,j) = (nn(1,j)-pnn(1,j))/dt - ...
                    (xi(1,j)*dldt/l-V(1,j)/l)*(nn(1,j+1)-nn(1,j))/dxi - ...
                        nn(1,j)*(km-kd) + nn(1,j)^2*(km-(1-delta)*kd);
            end        
        end
        % Right boundary
        km = mitosis(C(1,N),CC,m1);
        kd = death(C(1,N),BA,sigma,CD,m2);
        a(1,N) = 0.0;
        b(1,N) = 1.0;
        c(1,N) = 0.0;
        d(1,N) = -1.*( nn(1,N) - (km-kd)*exp(t*(km-kd))/...
            ( (km-kd) - (km-(1-delta)*kd)*(1.0-exp(t*(km-kd))) ) );
        delnn = thomas(N,a,b,c,d); % Solve discretized equations for n
        nn(1,:) = nn(1,:) + delnn(1,:); % Update solution for n

% Solve discretized equations for cell velocity v written in Tridiagonal form
        % Left boundary
        a(1,1) = 0.0;
        b(1,1) = 1.0;
        c(1,1) = 0.0;
        d(1,1) = -1*V(1,1);
        % Internal nodes
        for j = 2:N-1 
            km = mitosis(C(1,j),CC,m1);
            kd = death(C(1,j),BA,sigma,CD,m2);
            a(1,j) = -1 / (2.*dxi);
            b(1,j) = 2.0 / ( xi(1,j) );
            c(1,j) = 1 / (2.*dxi);
            d(1,j) = -(V(1,j+1)-V(1,j-1))/(2*dxi) - 2*V(1,j)/xi(1,j) + ...
                l*nn(1,j)*(km-(1-delta)*kd);
        end
        % Right boundary
        km = mitosis(C(1,N),CC,m1);
        kd = death(C(1,N),BA,sigma,CD,m2);
        a(1,N) = -1/dxi;
        b(1,N) = 2.0 + 1/dxi;
        c(1,N) = 0.0;
        d(1,N) = -(V(1,N)-V(1,N-1))/dxi - 2*V(1,N) + ...
            l*nn(1,N)*(km-(1-delta)*kd);
        delV = thomas(N,a,b,c,d); % Solve discretized equations for v
        V(1,:) = V(1,:) + delV(1,:); % Update solution for v
% Solve discretized equations for oxygen c written in Tridiagonal form
        % Left boundary
        a(1,1) = 0.0;
        b(1,1) = -1.0;
        c(1,1) = 1.0;
        d(1,1) = -1.0*(C(1,2)-C(1,1));
        % Internal nodes
        for j=2:N-1
            km = mitosis(C(1,j),CC,m1);
            kd = death(C(1,j),BA,sigma,CD,m2);
            a(1,j) = 1/(dxi^2)-1/(xi(1,j)*dxi);
            b(1,j) = -2.0/dxi^2-l*nn(1,j)*beta*m1*C(1,j)^(m1-1)*CC^m1/(C(1,j)^m1*CC^m1)^2;
            c(1,j) = 1/(dxi^2)+1/(xi(1,j)*dxi);
            d(1,j) = -1*(C(1,j+1)-2*C(1,j)+C(1,j-1))/(dxi^2)...
                -(C(1,j+1)-C(1,j-1))/(xi(1,j)*dxi)...
                    +l^2*nn(1,j)*beta*km;
        end
        % Right boundary 
        a(1,N) = 0.0;
        b(1,N) = 1.0;
        c(1,N) = 0.0;
        d(1,N) = -(C(1,N)-1.0);
        delO = thomas(N,a,b,c,d); % Solve discretized equations for c
        C(1,:) = C(1,:) + delO(1,:); % Update solution for c

        l = pl + dt*V(1,N); % UPdate spheroid radius
    end
    %----------------------------------------------------------------------
    % Output time and iteration counter every 1000 time steps
    if mod(i,1000) == 0
        fprintf('Time %d\n',t); 
        fprintf('Iteration %d\n',kk);
    end
    %----------------------------------------------------------------------
    pnn(1,:) = nn(1,:); % Update n at previous time steps
    pl = l; % Update l at previous time steps
    Lrecord(i+1,1) = l; % Store l at current time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numerical results visualisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any( abs(trecord(:)-t) < dt/10 )
        % Construct spatial coordinates in the oringinal system
        for i = 1:N
            xx(1,i) = l*xi(1,i);
        end
        %------------------------------------------------------------------
        subplot(2,2,2)
        hold on
        plot(xx,nn,'k','LineWidth',2)
        xlabel('Radius x')
        ylabel('Living cell density n')
        subplot(2,2,3)
        hold on
        plot(xx,C,'k','LineWidth',2)
        xlabel('Radius x')
        ylabel('Concentration c')
        subplot(2,2,4)
        hold on
        plot(xx,V,'k','LineWidth',2)
        xlabel('Radius x')
        ylabel('Velocity v')
    end
end
subplot(2,2,1)
plot(dt*(1:maxsteps+1),Lrecord,'k','LineWidth',2)
xlabel('Time t')
ylabel('Spheroid radius l(t)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Algorithm
function x = thomas(N,a,b,c,d)
    x = zeros(1,N);
    bb = b;
    dd = d;
    for i = 2:N
        ff = a(i) / bb(i-1);
        bb(i) = bb(i) - c(i-1)*ff;
        dd(i) = dd(i) - dd(i-1)*ff;
    end
    for i = 1:N-1
        x(N) = dd(N)/bb(N);    
        j = N - i;
        x(j) = ( dd(j)-c(j)*x(j+1) ) / bb(j);
    end
end

%Growth function
function x = mitosis(c,CC,m1)
    x = c^m1 / ( c^m1+CC^m1 );
end

%Death function
function x = death(c,BA,sigma,CD,m2)
    x = BA * ( 1 - sigma*c^m2/(c^m2+CD^m2) );
end

