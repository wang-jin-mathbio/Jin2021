%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the model for FUCCI spheroids using the Newton-Raphson method   %
% From paper manuscript "Mathematical model of tumour spheroid            %
%     experiments with real-time cell cycle imaging"                      %
% By Wang Jin and Matthew Simpson                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
figure
dxi = 0.001; % Set uniform spacing in the transformed system
dt = 0.001; % Set uniform time step in the transformed system
t = 0.0; % Initialse time t
tf = 2.4; % Final simulation time
trecord = [1.20 1.80 2.40]; % Time points for outputing numerical results
tol = 1e-6; % Set convergence tolerance for the Newton-Raphson method
l = 1.0; % Current spheroid radius          
pl = 1.0; % Spheroid radius at previous time step
BA = 1.5; % Parameter B
CA = 1.8; % Parameter C        
sigma = 1; % Parameter sigma
delta = 0.9; % Parameter delta
beta = 15; % Parameter beta
Cg = 0.1; % Parameter c_g      
Cr = 0.7; % Parameter c_r     
Cd = 0.5; % Parameter c_d
m1 = 10; % Parameter m_1          
m2 = 10; % Parameter m_2         
m3 = 10; % Parameter m_3
%--------------------------------------------------------------------------
maxsteps = round(tf/dt); % Total number of time stpes
N = round(1.0/dxi) + 1; % Total number of spatial nodes
NN1 = 1000; % Set initial cell population for r (parameter N_r)   
NN2 = 1000; % Set initial cell population for g (parameter N_g)
nn1 = NN1 / ( NN1+NN2 ) * ones(1,N); % Cell density r    
nn2 = NN2 / ( NN1+NN2 ) * ones(1,N); % Cell density g
pnn1 = nn1; % Cell density r at previous time step   
pnn2 = nn2; % Cell density g at previous time step
C = ones(1,N); % Oxygen density c       
V = zeros(1,N); % Cell velocity v
Lrecord = zeros(maxsteps+1,1); % Array for storing l at each time step
Lrecord(1,1) = l; % Initial condition for l
% Arrays for storing computing results directly from the Newton-Raphson method
delnn1 = zeros(1,N); % delta r    
delnn2 = zeros(1,N); % delta g
delc = zeros(1,N); % delta c    
delV = zeros(1,N); % delta v
%--------------------------------------------------------------------------
% Construct spatial coordinates for the transformed system
xx = zeros(1,N);     
xi = zeros(1,N);
for i=1:N
    xi(1,i)=0.0+(i-1)*dxi;
end
%--------------------------------------------------------------------------
% Initialise arrays for using the Newton-Raphson method
a = zeros(1,N);    b = zeros(1,N);    c = zeros(1,N);    d = zeros(1,N);
%--------------------------------------------------------------------------
for i = 1:maxsteps
    t = t + dt; % Forward time
    kk = 0; % Clear previous iteration counter number
    delnn1 = ones(1,N);    delnn2 = ones(1,N); % Clear previous computational results
    while norm(delnn1,Inf) > tol || norm(delnn2,Inf) > tol
        kk = kk + 1; % Update iteration counter
% Solve discretized equations for cell density r written in Tridiagonal form
        % Left boundary
        dldt = V(1,N);
        kg = Kg( C(1,1),Cg,m1 );
        kr = Kr( C(1,1),BA,Cr,m2 );
        kd = Kd( C(1,1),CA,sigma,Cd,m3 );
        dc = (1 - delta) * kd;    bc = kg - dc;
        a(1,1) =   0.0;
        b(1,1) = - 1.0 / dt ...
                 - kr - kd - bc * nn2(1,1) ...
                 + 2 * dc * nn1(1,1);
        c(1,1) =   0.0;
        d(1,1) = ( nn1(1,1)-pnn1(1,1) ) / dt - 2 * kg * nn2(1,1)...
                 + nn1(1,1) * ( kr + kd + bc * nn2(1,1) ) ...
                 - nn1(1,1)^2 * dc;
        % Internal nodes
        for j = 2:N-1
            kg = Kg( C(1,j),Cg,m1 );
            kr = Kr( C(1,j),BA,Cr,m2 );
            kd = Kd( C(1,j),CA,sigma,Cd,m3 );
            dc = (1 - delta) * kd;    bc = kg - dc;
            if ( V(1,j)/l - xi(1,j)*dldt/l ) >= 0 % Upwinding
                a(1,j) = ( V(1,j)/l - xi(1,j)*dldt/l ) / dxi;
                b(1,j) = - 1.0 / dt - ( V(1,j) - xi(1,j)*dldt ) / l / dxi ...
                         - bc * nn2(1,j) - kd - kr ...
                         + 2 * dc * nn1(1,j);
                c(1,j) =   0.0;
                d(1,j) = ( nn1(1,j)-pnn1(1,j) ) / dt - 2 * kg * nn2(1,j) ...
                         + ( V(1,j) - xi(1,j)*dldt ) * ( nn1(1,j) - nn1(1,j-1) ) / l / dxi ...
                         + nn1(1,j) * ( kr + kd + bc * nn2(1,j) ) ...
                         + nn1(1,j)^2 * dc;
            elseif ( V(1,j)/l - xi(1,j)*dldt/l ) < 0
                a(1,j) = 0.0;
                b(1,j) = - 1.0 / dt + ( V(1,j) - xi(1,j)*dldt ) / l / dxi ...
                         - bc * nn2(1,j) - kr - kd ...
                         + 2 * dc * nn1(1,j);
                c(1,j) = ( xi(1,j)*dldt/l - V(1,j)/l ) / dxi;
                d(1,j) = ( nn1(1,j)-pnn1(1,j) ) / dt - 2 * kg * nn2(1,j) ...
                         + ( V(1,j) - xi(1,j)*dldt ) * ( nn1(1,j+1) - nn1(1,j) ) / l / dxi ...
                         + nn1(1,j) * ( kr + kd + bc * nn2(1,j) ) ...
                         - nn1(1,j)^2 * dc;
            end        
        end
        % Right boundary
        kg = Kg( C(1,N),Cg,m1);
        kr = Kr( C(1,N),BA,Cr,m2 );
        kd = Kd( C(1,N),CA,sigma,Cd,m3 );
        dc = (1 - delta) * kd;    bc = kg - dc;
        a(1,N) = 0.0;
        b(1,N) = - 1.0 / dt ...
                 - bc * nn2(1,N) - kd - kr ...
                 + 2 * dc * nn1(1,N);
        c(1,N) = 0.0;
        d(1,N) = ( nn1(1,N)-pnn1(1,N) ) / dt - 2 * kg * nn2(1,N)...
                 + nn1(1,N) * ( kr + kd + bc * nn2(1,N) ) ...
                 - nn1(1,N)^2 * dc;
        delnn1 = thomas(N,a,b,c,d); % Solve discretized equations for r
        nn1(1,:) = nn1(1,:) + delnn1(1,:); % Update solution for r
% Solve discretized equations for cell density g written in Tridiagonal form        
        % Left boundary
        kg = Kg( C(1,1),Cg,m1 );
        kr = Kr( C(1,1),BA,Cr,m2 );
        kd = Kd( C(1,1),CA,sigma,Cd,m3 );
        dc = (1 - delta) * kd;    bc = kg - dc;
        a(1,1) =   0.0;
        b(1,1) = - 1.0 / dt ...
                 - kg - kd - 2 * bc * nn2(1,1) ...
                 + dc * nn1(1,1);
        c(1,1) =   0.0;
        d(1,1) = ( nn2(1,1)-pnn2(1,1) ) / dt - kr * nn1(1,1)...
                 + nn2(1,1) * ( kg + kd - dc * nn1(1,1) ) ...
                 + nn2(1,1)^2 * bc;
        % Internal nodes     
        for j = 2:N-1
            kg = Kg( C(1,j),Cg,m1 );
            kr = Kr( C(1,j),BA,Cr,m2 );
            kd = Kd( C(1,j),CA,sigma,Cd,m3 );
            dc = (1 - delta) * kd;    bc = kg - dc;
            if ( V(1,j)/l - xi(1,j)*dldt/l ) >= 0 % Upwinding
                a(1,j) = ( V(1,j)/l - xi(1,j)*dldt/l ) / dxi;
                b(1,j) = - 1.0 / dt - ( V(1,j) - xi(1,j)*dldt ) / l / dxi ...
                         + dc * nn1(1,j) - kd - kg ...
                         - 2 * bc * nn2(1,j);
                c(1,j) =   0.0;
                d(1,j) = ( nn2(1,j)-pnn2(1,j) ) / dt - kr * nn1(1,j) ...
                         + ( V(1,j) - xi(1,j)*dldt ) * ( nn2(1,j) - nn2(1,j-1) ) / l / dxi ...
                         + nn2(1,j) * ( kg + kd - dc * nn1(1,j) ) ...
                         + nn2(1,j)^2 * bc;
            elseif ( V(1,j)/l - xi(1,j)*dldt/l ) < 0
                a(1,j) = 0.0;
                b(1,j) = - 1.0 / dt + ( V(1,j) - xi(1,j)*dldt ) / l / dxi ...
                         + dc * nn1(1,j) - kg - kd ...
                         - 2 * bc * nn2(1,j);
                c(1,j) = ( xi(1,j)*dldt/l - V(1,j)/l ) / dxi;
                d(1,j) = ( nn2(1,j)-pnn2(1,j) ) / dt - kr * nn1(1,j) ...
                         + ( V(1,j) - xi(1,j)*dldt ) * ( nn2(1,j+1) - nn2(1,j) ) / l / dxi ...
                         + nn2(1,j) * ( kg + kd - dc * nn1(1,j) ) ...
                         + nn2(1,j)^2 * bc;
            end        
        end
        % Right boundary
        kg = Kg( C(1,N),Cg,m1);
        kr = Kr( C(1,N),BA,Cr,m2 );
        kd = Kd( C(1,N),CA,sigma,Cd,m3 );
        dc = (1 - delta) * kd;    bc = kg - dc;
        a(1,N) = 0.0;
        b(1,N) = - 1.0 / dt ...
                 + dc * nn1(1,N) - kd - kg ...
                 - 2 * bc * nn2(1,N);
        c(1,N) = 0.0;
        d(1,N) = ( nn2(1,N)-pnn2(1,N) ) / dt - kr * nn1(1,N)...
                 + nn2(1,N) * ( kg + kd - dc * nn1(1,N) ) ...
                 + nn2(1,N)^2 * bc;
        delnn2 = thomas(N,a,b,c,d); % Solve discretized equations for g
        nn2(1,:) = nn2(1,:) + delnn2(1,:); % Update solution for g
% Solve discretized equations for cell velocity v written in Tridiagonal form        
        % Left boundary
        a(1,1) = 0.0;
        b(1,1) = 1.0;
        c(1,1) = 0.0;
        d(1,1) = -1*V(1,1);
        % Internal nodes 
        for j = 2:N-1
            kg = Kg( C(1,j),Cg,m1 );
            kd = Kd( C(1,j),CA,sigma,Cd,m3 );
            dc = (1 - delta) * kd;    bc = kg - dc;
            a(1,j) =   1 / ( 2*dxi );
            b(1,j) = - 2 / xi(1,j);
            c(1,j) = - 1 / ( 2*dxi );
            d(1,j) = ( V(1,j+1)-V(1,j-1) ) / ( 2*dxi ) ...
                     + 2 * V(1,j) / xi(1,j) ...
                     - l * ( bc * nn2(1,j) - dc * nn1(1,j) );
        end
        % Right boundary
        kg = Kg( C(1,N),Cg,m1 );
        kd = Kd( C(1,N),CA,sigma,Cd,m3 );
        dc = (1 - delta) * kd;    bc = kg - dc;
        a(1,N) =  1 / dxi;
        b(1,N) = -2 - 1/dxi;
        c(1,N) =  0.0;
        d(1,N) = ( V(1,N)-V(1,N-1) ) / dxi ...
                 + 2 * V(1,N) ...
                 - l * ( bc * nn2(1,N) - dc * nn1(1,N) );

        delV = thomas(N,a,b,c,d); % Solve discretized equations for v
        V(1,:) = V(1,:) + delV(1,:); % Update solution for v
% Solve discretized equations for oxygen c written in Tridiagonal form
        % Left boundary
        a(1,1) =  0.0;
        b(1,1) = -1.0;
        c(1,1) =  1.0;
        d(1,1) = -( C(1,2) - C(1,1) );
        % Internal nodes 
        for j = 2:N-1 
            kg = Kg( C(1,j),Cg,m1 );
            a(1,j) = - ( xi(1,j-1) + xi(1,j) )^2 / 4 / xi(1,j)^2 / dxi^2;
            b(1,j) =   ( (xi(1,j-1) + xi(1,j))^2 + (xi(1,j) + xi(1,j+1))^2 ) / 4 / xi(1,j)^2 / dxi^2 ...
                     + l^2 * nn2(1,j) * ( beta*m1*Cg^m1*C(1,j)^(m1-1) / (C(1,j)^m1+Cg^m1)^2 );
            c(1,j) = - ( xi(1,j) + xi(1,j+1) )^2 / 4 / xi(1,j)^2 / dxi^2;
            d(1,j) =   ( ( xi(1,j)+xi(1,j+1) )^2 * (C(1,j+1)-C(1,j)) ...
                       - ( xi(1,j-1)+xi(1,j) )^2 * (C(1,j)-C(1,j-1)) ) / 4 / xi(1,j)^2 / dxi^2 ...
                     - l^2 * nn2(1,j) * beta * kg;
        end
        % Right boundary
        a(1,N) =  0.0;
        b(1,N) =  1.0;
        c(1,N) =  0.0;
        d(1,N) = -( C(1,N) - 1.0 );
        delc = thomas(N,a,b,c,d); % Solve discretized equations for c
        C(1,:) = C(1,:) + delc(1,:); % Update solution for c

        l = pl + dt * V(1,N); % UPdate spheroid radius
    end
    %----------------------------------------------------------------------
    % Output time and iteration counter every 100 time steps
    if mod(i,100)==0
        fprintf('Time %d\n',t); 
        fprintf('Iteration %d\n',kk);
    end
    %----------------------------------------------------------------------
    pl = l; % Update l at previous time steps                 
    Lrecord(i+1,1) = l; % Store l at current time step
    pnn1(1,:) = nn1(1,:); % Update r at previous time steps
    pnn2(1,:) = nn2(1,:); % Update g at previous time steps
    xx = l * xi; % Construct spatial coordinates in the oringinal system
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numerical results visualisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any( abs(trecord(:)-t) < dt/10 )
        subplot(2,3,1)
        hold on
        plot(xx,nn1,'r','LineWidth',2)
        plot(xx,nn2,'g','LineWidth',2)
        plot(xx,C,'b','LineWidth',2)
        ylim([0 1])
        xlim([0 l])
        xlabel('x')
        ylabel('Density')
        if abs(t-tf) < dt/10
            legend('r(x,t)','g(x,t)','c(x,t)')
        end
        pbaspect([2 1 1])
        
        subplot(2,3,2)
        hold on
        plot(xx,nn1+nn2,'k','LineWidth',2)
        plot(xx,nn1./(nn1+nn2),'r','LineWidth',2)
        ylim([0 1])
        xlim([0 l])
        xlabel('x')
        ylabel('Density')
        if abs(t-tf) < dt/10
            legend('n(x,t)','r(x,t)/n(x,t)')
        end
        pbaspect([2 1 1])
        
        subplot(2,3,3)
        hold on
        plot(xx,V,'k','LineWidth',2)
        xlim([0 l])
        xlabel('x')
        ylabel('v(x,t)')
        pbaspect([2 1 1])
    end
end

subplot(2,3,4)
plot(0:dt:tf,Lrecord,'k','LineWidth',2)
xlabel('t')
ylabel('l(t)')
pbaspect([2 1 1])

subplot(2,3,5)
hold on
c = 0:1e-3:1;
kg = Kg( c,Cg,m1 );
kr = Kr( c,BA,Cr,m2 );
kd = Kd( c,CA,sigma,Cd,m3 );
plot(c,kg,'g','LineWidth',2)
plot(c,kr,'r','LineWidth',2)
plot(c,kd,'k','LineWidth',2)
xlabel('c')
ylabel('Rates')
legend('Kg','Kr','Kd')
pbaspect([2 1 1])


% Function of phase transition -- S/G2/M -> G1
function x = Kg(c,Cg,m1)
    x = c.^m1 ./ ( c.^m1 + Cg^m1 );
end

% Function of phase transition -- G1 -> S/G2/M
function x = Kr(c,BA,Cr,m2)
    x = BA * c.^m2 ./ ( c.^m2 + Cr^m2 );
end

% Death function
function x = Kd(c,CA,sigma,Cd,m3)
x = CA * ( 1 - sigma * c.^m3 ./ ( c.^m3 + Cd^m3 ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Algorithm
function x = thomas(N,a,b,c,d)
x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end
    
    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end