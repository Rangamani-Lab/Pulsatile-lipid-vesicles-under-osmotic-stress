function swell_burst

clear all
format long
tic
randn('state',10) % initialize the random number generator


%% Definition of the Parameters %%

p.R0 = 14e-6 ; % radius of the vesicule in its unstreched state (m)
p.c0 = 0.2e3 ; % initial concentration of solute in the vesicle (mol/m3)
p.d= 3.5e-9;% membrane thickness (m)
p.rho_s = 1000 ; % density of the solvent (kg/m3)
p.nu_s = 18.04e-6; % solvent molar volume (m3/mol)
p.P = 2e-5; %m/s
p.kB=1.38e-23; % Boltzmann constant (J/K)
p.T=273+21 ;  %temperature in K
p.gamma =5e-12; % line tension at the pore (J/m)
p.eta_l= 5; % lipid monolayer viscosity (Pa s)
p.eta_s=1e-3; % solvent (water) viscosity (Pa s)
p.NA=6.023e23;% Avogadro number (mol-1)
p.D = 5e-10 ; % diffusivity of sucrose in water (m2/s) [Linder 1976]
p.C=2*pi; % viscous friction coefficient
p.kappa_eff=2e-3 ; % effectif stretching modulus (N/m)

p.A0 = 4*pi*p.R0^2 ; % reference area 

% Parameters related to fluctuations in stochastic model
p.Tfluct=1/150; % Period of the thermal fluctuation (s)

% Parameters related to critical tension in deterministic model
p.sc=p.kappa_eff*0.15; % critical tension
p.rp = p.gamma/p.sc; % intial pore radius at nucleation


%% Computation %%

X0=[p.R0 ; 0 ; p.c0]; t0=0; % Initial conditions
tf=t0+500; % Final time
dt=1e-4; % time step (in stochastic must be <= to p.Tfluct)

% Either stochastic
%[T, X] = swell_burst_rk1(@(X)swell_burst_system(X,p), t0, tf, X0, dt, p);

% Or deterministic
[T, X] = swell_burst_det(@(X)swell_burst_system(X,p), t0, tf, X0, dt, p); 

R=X(:,1); % results for GUV radius
r=X(:,2); % results for pore radius
c=X(:,3); % results for solute concentration

% save('pusatile_guv_00') % saves the results in .mat file 

toc


%% Visualisation %%

figure('Position',[0 1000 560 420*3/2])
subplot(3,1,1) % plots GUV radius as function of time
    plot(T,R*1e6,'LineWidth',2);
    hold on
    ylabel('R (\mum)')
    xlabel('t (s)')
    set(gca,'FontName', 'Times','FontSize', 16)
subplot(3,1,2) % plots pore radius as function of time
    plot(T,r*1e6,'LineWidth',2)
    ylabel('r (\mum)')
    xlabel('t (s)')
    set(gca,'FontName', 'Times','FontSize', 16)
subplot(3,1,3) % plots solute concentration as function of time
    plot(T,c,'LineWidth',2)
    ylabel('\Deltac (mM)')
    xlabel('t (s)')
    set(gca,'FontName', 'Times','FontSize', 16)

    
    
    
function [dX,dXf] = swell_burst_system(X,p)

R=1; r=2; c=3; % indices for the vector variable X

A =  4*pi*X(R)^2-pi*X(r)^2 ; % membrane area
sigma = p.kappa_eff*(A-p.A0)/p.A0 ; % effective tension

% ODEs 
dR =  1/(4*pi*X(R)^2) * (p.P*p.nu_s/(p.kB*p.T*p.NA) * ...
        ( p.kB*p.T*p.NA*X(c) - 2*sigma/X(R) ) * A ...
        - 2*sigma*X(r)^3 /(3*X(R)*p.eta_s)    );
dr =  1/(2*p.d*p.eta_l + p.C*p.eta_s*X(r))  *  (sigma*X(r) - p.gamma) ;   
dc =  1/(4/3*pi*X(R)^3) * ( pi*X(r)^2* (...
         - 2*sigma/X(R) * X(r)/(3*pi*p.eta_s) * X(c) ...
         - p.D*X(c)/X(R) ) - 4*pi*X(R)^2*X(c)*dR    );

dX = [dR ; dr ; dc];  % ODE system (drift)
dXf = [0 ; 1/(2*pi)*sqrt( 2*p.kB*p.T / (2*p.d*p.eta_l + p.C*p.eta_s*X(r)) ) ; 0]; % fluctuations (diffusion)



function [tout,Xout] = swell_burst_rk1(system, t0, tf, X0 , dt, p)
%Runge-Kutta scheme for the swell-burst cycles

disp(['Stochastic model started'])

Nsteps=round(tf/dt) ; 
Nout= min(Nsteps, 1e6); %maximum of output is 1e6 time points
%Nout= Nsteps;
X=X0;
Xout=zeros(Nout,numel(X0));  % Output variables matrix
Xout(1,:) = X0;
tout=[t0 ; zeros(Nout-1,1)]; % Output time vector
dW =sqrt(dt)*randn(1,Nsteps); % Wiener process
dtf=p.Tfluct; % period of the fluctuations
Nfluct = (tf-t0)/dtf; %% number of steps for fluctuations
N=round(dtf/dt);  %dtf = N*dt


av=0; % progress counter
j=2; % indice for output
for i=1: Nsteps-1    
    Xold=X;        
    [dX, dXf] = feval(system,Xold);    
    [dX2, dXf2] = feval(system,Xold+dXf*sqrt(dt));
    
    X=Xold + dX*dt;
    Xf=[0; 0; 0];  % Initiate drift

    if mod(i,round(Nsteps/Nfluct))==0  % Add diffusion at the fluctuation period
        Xf = sum(dW(i-N+1:i))*dXf + 0.5 * (dXf2 - dXf) .* (sum(dW(i-N+1:i)).^2 - dt) / sqrt(dt) ;
    end
    
    X=abs( max(0,X) + Xf );  % Insures non-negative variables 
    
    if mod(i,round(Nsteps/Nout))==0 % Saves only Nout data
        tout(j) = i*dt;
        Xout(j,:)  = X';
        j=j+1;
    end    
    
    if mod(i-1,round(Nsteps/10))==0 % Displays progress
        disp(['progress: ' num2str(av) '%'])
        av=av+10;
    end
end
disp(['progress: ' num2str(av) '%'])

  


function [tout,Xout] = swell_burst_det(system, t0, tf, X0, dt, p)
% Euler scheme for swell-burst cycles

disp(['Deterministic model started'])

Nsteps=round(tf/dt) ; 
Nout= min(Nsteps, 1e6); % maximum output is 1e6 timepoints

X=X0;
Xout=zeros(Nout,numel(X0));  % Output variables matrix
Xout(1,:) = X0;
tout=[t0 ; zeros(Nout-1,1)]; % Output time vector

av=0; % progress counter
j=2; % indice for output
for i=1: Nsteps -1    
    [dX, dXf] = feval(system,X);
    X = X + dX*dt;   
    X=max(0,X);  % Insures a non-negative variables
    
    A =  4*pi*X(1)^2-pi*X(2)^2 ; % membrane area
    sigma = p.kappa_eff*(A-p.A0)/p.A0 ; % effective tension
    
    if X(2)<=1e-9 && sigma >= p.sc  % Criteria for pore opening
       X(2)=p.rp; % Inital pore nucleation
    end  
    
    if mod(i,round(Nsteps/Nout))==0 % Saves only Nout data
        tout(j) = i*dt;
        Xout(j,:)  = X';
        j=j+1;
    end    
    if mod(i,round(Nsteps/10))==0 % Displays progress
        disp(['progress: ' num2str(av) '%'])
        av=av+10;
    end
end
disp(['progress: ' num2str(av) '%'])


    