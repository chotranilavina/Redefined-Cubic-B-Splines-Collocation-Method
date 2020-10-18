function [T,dh,dt] = CDR_CrankNicolson_Exp_InletMixing
%% 1D Convection-Diffusion-Reaction (CDR) equation 
% solves a d2u/dx2 + b du/dx +cu = du/dt 
tStart = tic; % Measuring computational time 
%% INPUT 
t_final=5 %: overall simulation time [s] 
n=10 %: Node number [-] 
h_tank=1; 
d_tank=0.58; 
vFR=11.75; 
v=((vFR*0.001)/60)/((pi*d_tank^2)/4); 
dh=h_tank/n; % [m] mesh size 
dt=1; % [s] time step 
maxk=round(t_final/dt,0); % Number of time steps 
% Constant CDR coefficients | a d2u/dx2 + b du/dx= du/dt 
a_const=0.02; 
b=-0.1;
c=0; %Coefficient for heat loss calculation 
alpha=1.17712434446770;
beta=-0.09;
eps_in=2;  
A_hyp=(eps_in-1)/(1-1/n); 
B_hyp=eps_in-A_hyp; 
Nsl=1:1:n; 
eps_eff=A_hyp./Nsl+B_hyp; 
a=a_const.*eps_eff; % Case 0 
% a=ones(n,1)*a_const; % Case 1 
% a=ones(n,1)*a_const*20; % Case 2 
T=zeros(n,maxk); 
T(:,1)=20; 
T(1,:)=60; 
%% Formation of Tridiagonal Matrices 
% Tridiagonal Matrix @Left-hand side 
subL(1:n-1)=-(2*dt*a(1:n-1)-dt*dh*b); % Sub diagonal - Coefficient of u_i-1,n+1 
maiL(1:n-0)=4*dh^2+4*dt*a-2*dh^2*dt*c; % Main diagonal - Coefficient of u_i,n+1 
supL(1:n-1)=-(2*dt*a(1:n-1)+dt*dh*b); % Super diagonal - Coefficient of u_i+1,n+1 
tdmL=diag(maiL,0)+diag(subL,-1)+diag(supL,1); 
% Tridiagonal Matrix @Right-hand side 
subR(1:n-1)=2*dt*a(1:n-1)-dt*dh*b; % Sub diagonal - Coefficient of u_i-1,n 
maiR(1:n-0)=4*dh^2-4*dt*a-2*dh^2*dt*c; % Main diagonal - Coefficient of u_i,n 
supR(1:n-1)=2*dt*a(1:n-1)+dt*dh*b; % Super diagonal - Coefficient of u_i+1,n 
tdmR=diag(maiR,0)+diag(subR,-1)+diag(supR,1); 

%% Boundary Condition - Matrices tdmL(1,1)=1; 
tdmL(1,2)=0; 
tdmL(end,end-1)=0; 
tdmL(end,end)=1; 
tdmR(1,1)=1; 
tdmR(1,2)=0; 
tdmR(end,end-1)=1; 
tdmR(end,end)=0; 
MMtr=tdmL\tdmR; 
%% Solution - System of Equations 
for j=2:maxk % Time Loop 
    Tpre=T(:,j-1); 
T(:,j)=MMtr*Tpre; 
if T(end,j)>=89.9 
        T(:,j+1:end)=[]; 
        finishedat=j*dt; 
        ChargingTime=sprintf('Charging time is %f [s]', finishedat) 
        tElapsed = toc(tStart); 
        SimulationTime=sprintf('Simulation time is %f [s]',tElapsed) 
return 
end
end
end