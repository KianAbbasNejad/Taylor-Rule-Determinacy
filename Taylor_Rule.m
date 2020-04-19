%% Taylor Rule Determinacy in Policy Space
% This code calculates the (in)determinacy of Taylor 
% rule monetary policy in log-linearised NK model 
% (Written by Kian Abbas Nejad)
% 
% Consdier a Discrete-time Calvo-lottery NK model with CRRA utility as a function 
% of consumption, leasure/labour and real money holdings (as described in Walsh 
% (2017) among many others). The log-linearised model here is summarised by 3 
% deterministic difference equations in output gap, inflation, and nominal 
% interst rate, $(y_t,\pi_t,i_t)$. 
%% Model Equations 
% The NKPC: 
% $$\pi_t=\beta E_t\pi_{t+1}+\kappa  y_t$$
%
% Where
% $$\kappa :=\frac{(1-\omega)(1-\beta \omega)(\sigma+\eta)}{\omega}$$
% 
% Where :
% $1-\omega=$Poisson arrival rate for intermediate firm price adjustment
% $\eta =$Coefficient of relative risk aversion w.r.t labour
% $\beta = $Household discount rate
% $\sigma =$Coefficient of relative risk aversion w.r.t consumption
% 
% The IS relation is: 
% $$y_t=E_t y_{t+1}-\frac{1}{\sigma}(i_t-E_t \pi_{t+1})$$
% 
% And Consider a forward-looking Taylor rule: 
% $$i_t=\chi_y E_t y_{t+1}+\chi_\pi E_t \pi_{t+1}$$
% 
% Substituting the Taylor rule in IS relation results in 2D system in  
% output gap and inflation with coefficient matrix A
% 
%% Parameter Values
% Here I use the parameters reported by Branch & McGough (2009), who use the 
% parameters estimated by Woodford (1999). 

clc
clearvars
close all
kappa = 0.024;
sigma = 0.157;
beta = 0.99;

%% Policy Space
% And take the monetary policy responses to inflation and output gap over the 
% reasonable range $[1,2]$: 

chi_pi = linspace(0,2,500);
chi_y  = linspace(0,2,500);
%% Determinacy
% Given that the system is in two jump variables, the Blanchard-Kahn condition 
% for determinacy states that $A$ must have two unstable eigenvalues 
% to prevent self-fullfilling equilibria. 

[PI Y]=meshgrid(chi_pi,chi_y);
det=ones(numel(chi_y),numel(chi_y));
cmap = colormap('gray'); %defining the colour map
cmap(1:20,:)=[]; %removing too dark colours
for i=1:numel(Y)
    A=inv([sigma-Y(i) 1-PI(i); 0 beta])*[sigma 0; -kappa 1];
    x = abs(eig(A)); 
    if min(x)>1
        det(i)=2; % value 2 if determinate
    elseif max(x)<1
        det(i)=0; %value 0 if order 2 indeterminate
    end  
end

fig1= figure(1);
[M,c]=contourf(chi_pi,chi_y,det);
fig1.Name = 'Determinacy';
title('Determinacy of Monetary Policy','Interpreter',"latex");
set(gca,'FontSize',14);
ylabel('$\chi_y$','Interpreter',"latex");
xlabel('$\chi_\pi$','Interpreter',"latex");
colormap(cmap);
hold on
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'o','Color',cmap(end,:));
h(2) = plot(NaN,NaN,'o','Color',cmap(end/2,:));
h(3) = plot(NaN,NaN,'o','Color',cmap(1,:));
legend(h, 'Determinate','Oder 1 Indet.','Order 2 Indet.','interpreter','latex');







