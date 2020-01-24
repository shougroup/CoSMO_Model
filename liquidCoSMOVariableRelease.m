%% this plots CoSMO dynamics of Model iii, and compares growth rates from simulations and calculations.
close all;
clear all;
global gL gA dL dA cA cL rA KL KA nL nA

LordRam
%L-A+ parameters
gL =0.51 ; %% max birth of L-A+ (/hr)
KL = 2.1;    %% Monod constant of L-A+ (uM)
nL=3.2; %% birth cooperativity
dL = 0.0024;  %% death rate of L-A+ (/hr)
cL = 5.4;    %% lysine consumption per new cell (fmole/cell) 
rA = 0.27;   %% fmole/cell/hr

% A-L+ parameters
gA =0.44 ; %% max birth of A-L+ (/hr)
KA = 1.3;   %% Monod constant of A-L+ (uM)
nA=1.5; %% birth cooperativity
dA =0.015 ;  %%  death rate of A-L+ (/hr)
cA = 3.1;    %% hyp consumption per new cell (fmole/cell)

%variable rL is in ODE function.

%%%Sam's experiment 20150922

shour=[0.00	5.83	22.67	29.50	46.50	53.50 74.50]; % time (hr)
sred=[1.75E+06	2.83E+06	4.95E+06	1.60E+07	3.03E+08	4.38E+08	1.69E+09
    1.75E+06	2.85E+06	4.88E+06	2.13E+07	2.85E+08	3.69E+08	1.19E+09
1.75E+06	2.80E+06	4.72E+06	2.06E+07	2.90E+08	3.76E+08 1.55E+09]; %Live L-A+ density (/ml), one experiment/row

sgreen=[1.75E+06	1.81E+06	9.88E+06	8.89E+06	2.47E+07	2.78E+08	1.94E+09
1.75E+06	1.80E+06	5.63E+06	5.95E+06	4.80E+07	3.54E+08	1.27E+09
1.75E+06	1.86E+06	5.07E+06	5.69E+06	4.84E+07	3.61E+08 1.73E+09]; %Live A-L+ density (/ml), one experiment/row

stotal=[3.50E+06	4.66E+06	1.50E+07	2.51E+07	3.28E+08	7.19E+08	3.64E+09
3.50E+06	4.70E+06	9.86E+06	2.64E+07	3.39E+08	7.41E+08	3.30E+09
3.50E+06	4.67E+06	1.06E+07	2.74E+07	3.34E+08	7.26E+08	2.47E+09]; %total Live density (/ml), one experiment/row

tsstep=0.5; % time step for evaluating ODE solution
tspan=0:tsstep:80;% duration for evaluating ODE solution
c0=[1.75E+06/1000000;1.75E+06/1000000;0;0];% initial condition
[t x]=ode23s(@cosmofuncVariableRelease,tspan, c0);


% this plots expt/simulation comparison
marksize=9;
figure (1)
plot(t, x(:,1), 'm:')
hold on
plot(t, x(:,2), ':','color', [0 0.75 0])
plot(t, x(:,1)+x(:,2), 'k:')
e1=errorbar(shour, mean(sred, 1)/1000000,2*std(sred, 1)/1000000,'ms-'); 
e2=errorbar(shour, mean(sgreen, 1)/1000000,2*std(sgreen, 1)/1000000,'s-', 'color', [0 0.5 0]); 
e3=errorbar(shour, mean(stotal, 1)/1000000,2*std(stotal, 1)/1000000,'ks-'); 
hold off
set(gca, 'yscale', 'log', 'ticklength', [0.025 0.05],'ytick', [1 10 1E2 1E3 1E4],'yticklabel',{'1','10','10^2','10^3','10^4'})
ylim([1 5E4])
xlim([0 80])
xlabel('Time (hours)')
ylabel('Accum. cell density (/ml)')


