% Script for estimating the parameters in the multi-stage model using proliferation data of FUCCI-expressing cells
%
% Author: Sean T. Vittadello
%         sean.vittadello@qut.edu.au
%         School of Mathematical Sciences
%         Queensland University of Technology
%
% Last update: 3 June 2019
%
% Associated function: MultiStageModel_fit_function.m


%% General input
% Text interpreter
set(0,'defaultTextInterpreter','latex');

% Microscopy time-series images
dt = 15/60; % Time interval between images in time-series stack (h).
slices = 192; % Number of slices in time-series stack.
T = dt*(slices-1); % Total experiment duration (first slice is at t = 0)


%% Input parameters
% Cell line
Cdata = 18; % Experimental mean cell-cycle duration (h)
            % C8161 cell line - 18 h
            % WM983C cell line - 27 h
            % 1205Lu cell line - 36 h

% Multi-stage model
Sr = 18; % Number of stages in G1.
Sy = 18; % Number of stages in eS.
Sg = 18; % Number of stages in S/G2/M.
k = Sr + Sy + Sg; % Total number of stages in the cell cycle.
Ir = 9;  % Number of consecutive stages in G1 for which the initial
         % cell numbers are equal. **Must be a divisor of Sr.**
Iy = 9;  % Number of consecutive stages in eS for which the initial
         % cell numbers are equal. **Must be a divisor of Sy.**
Ig = 9;  % Number of consecutive stages in S/G2/M for which the initial
         % cell numbers are equal. **Must be a divisor of Sg.**

% Objective function
param = [0.1+(1-0.1)*rand(1,round(Sr/Ir+Sy/Iy+Sg/Ig)) 4+(8-4)*rand(1,3)]; % Specify the starting parameters for the objective function [Sr/Ir parameters for cells in the G1 stages, Sy/Iy parameters for cells in the eS stages, Sg/Ig parameters for cells in the S/G2/M stages, Mean duration of cells in G1, Mean duration of cells in eS, Mean duration of cells in S/G2/M]
weights = [0.0001 0.0001 0.0001 0.001 0 0]; % Specify the weights [w2 w3 w4 w5 w6 w7]

% Location of data file
path = '..\C8161_data_Figure1.csv';


%% Create required data sets
Cryg = csvread(path,1,1);
Cryg = Cryg';

RatioData = Cryg(1,:)./(Cryg(2,:) + Cryg(3,:)); % Ratio of the number of cells in G1 to the number of cells in eS and S/G2/M, at each time point
Data_r = Cryg(1,:); % Number of cells in G1 at each time point
Data_y = Cryg(2,:); % Number of cells in eS at each time point
Data_g = Cryg(3,:); % Number of cells in S/G2/M at each time point
Data_tot = sum(Cryg(:,:),1); % Total number of cells at each time point


%% Optimisation of the multi-stage model
func = @(x)MultiStageModel_fit_function(x,RatioData,Data_r,Data_y,Data_g,Data_tot,Cdata,weights,T,Sr,Sy,Sg,k,Ir,Iy,Ig);
options = optimoptions(@lsqnonlin,'MaxIter',1000000000000);
[optim_param,resnorm,exitflag] = lsqnonlin(func,param,zeros(1,size(param,2)),[],options);


%% Multi-stage model evaluation with optimised starting parameters
% Initial parameters
Lr = optim_param(end-2); % Mean duration of G1 phase
Ly = optim_param(end-1); % Mean duration of eS phase
Lg = optim_param(end); % Mean duration of S/G2/M phase
C = Lr + Ly + Lg; % Mean cell cycle duration
lambda = [(Sr/Lr)*ones(1,Sr),(Sy/Ly)*ones(1,Sy),(Sg/Lg)*ones(1,Sg)]; % Transition rates between stages; set equal within each phase

I_optim = zeros(1,k); % Initial cell density in each stage at time t
for i = 1 : round(Sr/Ir)
    I_optim(1,((i-1)*Ir+1):i*Ir) = optim_param(i);
end
for i = 1 : round(Sy/Iy)
    I_optim(1,(Sr+(i-1)*Iy+1):(Sr+i*Iy)) = optim_param(round(Sr/Ir)+i);
end
for i = 1 : round(Sg/Ig)
    I_optim(1,((Sr+Sy+(i-1)*Ig+1)):(Sr+Sy+i*Ig)) = optim_param(round(Sr/Ir + Sy/Iy)+i);
end

% Numerical solution by forward Euler, using the optimised initial parameters
h = 0.05; % Step size
N = T/h; % Number of steps

M = zeros(N+1,k); % Cell numbers in each stage at time t
M(1,:) = I_optim; % Initial cell numbers

for i = 1:N
    M(i+1,1) = M(i,1) + h*(2*lambda(1,k)*M(i,k)-lambda(1,1)*M(i,1));
    for j = 2:k-1
        M(i+1,j) = M(i,j) + h*(lambda(1,j-1)*M(i,j-1)-lambda(1,j)*M(i,j));
    end
    M(i+1,k) = M(i,k) + h*(lambda(1,k-1)*M(i,k-1)-lambda(1,k)*M(i,k));
end

% Display the optimisation residual, estimated initial numbers of cells in each phase, the estimated phase durations, and estimated cell cycle duration
[resnorm sum(M(1,1:Sr),2) sum(M(1,(Sr+1):(Sr+Sy)),2) sum(M(1,(Sr+Sy+1):k),2) Lr Ly Lg C]


%% Plot the data and the optimised model solution
% Ratio Q(t) - data
figure(1)
tdata = 0:dt:(slices-1)*dt;
plot(tdata,Cryg(1,:)./(Cryg(2,:)+Cryg(3,:)),'Marker','.','MarkerSize',10,'Color','b');
xlim([0 48]);
xticks(0:12:48);
xlabel('Time, $t$ (h)');
% ylim([0.3 1.2]);
% yticks([0.5 1]);
ylabel('$Q(t)$');
set(gca,'FontSize',28)
hold on
% Ratio Q(t) - optimised model solution
t = 0:h:T;
plot(t,(sum(M(:,1:Sr),2)./sum(M(:,Sr+1:k),2))','Marker','.','MarkerSize',1,'LineWidth',1,'Color','r');
legend({'Data','Multi-stage model'},'Location','SouthEast','FontSize',16)

% Subpopulations R(t), Y(t) and G(t) - data
figure(2)
plot(tdata,Cryg(1,:),'Marker','.','MarkerSize',10,'Color','r');
xlim([0 48]);
xticks(0:12:48);
xlabel('Time, $t$ (h)');
% ylim([0 1300]);
% yticks([0,500,1000]);
ylabel('$R(t)$, $Y(t)$, $G(t)$');
set(gca,'FontSize',28)
hold on
plot(tdata,Cryg(2,:),'Marker','.','MarkerSize',10,'Color','y');
plot(tdata,Cryg(3,:),'Marker','.','MarkerSize',10,'Color','g');
% Subpopulations R(t), Y(t) and G(t) - optimised model solution
plot(t,sum(M(:,1:Sr),2),'LineWidth',2,'Color','r');
plot(t,sum(M(:,(Sr+1):(Sr+Sy)),2),'LineWidth',2,'Color','y');
plot(t,sum(M(:,(Sr+Sy+1):k),2),'LineWidth',2,'Color','g');
legend({'Data - G1','Data - eS','Data - S/G2/M','Multi-stage model - G1','Multi-stage model - eS','Multi-stage model - S/G2/M'},'Location','NorthWest','FontSize',16)

% Total population M(t) - data
figure(3)
plot(tdata,sum(Cryg(:,:),1),'Marker','.','MarkerSize',10,'Color','b');
hold on
% Total population M(t) - optimised model solution
plot(t,sum(M(1:end,:),2),'LineWidth',2,'Color','k');
% Total population M(t) - exponential model fit
exp_fit = fit(tdata',sum(Cryg(:,:),1)','exp1');
exp_plot = plot(exp_fit);
set(exp_plot,'LineWidth',2);
xlim([0 48]);
xticks(0:12:48);
xlabel('Time, $t$ (h)');
% ylim([0 3000]);
% yticks([0 1500 3000]);
ylabel('$M(t)$');
set(gca,'FontSize',28)
legend({'Data','Multi-stage model','Exponential model'},'Location','NorthWest','FontSize',16);
% R-squared from linear regression of ln(M(t)) versus t
lnMt = log(sum(Cryg(:,:),1));
line_coeff = polyfit(tdata,lnMt,1);
line_values = polyval(line_coeff,tdata);
yresid = lnMt - line_values;
SSresid = sum(yresid.^2);
SStotal = (length(lnMt)-1) * var(lnMt);
rsq = 1 - SSresid/SStotal; % R-squared
text(0.85,0.1,['$R^2$ = ',num2str(round(rsq,2))],'Units','Normalized','Interpreter','Latex','FontSize',14);
