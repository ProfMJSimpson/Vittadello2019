% Function for estimating the parameters in the multi-stage model using proliferation data of FUCCI-expressing cells
%
% Author: Sean T. Vittadello
%         sean.vittadello@qut.edu.au
%         School of Mathematical Sciences
%         Queensland University of Technology
%
% Last update: 3 June 2019
%
% Associated script: MultiStageModel_fit.m


%% F is the objective function to optimise
function F = MultiStageModel_fit_function(param,RatioData,Data_r,Data_y,Data_g,Data_tot,Cdata,weights,T,Sr,Sy,Sg,k,Ir,Iy,Ig)


%% Inputs
% param - Starting parameters for the objective function
% RatioData - Ratio of the number of cells in G1 to the number of cells in eS and S/G2/M, at each time point
% Data_r - Number of cells in G1 at each time point
% Data_y - Number of cells in eS at each time point
% Data_g - Number of cells in S/G2/M at each time point
% Data_tot - Total number of cells at each time point
% Cdata - Experimental mean cell-cycle duration (h)
% weights - Weights for the components of the objective function F
% T - Total experiment duration
% Sr - Number of stages in G1
% Sy - Number of stages in eS
% Sg - Number of stages in S/G2/M
% k - Total number of stages in the cell cycle
% Ir - Number of consecutive stages in G1 for which the initial cell numbers are equal
% Iy - Number of consecutive stages in eS for which the initial cell numbers are equal
% Ig - Number of consecutive stages in S/G2/M for which the initial cell numbers are equal


%% Outputs
% F - Objective function evaluated with the multi-stage model solution from the starting parameters given in param

%% Multi-stage model solution using starting parameters
% Initial parameters
Lr = param(end-2); % Mean duration of G1 phase
Ly = param(end-1); % Mean duration of eS phase
Lg = param(end); % Mean duration of S/G2/M phase
lambda = [(Sr/Lr)*ones(1,Sr),(Sy/Ly)*ones(1,Sy),(Sg/Lg)*ones(1,Sg)]; % Transition rates between stages

% Initial cell numbers in each stage
M_init = zeros(1,k);
for i = 1 : round(Sr/Ir)
    M_init(1,((i-1)*Ir+1):i*Ir) = param(i);
end
for i = 1 : round(Sy/Iy)
    M_init(1,(Sr+(i-1)*Iy+1):(Sr+i*Iy)) = param(round(Sr/Ir)+i);
end
for i = 1 : round(Sg/Ig)
    M_init(1,((Sr+Sy+(i-1)*Ig+1)):(Sr+Sy+i*Ig)) = param(round(Sr/Ir + Sy/Iy)+i);
end

% Numerical solution by forward Euler
h = 0.05; % Step size
N = T/h; % Number of steps

M = zeros(N+1,k); % Cell numbers in each stage at time t
M(1,:) = M_init; % Initial cell numbers

for i = 1:N
    M(i+1,1) = M(i,1) + h*(2*lambda(1,k)*M(i,k)-lambda(1,1)*M(i,1));
    for j = 2:k-1
        M(i+1,j) = M(i,j) + h*(lambda(1,j-1)*M(i,j-1)-lambda(1,j)*M(i,j));
    end
    M(i+1,k) = M(i,k) + h*(lambda(1,k-1)*M(i,k-1)-lambda(1,k)*M(i,k));
end


%% Evaluation of the objective function
% Create required variables from the multi-stage model solution
Model_r = sum(M(:,1:Sr),2)'; % Number of cells in G1 at each time point
Model_y = sum(M(:,(Sr+1):(Sr+Sy)),2)'; % Number of cells in eS at each time point
Model_g = sum(M(:,(Sr+Sy+1):k),2)'; % Number of cells in S/G2/M at each time point
RatioModel = Model_r./(Model_y + Model_g); % Ratio of the number of cells in G1 to the number of cells in eS and S/G2/M, at each time point
Model_tot = sum(M(:,:),2)'; % Total number of cells at each time point

% Evaluate the components of the objective function
f1 = RatioData - RatioModel(1,1:5:end); % Difference between the data and model values of the ratio
f2 = Data_r - Model_r(1,1:5:end); % Difference between the data and model values of the number of cells in G1
f3 = Data_y - Model_y(1,1:5:end); % Difference between the data and model values of the number of cells in eS
f4 = Data_g - Model_g(1,1:5:end); % Difference between the data and model values of the number of cells in S/G2/M
f5 = Data_g(1,end) - Model_g(1,end); % Difference between the data and model values of the number of cells in S/G2/M at the final time point
f6 = Cdata-Lr-Ly-Lg; % Difference between the data and model values of the cell cycle duration
f7 = Data_tot - Model_tot(1,1:5:end); % Difference between the data and model values of the total number of cells

% Objective function
F = [f1 weights(1)*f2 weights(2)*f3 weights(3)*f4 weights(4)*f5 weights(5)*f6 weights(6)*f7];


end