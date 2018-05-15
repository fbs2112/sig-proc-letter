clear;
clc;
close all;


maxRuns = 7000; % max runs in a single independent trial
maxIt = 1000;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

alpha = 0.01;      %forgetting factor of the correlation matrix in SML case

K = 2;             %number of products in the SML case
M = 10;            %length of the adaptiv filter in SML case
mu = 0.8;         %step size

% h(:,1) = [1 -2.5 1 0 2 0 0 0.7 0].';
% h(:,2) = [0.5 3 0 0.5 0.001 0.3 0 0 0].';

h(:,1) = [0.5 3 5 0 0.3 0 0 1.2 0].';
% h(:,2) = [0.5 3 0 0.5 0.001 0.3 0 0 0].';
h(:,2) = [2 0 0 0.2 0.3 -0.7 0 0 0].'; 

% h(:,2) = [2 0 2 0 -0.7 0 0 3 0].';


h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2); %unknown system in SML case


%-------------------------------------------------------------------------%
%Volterra set-membership

N = 7:9;

l1 = cell(length(N),1);
l2 = cell(length(N),1);
adapFiltLength = zeros(length(N),1);

for i = 1:length(N)
    auxMatrix = triu(ones(N(i)));
    [l1{i},l2{i}] = find(auxMatrix);
    adapFiltLength(i) = (N(i)^2+N(i))/2 + N(i);
end

kappa = 0.5;
gamma = 1e-12;

barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes

numberOfBits = 2;

pamOrder = 2^numberOfBits;

memoryChannelLength = 3;

auxMatrix = triu(ones(memoryChannelLength));
[l1Pilot,l2Pilot] = find(auxMatrix);


changingIteration = 10000;


SNR = db2pow(30);

save(['.' filesep 'simParameters' filesep 'paramEq.mat']);