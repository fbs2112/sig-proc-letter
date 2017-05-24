clear;
clc;
close all;


maxRuns = 5000; % max runs in a single independent trial
maxIt = 20;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

alpha = 0.01;      %forgetting factor of the correlation matrix in SML case

K = 2;             %number of products in the SML case
M = 10;            %length of the adaptiv filter in SML case
mu = 0.1;         %step size

h = [1 0.2 -0.3];

h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2); %unknown system in SML case


%-------------------------------------------------------------------------%
%DFE set-membership

kappa = 0.5;
gamma = 1e-12;

feedforwardLength = 4;
feedbackLength = 2;
   
adapFiltLength = feedforwardLength + feedbackLength;


barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes


numberOfBits = 2;
SNR = db2pow(30);

save(['.' filesep 'simParameters' filesep 'paramDFE.mat']);