clear;
clc;
close all;


maxRuns = 8000; % max runs in a single independent trial
maxIt = 50;    %number of independent trial
signalPower = 1;    %desired input signal power
noisePower = 1e-3;  %desired measurement noise power

alpha = 0.01;      %forgetting factor of the correlation matrix in SML case

K = 2;             %number of products in the SML case
M = 10;            %length of the adaptiv filter in SML case
mu = 0.1;         %step size

h(:,1) = [1 -2.5 1 0 2 0 0 0.7 0].';
h(:,2) = [0.5 3 0 0.5 0.001 0.3 0 0 0].';

h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2); %unknown system in SML case


%-------------------------------------------------------------------------%
%DFE set-membership

kappa = 0.5;
gamma = 1e-12;


memoryChannelLength = 3;

volterraFFFlag = 1;
volterraFBFlag = 0;


feedforwardLength = 6;
feedbackLength = 6;

adaptfiltFF = (feedforwardLength^2+feedforwardLength)/2 + feedforwardLength;
adaptfiltFB = (feedbackLength^2+feedbackLength)/2 + feedbackLength;

adaptfilt = adaptfiltFF + adaptfiltFB;

auxMatrix = triu(ones(feedforwardLength));
[l1FF,l2FF] = find(auxMatrix);

auxMatrix = triu(ones(feedbackLength));
[l1FB,l2FB] = find(auxMatrix);

auxMatrix = triu(ones(memoryChannelLength));
[l1Pilot,l2Pilot] = find(auxMatrix);


if ~volterraFFFlag
    adaptfiltFF = feedforwardLength;
end

if ~volterraFBFlag
    adaptfiltFB = feedbackLength;
end

   
adapFiltLength = adaptfiltFF + adaptfiltFB;


barGamma = 4*sqrt(5*noisePower); %threshold for set-membership purposes


numberOfBits = 2;
changingIteration = 4000;

pamOrder = 2^numberOfBits;

SNR = db2pow(30);

save(['.' filesep 'simParameters' filesep 'paramDFE_FF.mat']);