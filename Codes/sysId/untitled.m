clear;
clc;
close all;




M = 4;
channel = [1 1];

SNR = 30; 
signal = randi([0 M-1],1000,1);

sentSignal = qammod(signal,M,0);

plot(sentSignal,'o')

figure;

scatterplot(sentSignal)

signalVar = var(sentSignal);

noise = randn(1000,1) + randn(1000,1)*1j;

noisePower = signalVar/db2pow(SNR);

noise = noise*sqrt(noisePower/var(noise));

corruptedSignal = filter(channel,1,sentSignal) + noise;

figure

scatterplot(corruptedSignal)