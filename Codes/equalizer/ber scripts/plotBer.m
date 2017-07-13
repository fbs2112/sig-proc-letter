clear;
clc;
close all;


addpath(['.' filesep 'results']);



load resultsBER02.mat;

for i = 1:size(ber,1)
    for j = 1:size(ber,2)
        
        figure
        semilogy(SNR,squeeze(ber(i,j,:)))
    end
end

