%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers



clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'simParameters']);


load paramEq.mat;



load testSM_PAPA_Volterra.mat;

numberOfSymbols = 1000;

blockLength = numberOfSymbols*numberOfBits;

monteCarloLoops = 1000;

berAux = zeros(monteCarloLoops,1);

SNR = 0:5:30;

ber = zeros(length(SNR),1);

chosenDelay = 1;

equalyzerFilter = squeeze(w3{chosenDelay}(:,end));


for SNRIndex = 1:length(SNR)

    for index = 1:monteCarloLoops
        index
        equalyzedSignal = zeros(numberOfSymbols,1);

        binaryInputData = randi([0,1],blockLength + 100,1);
        binaryInputData = reshape(binaryInputData,[],numberOfBits);
        deciInputData = bi2de(binaryInputData);    
        pilot = qammod(deciInputData,2^numberOfBits,0,'gray');
        lengthAux = length(pilot);
        
        xAux2 = filter(h,1,pilot);

        n = randn(lengthAux,1) + randn(lengthAux,1)*1i;
        powerSignal = xAux2'*xAux2./(lengthAux);
        powerNoiseAux = n'*n/(lengthAux);
        powerNoise = (powerSignal/db2pow(SNR(SNRIndex)));
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;

        xAux = [zeros(N-1,1);xAux];

        for k = N:length(pilot)
            
            xFlip = xAux(k:-1:k-N+1);
            xTDLAux = zeros(adapFiltLength - N,1);

            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = xFlip(l1(lIndex),1)*(xFlip(l2(lIndex),1));
            end

            xAP = [xFlip;xTDLAux];

            equalyzedSignal(k,1) =  (equalyzerFilter)'*xAP;

        end
        [corr,lags] = xcorr(equalyzedSignal,xAux(N:end));
        [~,idx] = max(abs(corr));
        delay = abs(lags(idx));

        decDemodSignal = qamdemod(equalyzedSignal,2^numberOfBits,0,'gray');

        binaryOutputData = de2bi(decDemodSignal,numberOfBits);

        berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
    end

    ber(SNRIndex) = mean(berAux);

end


save(['.' filesep 'results' filesep 'testSM_PAPA_VolterraEqBer.mat'],'SNR','ber');

rmpath(['..' filesep 'simParameters']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);


