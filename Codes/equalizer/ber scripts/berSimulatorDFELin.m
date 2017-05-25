%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers



clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'simParameters']);


load paramDFE.mat;


load testSM_PAPA_DFE_LinEq.mat;

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
        pilot = pammod(input,pamOrder,0,'gray');
        lengthAux = length(pilot);
        
        xAux2 = filter(h,1,pilot);

        n = randn(lengthAux,1) + randn(lengthAux,1)*1i;
        powerSignal = xAux2'*xAux2./(lengthAux);
        powerNoiseAux = n'*n/(lengthAux);
        powerNoise = (powerSignal/db2pow(SNR(SNRIndex)));
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;

        xAux = [zeros(feedforwardLength-1,1);xAux];
        outputFB = zeros(length(pilot),1);
        outputFF = zeros(length(pilot),1);
        inputFB = zeros(feedbackLength,length(pilot));
        
        for k = feedforwardLength:length(pilot)

            xFlip = xAux(k:-1:k-feedforwardLength+1);
            
                
            outputFF(k) =  (equalyzerFilter(1:feedforwardLength))'*xFlip;
            equalyzedSignal(k,1) = signal2Symb(outputFF(k) + outputFB(k-1),signalPower);
            inputFB(:,k) = equalyzedSignal(k:-1:k-feedbackLength + 1); 
            outputFB(k) = (equalyzerFilter(feedforwardLength+1:end))'*inputFB(:,k);
   
        end
        [corr,lags] = xcorr(equalyzedSignal,xAux(feedforwardLength:end));
        [~,idx] = max(abs(corr));
        delay = abs(lags(idx));

        decDemodSignal = pamdemod(equalyzedSignal,pamOrder,0,'gray');

        binaryOutputData = de2bi(decDemodSignal,numberOfBits);

        berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
    end

    ber(SNRIndex) = mean(berAux);

end


save(['.' filesep 'results' filesep 'testSM_PAPA_DFELinBer.mat'],'SNR','ber');

rmpath(['..' filesep 'simParameters']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);


