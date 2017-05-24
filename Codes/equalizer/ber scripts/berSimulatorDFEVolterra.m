%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers



clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'simParameters']);


load paramDFE_FF_FB.mat;


load testSM_PAPA_DFE_Volterra.mat;

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

        xAux = [zeros(feedforwardLength-1,1);xAux];
        outputFB = zeros(length(pilot),1);
        outputFF = zeros(length(pilot),1);
        inputFB = zeros(feedbackLength,length(pilot));
        
        for k = feedforwardLength:length(pilot)

            xFlip = xAux(k:-1:k-feedforwardLength+1);
            
            if volterraFFFlag
                    
                aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                for lIndex = 1:length(l1FF)
                    aux(lIndex,1) = xFlip(l1FF(lIndex),1)*(xFlip(l2FF(lIndex),1));
                end
                xConc = [xFlip(:,1);aux];
            else
                xConc = xFlip(:,1);
            end


            if ~volterraFFFlag && ~volterraFBFlag 
               xConc = x(:,k);
            end

                
            outputFF(k) =  (equalyzerFilter(1:adaptfiltFF))'*xConc;
            equalyzedSignal(k,1) = signal2Symb(outputFF(k) + outputFB(k-1),signalPower);
            inputFB(:,k) = equalyzedSignal(k:-1:k-feedbackLength + 1); 
            
            
            if volterraFBFlag
                aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                for lIndex = 1:length(l1FB)
                    aux(lIndex,1) = inputFB(l1FB(lIndex),k)*(inputFB(l2FB(lIndex),k));
                end

                yHatConc = [inputFB(:,k);aux];
            else
                yHatConc = inputFB(:,k);
            end

            if ~volterraFFFlag && ~volterraFBFlag 
                yHatConc = inputFB(:,k);
            end

            
            outputFB(k) = (equalyzerFilter(adaptfiltFF+1:end))'*yHatConc;
   
        end
        [corr,lags] = xcorr(equalyzedSignal,xAux(feedforwardLength:end));
        [~,idx] = max(abs(corr));
        delay = abs(lags(idx));

        decDemodSignal = qamdemod(equalyzedSignal,2^numberOfBits,0,'gray');

        binaryOutputData = de2bi(decDemodSignal,numberOfBits);

        berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
    end

    ber(SNRIndex) = mean(berAux);

end


save(['.' filesep 'results' filesep 'testSM_PAPA_DFEVolterraBer.mat'],'SNR','ber');

rmpath(['..' filesep 'simParameters']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);


