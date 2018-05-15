%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers

clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'berParameters']);
addpath(['..' filesep 'simParameters']);
load paramEq.mat;
load param_feedforwardEq.mat;
load results53.mat;

changingIteration = (blockLength + 100)/2;

ber = zeros(size(w4,1),length(SNR));

for SNRIndex = 1:length(SNR)
    for NIndex = 3:size(w4,1)
        equalizerFilter = [];
        berAux = zeros(monteCarloLoops,1);
        equalizerFilter(:,1) = squeeze(w4{NIndex}{1}(:,end));
        equalizerFilter(:,2) = squeeze(w4{NIndex}{1}(:,end-4));
       
        for index = 1:monteCarloLoops
            index
            equalizedSignal = zeros(numberOfSymbols,1);

            binaryInputData = randi([0,1],blockLength + 100,1);
            binaryInputData = reshape(binaryInputData,[],numberOfBits);
            deciInputData = bi2de(binaryInputData);    
            pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
            pilotVar = var(pilot);
            pilot = pilot.*sqrt(signalPower/var(pilot));
            
            pilotAux = pilot;
%             pilotAux = [zeros(memoryChannelLength-1,1);pilot];
            xAux = zeros(length(pilot),size(h,2));
            pilot2 = pilot;
            globalLength = length(pilot2);
%
            for channelIndex = 1:size(h,2)
                aux2 = zeros(length(l1Pilot),1);
                xAux2 = zeros(length(pilot),1);

                for i = memoryChannelLength:length(pilot) %Channel 1
                   xPilot = (pilot2(i:-1:i-memoryChannelLength+1));
                   for lIndex = 1:length(l1Pilot)
                      aux2(lIndex,1) = xPilot(l1Pilot(lIndex),1)*(xPilot(l2Pilot(lIndex),1));
                   end
                   xConc = [xPilot;(aux2)];
                   xAux2(i-memoryChannelLength+1,1) = xConc.'*h(:,channelIndex);
                end


    %         n = randn(globalLength,1) + randn(globalLength,1)*1i;

                n = randn(globalLength,1);
                powerSignal = xAux2'*xAux2./(globalLength);
                powerNoiseAux = n'*n/(globalLength);
                powerNoise = (powerSignal/db2pow(SNR(SNRIndex)));
                n = n.*sqrt(powerNoise/powerNoiseAux);

                xAux(:,channelIndex) = xAux2 + n;

            end

            channelIndex = 1;

            xAux = [zeros(N(NIndex)-1,2);xAux];
    
            for k = N(NIndex):length(pilot) 
                
                x = xAux(k:-1:k-N(NIndex)+1,channelIndex);

                xTDLAux = zeros(length(l1{NIndex}),1);

                for lIndex = 1:length(l1{NIndex})
                    xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex))*(x(l2{NIndex}(lIndex)));
                end

                xAP = [x;xTDLAux];

                equalizedSignal(k,1) =  equalizerFilter(:,channelIndex)'*xAP;

            end
            [corr,lags] = xcorr(equalizedSignal,xAux(N(NIndex):end,1));
            [~,idx] = max(abs(corr));
            delay = abs(lags(idx));
            decDemodSignal = pamdemod(equalizedSignal,pamOrder,0,'gray');
            binaryOutputData = de2bi(decDemodSignal,numberOfBits);

            berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
        end

        ber(NIndex,SNRIndex) = mean(berAux);
    end

end


save(['.' filesep 'results' filesep 'resultsBER17.mat'],'SNR','ber');

rmpath(['..' filesep 'berParameters']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);
rmpath(['..' filesep 'simParameters']);


