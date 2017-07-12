%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers

clear;
clc;
close all;

addpath(['..' filesep 'learning scripts' filesep 'results']);
addpath(['..' filesep 'berParameters']);
addpath(['..' filesep 'simParameters']);

load paramDFE_FF.mat;
load param_feedforwardEq.mat;
load results47.mat;

ber = zeros(size(w4,1),size(w4,2),length(SNR));

for SNRIndex = 1:length(SNR)
   for FFIndex = 1:length(feedforwardLength)
       FFIndex
       for FBIndex = 1:length(feedbackLength)
           FBIndex
            equalyzerFilter = [];
            berAux = zeros(monteCarloLoops,1);
            equalyzerFilter(:,1) = squeeze(w4{FFIndex,FBIndex}{1}(:,4996));
            equalyzerFilter(:,2) = squeeze(w4{FFIndex,FBIndex}{1}(:,end-4));
            for index = 1:monteCarloLoops
                index
                equalyzedSignal = zeros(numberOfSymbols,1);

                binaryInputData = randi([0,1],blockLength + 100,1);
                binaryInputData = reshape(binaryInputData,[],numberOfBits);
                deciInputData = bi2de(binaryInputData);    
                pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
                pilotVar = var(pilot);
                pilot = pilot.*sqrt(signalPower/var(pilot));

                pilotAux = pilot;
    %             pilotAux = [zeros(memoryChannelLength-1,1);pilot];
                xAux = zeros(length(pilot),size(h,2));
                for channelIndex = 1:size(h,2)

                    aux2 = zeros(length(l1Pilot),1);
                    xAux2 = zeros(length(pilot),1);

                    for i = memoryChannelLength:length(pilotAux) %Channel 1
                       xPilot = (pilotAux(i:-1:i-memoryChannelLength+1));
                       for lIndex = 1:length(l1Pilot)
                          aux2(lIndex,1) = xPilot(l1Pilot(lIndex),1)*(xPilot(l2Pilot(lIndex),1));
                       end
                       xConc = [xPilot;(aux2)];
    %                    xAux2(i - memoryChannelLength + 1,1) = xConc.'*h(:,channelIndex);
                       xAux2(i - memoryChannelLength + 1,1) = xConc.'*h(:,channelIndex);
                    end

                    n = randn(length(pilot),1);
                    powerSignal = xAux2'*xAux2./(length(pilot));
                    powerNoiseAux = n'*n/(length(pilot));
                    powerNoise = (powerSignal/db2pow(SNR(SNRIndex)));
                    n = n.*sqrt(powerNoise/powerNoiseAux);

                    xAux(:,channelIndex) = xAux2 + n;

                end

                channelIndex = 1;

                xAux = [zeros(feedforwardLength(FFIndex)-1,2);xAux];
                outputFF = zeros(length(pilot),1);
                outputFB = zeros(length(pilot),1);


                for k = feedforwardLength(FFIndex) + feedbackLength(FBIndex) + 1:length(pilot) 

                    if k >= changingIteration

                        channelIndex = 2;
                    else
                        channelIndex = 1;
                    end

                    x = xAux(k:-1:k-feedforwardLength(FFIndex)+1,channelIndex);

%                     yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength(FBIndex) - 1 + 1));

                   if volterraFFFlag

                        aux = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);

                        for lIndex = 1:length(l1FF{FFIndex})
                            aux(lIndex,1) = x(l1FF{FFIndex}(lIndex),1)*(x(l2FF{FFIndex}(lIndex),1));
                        end
                        xConc = [x(:,1);aux];
                    else
                        xConc = x(:,1);
                   end
                  outputFF(k) =  (equalyzerFilter(1:adaptfiltFF(FFIndex),channelIndex))'*xConc;
                  equalyzedSignal(k,1) = pamHardThreshold(outputFF(k) + outputFB(k-1)); 

                  inputFB = equalyzedSignal(k:-1:k-feedbackLength(FBIndex) + 1); 

                    if volterraFBFlag
                        aux = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                        for lIndex = 1:length(l1FB{FBIndex})
                            aux(lIndex,1) = inputFB(l1FB{FBIndex}(lIndex),1)*(inputFB(l2FB{FBIndex}(lIndex),1));
                        end

                        yHatConc = [inputFB(:,1);aux];
                    else
                        yHatConc = inputFB(:,1);
                    end

                    if ~volterraFFFlag && ~volterraFBFlag 
                        xConc = x(:,1);
                        yHatConc = inputFB;
                    end
                    outputFB(k) = (equalyzerFilter(adaptfiltFF(FFIndex)+1:end,channelIndex))'*yHatConc;


                end
                [corr,lags] = xcorr(equalyzedSignal,xAux(feedforwardLength(FFIndex):end,1));
                [~,idx] = max(abs(corr));
                delay = abs(lags(idx));

                decDemodSignal = pamdemod(equalyzedSignal,pamOrder,0,'gray');

                binaryOutputData = de2bi(decDemodSignal,numberOfBits);

                berAux(index) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
            end

            ber(FFIndex,FBIndex,SNRIndex) = mean(berAux);
        end
   end
end


save(['.' filesep 'results' filesep 'resultsBER11.mat'],'SNR','ber');

rmpath(['..' filesep 'berParameters']);
rmpath(['..' filesep 'learning scripts' filesep 'results']);
rmpath(['..' filesep 'simParameters']);


