%Volterra RLS DFE FF and FB


clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);

load paramDFE_FF.mat;

lambda = 0.95;

numberOfSymbols = 2^numberOfBits;

e4 = cell(length(feedforwardLength),length(feedbackLength));
w4 = cell(length(feedforwardLength),length(feedbackLength));
% maxIt = 20;

for FFIndex = 1:length(feedforwardLength)
    FFIndex
    for FBIndex = 1:length(feedbackLength)
        FBIndex
        
        delayVector = feedforwardLength(FFIndex)+1;
        
        
        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        
        
        for delay = 1:length(delayVector)
            
            delayVector2 = [feedforwardLength(FFIndex)+1 feedforwardLength(FFIndex)-2];
            globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + max(delayVector2) - 1;
            wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);
            
            for index = 1:maxIt
                index
                
                d = zeros(globalLength,1);
                Sd = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);
                Sd(:,:,(adapFiltLength(FFIndex,FBIndex) + max(delayVector2)) -1) = eye(adapFiltLength(FFIndex,FBIndex))*(1-lambda)*signalPower;
                psi = zeros(adapFiltLength(FFIndex,FBIndex),globalLength);
                
                e = zeros(globalLength,1);
                
                xFiltered = zeros(globalLength,1);
                x = zeros(feedforwardLength(FFIndex),globalLength);
                yHat = zeros(feedbackLength(FBIndex),1);
                
                input = randi([0,numberOfSymbols-1],globalLength,1);
                
                pilot = pammod(input,pamOrder,0,'gray');
                
                
                pilot2 = pilot.*sqrt(signalPower/var(pilot));
                
                
                
                xAux = zeros(length(pilot),size(h,2));
                
                for channelIndex = 1:size(h,2)
                    aux2 = zeros(length(l1Pilot),1);
                    xAux2 = zeros(length(pilot),1);
                    
                    for i = memoryChannelLength:length(pilot) %Channel 1
                        xPilot = (pilot2(i:-1:i-memoryChannelLength+1));
                        for lIndex = 1:length(l1Pilot)
                            aux2(lIndex,1) = xPilot(l1Pilot(lIndex),1)*(xPilot(l2Pilot(lIndex),1));
                        end
                        xConc = [xPilot;(aux2)];
                        xAux2(i,1) = xConc.'*h(:,channelIndex);
                    end
                    
                    
                    %         n = randn(globalLength,1) + randn(globalLength,1)*1i;
                    
                    n = randn(globalLength,1);
                    powerSignal = xAux2'*xAux2./(globalLength);
                    powerNoiseAux = n'*n/(globalLength);
                    powerNoise = (powerSignal/SNR);
                    n = n.*sqrt(powerNoise/powerNoiseAux);
                    
                    xAux(:,channelIndex) = xAux2 + n;
                    
                end
                
                w = zeros(adapFiltLength(FFIndex,FBIndex),globalLength);
                
                channelIndex = 1;
                
                for k = (adapFiltLength(FFIndex,FBIndex) + max(delayVector2)):globalLength
                    
                    if k >= changingIteration
                        if feedforwardLength(FFIndex) > 2
                            delayVector = delayVector2(2);
                        else
                            delayVector = delayVector2(2) + 2;
                        end
                        channelIndex = 2;
                    else
                        delayVector = delayVector2(1);
                        channelIndex = 1;
                    end
                    
                    x(:,k) = xAux(k:-1:k-feedforwardLength(FFIndex)+1,channelIndex);
                    
                    yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength(FBIndex) - 1 + 1));
                    
                    if volterraFFFlag
                        
                        aux = zeros((feedforwardLength(FFIndex)^2+feedforwardLength(FFIndex))/2,1);
                        
                        for lIndex = 1:length(l1FF{FFIndex})
                            aux(lIndex,1) = x(l1FF{FFIndex}(lIndex),k)*(x(l2FF{FFIndex}(lIndex),k));
                        end
                        xConc = [x(:,k);aux];
                    else
                        xConc = x(:,k);
                    end
                    
                    
                    if volterraFBFlag
                        aux = zeros((feedbackLength(FBIndex)^2+feedbackLength(FBIndex))/2,1);
                        for lIndex = 1:length(l1FB{FBIndex})
                            aux(lIndex,1) = yHat(l1FB{FBIndex}(lIndex),k)*(yHat(l2FB{FBIndex}(lIndex),k));
                        end
                        
                        yHatConc = [yHat(:,k);aux];
                    else
                        yHatConc = yHat(:,k);
                    end
                    
                    if ~volterraFFFlag && ~volterraFBFlag
                        xConc = x(:,k);
                        yHatConc = yHat(:,k);
                    end
                    
                    
                    z = [xConc;yHatConc];
                    
                    d(k) = (pilot(-delayVector(delay) + k + 1));
                    
                    e(k) = d(k) - w(:,k-1)'*z;
                    
                    psi(:,k) = Sd(:,:,k-1)*z;
                    
                    Sd(:,:,k) = (1/lambda)*(Sd(:,:,k-1)-((psi(:,k)*psi(:,k).')/(lambda+psi(:,k).'*z)));
                    
                    w(:,k) = w(:,k-1) + conj(e(k))*Sd(:,:,k)*z;
                    
                end
                wIndex(:,:,index) = conj(w(:,1:globalLength));
                e2(:,index) = abs(e).^2;
            end
            
            w3{delay} = mean(wIndex,3);
            
            e3{delay} = mean(e2,2);
            
        end
        w4{FFIndex,FBIndex} = w3;
        e4{FFIndex,FBIndex} = e3;
    end
end

save(['.' filesep 'results' filesep 'results51.mat'],'w4','e4');

rmpath(['..' filesep 'simParameters' filesep]);
