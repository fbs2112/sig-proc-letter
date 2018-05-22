%Volterra Data Reuse Equalyzer

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;

eta = 0.1:0.1:0.3;

numberOfSymbols = 2^numberOfBits;

e4 = cell(length(N),length(eta));
w4 = cell(length(N),length(eta));
meanCount2 = cell(length(N),length(eta));
windowLength = 100;

blindIt = zeros(maxIt,1,length(N),length(eta));
blindIt2 = zeros(maxIt,1,length(N),length(eta));

for etaIndex = 1:length(eta)
    
    for NIndex = 3:length(N)
        NIndex
        
        delayVector = N(NIndex)+1;%adapFiltLength + 10;
        delayVector2 = [N(NIndex)+1 N(NIndex)-2];
        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector),1);
        
        for delay = 1:length(delayVector)
            
            
            globalLength = maxRuns + adapFiltLength(NIndex) + max(delayVector2) - 1;
            
            wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);
            
            count = zeros(globalLength,maxIt);
            
            for index = 1:maxIt
                index
                
                
                mu = zeros(globalLength,1);
                d = zeros(globalLength,1);
                e = zeros(globalLength,1);
                medianAux = zeros(globalLength,1);
                G = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
                
                x = zeros(N(NIndex),globalLength);
                
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
                    
                    n = randn(globalLength,1);
                    powerSignal = xAux2'*xAux2./(globalLength);
                    powerNoiseAux = n'*n/(globalLength);
                    powerNoise = (powerSignal/SNR);
                    n = n.*sqrt(powerNoise/powerNoiseAux);
                    
                    xAux(:,channelIndex) = xAux2 + n;
                    
                end
                
                w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;
                
                channelIndex = 1;
                blindFlag = 0;
                
                for k = (adapFiltLength(NIndex) + max(delayVector2)):globalLength
                    
                    if k >= changingIteration
                        if N(NIndex) > 2
                            delayVector = delayVector2(2);
                        else
                            delayVector = delayVector2(2) + 2;
                        end
                        channelIndex = 2;
                    else
                        delayVector = delayVector2(1);
                        channelIndex = 1;
                    end
                    
                    x(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);
                    
                    xTDLAux = zeros(length(l1{NIndex}),1);
                    
                    for lIndex = 1:length(l1{NIndex})
                        xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                    end
                    
                    xAP = [x(:,k);xTDLAux];
                    
                    y = w(:,k)'*xAP(:,1);
                    
                    if (k > (adapFiltLength(NIndex) + max(delayVector2)) + windowLength && k < changingIteration) || k > changingIteration + windowLength
                        medianAux(k) = median(abs(e(k-1 - windowLength:k-1)));
                        if medianAux(k) <= 2*eta(etaIndex) || blindFlag == 1
                            d(k) = pamHardThreshold(y);
                            
                            if ~blindFlag && k < changingIteration
                                blindIt(index,delay,NIndex,etaIndex) = k;
                            elseif ~blindFlag && k > changingIteration
                                blindIt2(index,delay,NIndex,etaIndex) = k;
                            end
                            blindFlag = 1;
                            
                        else
                            d(k) = (pilot(-delayVector(delay) + k + 1));
                        end
                        
                    else
                        d(k) = (pilot(-delayVector(delay) + k + 1));
                    end
                    
                    if k == changingIteration
                        blindFlag = 0;
                    end
                    
                    e(k) = d(k) - y;
                    
                    absoluteValueError = abs(e(k));
                    
                    if absoluteValueError > barGamma
                        mu(k) = 1 - barGamma/absoluteValueError;
                        G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(NIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                        w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                        count(k,index) = 1;
                    else
                        mu(k) = 0;
                        w(:,k+1) = w(:,k);
                        G(:,:,k) = eye(adapFiltLength(NIndex));
                    end
                    
                end
                wIndex(:,:,index) = conj(w(:,1:globalLength));
                e2(:,index) = abs(e).^2;
            end
            
            meanCount{delay} = mean(count,2);
            
            w3{delay} = mean(wIndex,3);
            
            e3{delay} = mean(e2,2);
            
        end
        
        meanCount2{NIndex,etaIndex} = meanCount;
        w4{NIndex,etaIndex} = w3;
        e4{NIndex,etaIndex} = e3;
        
    end
end

save(['.' filesep 'results' filesep 'results57.mat'],'w4','e4','meanCount2','blindIt', 'blindIt2');
    
rmpath(['..' filesep 'simParameters' filesep]);

