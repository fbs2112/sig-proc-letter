%Teste Volterra SM-NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;

eta = 0.1:0.1:0.3;

numberOfSymbols = 2^numberOfBits;
windowLength = 100;

e4 = cell(length(N),length(eta));
w4 = cell(length(N),length(eta));
meanCount2 = cell(length(N),length(eta));
blindIt = zeros(maxIt,1,length(N),length(eta));
blindIt2 = zeros(maxIt,1,length(N),length(eta));

for etaIndex = 1:length(eta)
    
    for NIndex = 3:length(N)
        delayVector2 = [N(NIndex)+1 N(NIndex)-2];
        
        %     delayVector = 1:N(NIndex)+length(h);%adapFiltLength + 10;
        delayVector = N(NIndex) + 1;
        
        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector));
        
        
        for delay = 1:length(delayVector)
            delay
            globalLength = maxRuns + adapFiltLength(NIndex) + max(delayVector2) - 1;
            
            wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);
            
            count = zeros(globalLength,maxIt);
            
            w2 = zeros(adapFiltLength(NIndex),globalLength,maxIt);
            for index = 1:maxIt
                index
                
                d = zeros(globalLength,1);
                P = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
                P(:,:,adapFiltLength(NIndex) + max(delayVector2)) = eye(adapFiltLength(NIndex))*1e-6;
                sigma = zeros(globalLength,1);
                medianAux = zeros(globalLength,1);
                sigma(adapFiltLength(NIndex) + delayVector(delay)) = 1;
                delta = zeros(globalLength,1);
                lambda = zeros(globalLength,1);
                G = zeros(globalLength,1);
                
                x = zeros(N(NIndex),globalLength);
                xAP = zeros(adapFiltLength(NIndex),globalLength);
                
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
                
                theta = zeros(adapFiltLength(NIndex),globalLength);
                
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
                    
                    xAP(:,k) = [x(:,k);xTDLAux];
                    
                    
                    xAP = fliplr(xAP);
                    
                    y = theta(:,k)'*xAP(:,1);
                    
                    if (k > (adapFiltLength(NIndex) + max(delayVector2)) + windowLength && k < changingIteration) || k > changingIteration + windowLength
                        medianAux(k) = median(abs(delta(k-1 - windowLength:k-1)));
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
                    
                    delta(k) = d(k) - y;
                    
                    if abs(delta(k)) > barGamma
                        G(k) = xAP(:,k).'*P(:,:,k)*conj(xAP(:,k));
                        lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                        lambda2 = 1/lambda(k);
                        
                        %                     P(:,:,k+1) = P(:,:,k) - (lambda(k)*P(:,:,k)*conj(xAP(:,k))*xAP(:,k).'*P(:,:,k))/(1+lambda(k)*G(k));
                        
                        %                     P(:,:,k+1) = P(:,:,k) - (P(:,:,k)*lambda(k)*(conj(xAP(:,k))*xAP(:,k).')*P(:,:,k))/(1+(abs(delta(k))/barGamma)- 1);
                        
                        
                        P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(xAP(:,k))*xAP(:,k).'*P(:,:,k))/(lambda2+G(k)));
                        
                        
                        theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(xAP(:,k))*delta(k);
                        
                        sigma(k+1) = sigma(k) - (lambda(k)*delta(k)^2)/(1+lambda(k)*G(k)) + lambda(k)*delta(k)^2;
                        
                        count(k,index) = 1;
                    else
                        lambda(k) = 0;
                        P(:,:,k+1) = P(:,:,k);
                        theta(:,k+1) = theta(:,k);
                        sigma(k+1) = sigma(k);
                    end
                    
                end
                wIndex(:,:,index) = (theta(:,1:globalLength));
                e2(:,index) = abs(delta).^2;
            end
            
            meanCount{delay} = mean(count,2);
            
            w3{delay} = mean(wIndex,3);
            
            e3{delay} = mean(e2,2);
        end
        w4{NIndex,etaIndex} = w3;
        e4{NIndex,etaIndex} = e3;
        meanCount2{NIndex,etaIndex} = meanCount;
        
    end
end


save(['.' filesep 'results' filesep 'results58.mat'],'e4','w4','meanCount2', 'blindIt', 'blindIt2');

rmpath(['..' filesep 'simParameters' filesep]);

