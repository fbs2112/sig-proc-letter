%Volterra Data Reuse Equalyzer 

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;

% h([5 8],1) = 0.5;


numberOfSymbols = 2^numberOfBits;

delayVector = 1:N+length(h);%adapFiltLength + 10;

delayVector = 1:N;%adapFiltLength + 10;

e3 = cell(length(delayVector),1);
w3 = cell(length(delayVector),1);
meanCount = zeros(length(delayVector),1);


for delay = 1:length(delayVector)
    
    
    globalLength = maxRuns + adapFiltLength + delayVector(delay) - 1;

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);

    count = zeros(globalLength,maxIt);
            
    for index = 1:maxIt
        index

        
        mu = zeros(globalLength,1);
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        x = zeros(N,globalLength);
        
        input = randi([0,numberOfSymbols-1],globalLength,1);

        pilot = pammod(input,pamOrder,0,'gray');

        pilot = pilot.*sqrt(signalPower/var(pilot));
        
        
        xAux = zeros(length(pilot),size(h,2));
        
        for channelIndex = 1:size(h,2)
            aux2 = zeros(length(l1Pilot),1);
            xAux2 = zeros(length(pilot),1);
        
            for i = memoryChannelLength:length(pilot) %Channel 1
               xPilot = (pilot(i:-1:i-memoryChannelLength+1));
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

        w = zeros(adapFiltLength,maxRuns) + 1e-6;
        
        channelIndex = 1;

        for k = (adapFiltLength + delayVector(delay)):globalLength
            
            if k >= changingIteration
                channelIndex = 2;
            end

            x(:,k) = xAux(k:-1:k-N+1,channelIndex);
            
            xTDLAux = zeros(length(l1),1);
            
            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = x(l1(lIndex),k)*(x(l2(lIndex),k));
            end

            xAP = [x(:,k);xTDLAux];

            d(k) = (pilot(-delayVector(delay) + k + 1)); 

            e(k) = d(k) - w(:,k)'*xAP(:,1);

            absoluteValueError = abs(e(k));

            if absoluteValueError > barGamma
                mu(k) = 1 - barGamma/absoluteValueError;
                G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));
                count(k,index) = 1;
            else
                mu(k) = 0;
                w(:,k+1) = w(:,k);
                G(:,:,k) = eye(adapFiltLength);
            end

        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end

    meanCount(delay) = mean(count,2);

    w3{delay} = mean(wIndex,3);

    e3{delay} = mean(e2,2);

end


save(['.' filesep 'results' filesep 'results05.mat'],'w3','e3','meanCount');

rmpath(['..' filesep 'simParameters' filesep]);


