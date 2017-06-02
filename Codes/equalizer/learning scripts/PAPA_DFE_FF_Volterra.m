%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);

load paramDFE_FF.mat;

% h = [1 -2.5 1 0 0 0 0 0 0].';
% h(4:end,:) = 0;
% h([5 8],1) = 0.75;

numberOfSymbols = 2^numberOfBits;

delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

delayVector = 6;

e3 = cell(length(delayVector),1);
w3 = cell(length(delayVector),1);


for delay = 1:length(delayVector)
    
    globalLength = maxRuns + adapFiltLength + delayVector(delay) - 1 + length(h(:,1));

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);

    for index = 1:maxIt
        index

        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);

        xFiltered = zeros(globalLength,1);
        x = zeros(feedforwardLength,globalLength);
        yHat = zeros(feedbackLength,1);

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

            x(:,k) = xAux(k:-1:k-feedforwardLength+1,channelIndex);

            yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength - 1 + 1));

            if volterraFFFlag

                aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                for lIndex = 1:length(l1FF)
                    aux(lIndex,1) = x(l1FF(lIndex),k)*(x(l2FF(lIndex),k));
                end
                xConc = [x(:,k);aux];
            else
                xConc = x(:,k);
            end


            if volterraFBFlag
                aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                for lIndex = 1:length(l1FB)
                    aux(lIndex,1) = yHat(l1FB(lIndex),k)*(yHat(l2FB(lIndex),k));
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

            e(k) = d(k) - w(:,k)'*z;

            G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
            w(:,k+1) = w(:,k) + mu*G(:,:,k)*z*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));

        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end

    w3{delay} = mean(wIndex,3);

    e3{delay} = mean(e2,2);

end

save(['.' filesep 'results' filesep 'results03.mat'],'w3','e3');

rmpath(['..' filesep 'simParameters' filesep]);

