%Teste Volterra SM-NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;


numberOfSymbols = 2^numberOfBits;

e4 = cell(length(N),1);
w4 = cell(length(N),1);
meanCount2 = cell(length(N),1);
maxIt = 20;
signalPower = 1;

% h(:,1) = [0.5 3 5 0 0.4 0 0 1.3 0].';
% h(:,1) = h(:,2);
% h(:,2) = [1 -2.5 0 0.01 0.007 0.2 0 0 0].';

% c = h(:,2);
% h(:,2) = h(:,1);
% h(:,1) = c;

% h(:,2) = randn(9,1);

% barGamma = sqrt(5*noisePower);
for NIndex = 2:2%length(N)


%     delayVector = 1:N(NIndex)+length(h);%adapFiltLength + 10;
    delayVector = N(NIndex)+1;

    e3 = cell(length(delayVector),1);
    w3 = cell(length(delayVector),1);
    meanCount = cell(length(delayVector));
    numberOfSymbols = 4;
    pamOrder = 4;
%       = 0.05;
%     signalPower = 0.9;
    for delay = 1:length(delayVector)
        delay
        globalLength = maxRuns + adapFiltLength(NIndex) + delayVector(delay) - 1;

        wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
        e2 = zeros(globalLength,maxIt);

        count = zeros(globalLength,maxIt);

        w2 = zeros(adapFiltLength(NIndex),globalLength,maxIt);
        for index = 1:maxIt
            index

            d = zeros(globalLength,1);
            P = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);
            P(:,:,adapFiltLength(NIndex) + delayVector(delay)) = eye(adapFiltLength(NIndex))*1e3;
            sigma = zeros(globalLength,1);
            sigma(adapFiltLength(NIndex) + delayVector(delay)) = 1;
            delta = zeros(globalLength,1);
            lambda = zeros(globalLength,1);
            G = zeros(globalLength,1);

            x = zeros(N(NIndex),globalLength);
            xAP = zeros(adapFiltLength(NIndex),globalLength);

            input = randi([0,numberOfSymbols-1],globalLength,1);

            pilot = pammod(input,pamOrder,0,'gray');

            pilot = pilot.*sqrt(signalPower/var(pilot));
%             pilot = pilot./max(abs(pilot));

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
%                 xAux(:,channelIndex) = xAux(:,channelIndex) - mean(xAux(:,channelIndex));
%                 xAux(:,channelIndex) = xAux(:,channelIndex).*sqrt(1/var(xAux(:,channelIndex)));

            end

%             theta = randn(adapFiltLength(NIndex),globalLength);
            theta = zeros(adapFiltLength(NIndex),globalLength);

            channelIndex = 1;

            for k = (adapFiltLength(NIndex) + delayVector(delay)):globalLength

                if k >= changingIteration
                    channelIndex = 2;
                end

                x(:,k) = xAux(k:-1:k-N(NIndex)+1,channelIndex);

                xTDLAux = zeros(length(l1{NIndex}),1);

                for lIndex = 1:length(l1{NIndex})
                    xTDLAux(lIndex,1) = x(l1{NIndex}(lIndex),k)*(x(l2{NIndex}(lIndex),k));
                end
                
                xAP(:,k) = [x(:,k);xTDLAux];
%                 xAP(:,k) = xAP(:,k)./max(abs(xAP(:,k)));

                d(k) = (pilot(-delayVector(delay) + k + 1)); 

                delta(k) = d(k) - theta(:,k).'*xAP(:,k);
                G(k) = xAP(:,k).'*P(:,:,k)*conj(xAP(:,k));
                
               
                
                if abs(delta(k)) > barGamma || k <= adapFiltLength(NIndex) + delayVector(delay) + adapFiltLength(NIndex)
                    lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);

%                     P(:,:,k+1) = P(:,:,k) - (lambda(k)*P(:,:,k)*conj(xAP(:,k))*xAP(:,k).'*P(:,:,k))/(1+lambda(k)*G(k));
                    condV(k) = cond(P(:,:,k));

                    P(:,:,k+1) = P(:,:,k) - (P(:,:,k)*lambda(k)*conj(xAP(:,k))*xAP(:,k).'*P(:,:,k))/(1+(abs(delta(k))/barGamma)- 1);
                    if ~issymmetric(P(:,:,k+1))
                        P(:,:,k+1);
                    end
                    
                    if isnan(P(:,:,k+1))
                            P(:,:,k+1);
                    end

                    theta(:,k+1) = theta(:,k) + lambda(k)*P(:,:,k+1)*conj(xAP(:,k))*delta(k);

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
    w4{NIndex} = w3;
    e4{NIndex} = e3;
    meanCount2{NIndex} = meanCount;
    
end

% save(['.' filesep 'results' filesep 'results33.mat'],'e4','w4','meanCount2');

rmpath(['..' filesep 'simParameters' filesep]);

