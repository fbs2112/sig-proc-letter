%Volterra Data Reuse Equalyzer 

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;

numberOfSymbols = 2^numberOfBits;

e4 = cell(length(N),1);
w4 = cell(length(N),1);
% maxIt = 20;
% delayVector = 1:N+length(h);%adapFiltLength + 10;

for NIndex = 1:length(N)
    NIndex
    
    delayVector = N(NIndex)+1;%adapFiltLength + 10;
    delayVector2 = [N(NIndex)+1 N(NIndex)-2];
    e3 = cell(length(delayVector),1);
    w3 = cell(length(delayVector),1);
    


    for delay = 1:length(delayVector)


        globalLength = maxRuns + adapFiltLength(NIndex) + max(delayVector2) - 1;

        wIndex = zeros(adapFiltLength(NIndex),globalLength,maxIt);
        e2 = zeros(globalLength,maxIt);


        for index = 1:maxIt
            index

            d = zeros(globalLength,1);
            e = zeros(globalLength,1);
            G = zeros(adapFiltLength(NIndex),adapFiltLength(NIndex),globalLength);

            x = zeros(N(NIndex),globalLength);

            input = randi([0,numberOfSymbols-1],globalLength,1);

            pilot = qammod(input,numberOfSymbols,0,'gray');

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

            w = zeros(adapFiltLength(NIndex),globalLength) + 1e-6;

            channelIndex = 1;

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

                d(k) = (pilot(-delayVector(delay) + k + 1)); 

                e(k) = d(k) - w(:,k)'*xAP(:,1);


                G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength(NIndex)) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
                w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));    


            end
            wIndex(:,:,index) = conj(w(:,1:globalLength));
            e2(:,index) = abs(e).^2;
        end

        w3{delay} = mean(wIndex,3);

        e3{delay} = mean(e2,2);

    end
    w4{NIndex} = w3;
    e4{NIndex} = e3;
    
end

save(['.' filesep 'results' filesep 'results37.mat'],'w4','e4');

rmpath(['..' filesep 'simParameters' filesep]);


