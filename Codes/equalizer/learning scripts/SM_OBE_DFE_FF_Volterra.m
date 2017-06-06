%Volterra SM OBE DFE

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);

load paramDFE_FF.mat;


numberOfSymbols = 2^numberOfBits;


e4 = cell(length(feedforwardLength),length(feedbackLength));
w4 = cell(length(feedforwardLength),length(feedbackLength));
meanCount2 = cell(length(feedforwardLength),length(feedbackLength));


for FFIndex = 1:length(feedforwardLength)
    FFIndex
    for FBIndex = 1:length(feedbackLength)
         FBIndex
%         delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

        delayVector = 1:feedforwardLength(FFIndex)+length(h);%adapFiltLength + 10;
        
        
        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector),1);
        

        for delay = 1:length(delayVector)     

            delay
            globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + delayVector(delay) - 1;

            wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);

            count = zeros(globalLength,maxIt);

            for index = 1:maxIt
                index

                d = zeros(globalLength,1);
                P = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),maxIt);
                P(:,:,adapFiltLength(FFIndex,FBIndex) + delayVector(delay)) = eye(adapFiltLength(FFIndex,FBIndex))*1e-6;
                sigma = zeros(globalLength,1);
                sigma(adapFiltLength(FFIndex,FBIndex) + delayVector(delay)) = 1;
                delta = zeros(globalLength,1);
                lambda = zeros(globalLength,1);
                G = zeros(globalLength,1);

                x = zeros(feedforwardLength(FFIndex),globalLength);
                yHat = zeros(feedbackLength(FBIndex),1);

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

                theta = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

                channelIndex = 1;

                for k = (adapFiltLength(FFIndex,FBIndex) + delayVector(delay)):globalLength

                    if k >= changingIteration
                        channelIndex = 2;
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

                    delta(k) = d(k) - theta(:,k).'*z;

                    G(k) = z.'*P(:,:,k)*conj(z);

                    if abs(delta(k)) > barGamma
                        lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);

                        P(:,:,k+1) = P(:,:,k) - (lambda(k)*P(:,:,k)*conj(z)*z.'*P(:,:,k))/(1+lambda(k)*G(k));

                        theta(:,k+1) = theta(:,k) + lambda(k)*P(:,:,k+1)*conj(z)*delta(k);

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
        
        meanCount2{FFIndex,FBIndex} = meanCount;
        w4{FFIndex,FBIndex} = w3;
        e4{FFIndex,FBIndex} = e3;
        
    end
    
end

save(['.' filesep 'results' filesep 'results35.mat'],'w4','e4','meanCount2');

rmpath(['..' filesep 'simParameters' filesep]);
