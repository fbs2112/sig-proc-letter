%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);

load paramDFE_FF_FB.mat;

numberOfSymbols = 2^numberOfBits;

e4 = cell(length(feedforwardLength),length(feedforwardLength));
w4 = cell(length(feedforwardLength),length(feedforwardLength));
meanCount2 = cell(length(feedforwardLength),length(feedforwardLength));


maxIt = 10;

for FFIndex = 1:length(feedforwardLength)
    FFIndex
    for FBIndex = 1:length(feedforwardLength)
         FBIndex
%         delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

        delayVector = feedforwardLength(FFIndex)-2:feedforwardLength(FFIndex)+2;%adapFiltLength + 10;

        e3 = cell(length(delayVector),1);
        w3 = cell(length(delayVector),1);
        meanCount = cell(length(delayVector),1);
        
        for delay = 1:length(delayVector)

            globalLength = maxRuns + adapFiltLength(FFIndex,FBIndex) + delayVector(delay) - 1;

            wIndex = zeros(adapFiltLength(FFIndex,FBIndex),globalLength,maxIt);
            e2 = zeros(globalLength,maxIt);
            count = zeros(globalLength,maxIt);

            for index = 1:maxIt
                index

                mu = zeros(globalLength,1);
                d = zeros(globalLength,1);
                e = zeros(globalLength,1);
                G = zeros(adapFiltLength(FFIndex,FBIndex),adapFiltLength(FFIndex,FBIndex),globalLength);

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

                w = zeros(adapFiltLength(FFIndex,FBIndex),globalLength) + 1e-6;

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

                    e(k) = d(k) - w(:,k)'*z;

                    absoluteValueError = abs(e(k));

                    if absoluteValueError > barGamma
                        mu(k) = 1 - barGamma/absoluteValueError;
                        G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength(FFIndex,FBIndex)) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                        w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*z*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
                        count(k,index) = 1;
                    else
                        mu(k) = 0;
                        w(:,k+1) = w(:,k);
                        G(:,:,k) = eye(adapFiltLength(FFIndex,FBIndex));
                    end

                end
                wIndex(:,:,index) = conj(w(:,1:globalLength));
                e2(:,index) = abs(e).^2;
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

save(['.' filesep 'results' filesep 'testLen.mat'],'w4','e4','meanCount2');

rmpath(['..' filesep 'simParameters' filesep]);

