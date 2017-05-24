%Volterra NLMS DFE

clear;
clc;
close all;

addpath(['..' filesep 'simParameters' filesep]);

load paramDFE.mat;


numberOfSymbols = 2^numberOfBits;

delayVector = 1:feedforwardLength+length(h);%adapFiltLength + 10;

e3 = cell(length(delayVector),1);
w3 = cell(length(delayVector),1);


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

        x = zeros(feedforwardLength,1);
        yHat = zeros(feedbackLength,1);

        input = randi([0,numberOfSymbols-1],globalLength,1);

        pilot = qammod(input,numberOfSymbols,0,'gray');

        pilot = pilot.*sqrt(signalPower/var(pilot));

        xAux2 = filter(h,1,pilot);

        xAux2 = xAux2;

        n = randn(globalLength,1) + randn(globalLength,1)*1i;
        powerSignal = xAux2'*xAux2./(globalLength);
        powerNoiseAux = n'*n/(globalLength);
        powerNoise = (powerSignal/SNR);
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;

        w = zeros(adapFiltLength,maxRuns) + 1e-6;

        for k = (adapFiltLength + delayVector(delay)):globalLength

            x(:,k) = xAux(k:-1:k-feedforwardLength+1);

            yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength - 1 + 1));
            
            z = [x(:,k);yHat(:,k)];

            d(k) = (pilot(-delayVector(delay) + k + 1));

            e(k) = d(k) - w(:,k)'*z;

            absoluteValueError = abs(e(k));

            if absoluteValueError > barGamma
                mu(k) = 1 - barGamma/absoluteValueError;
                G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*z*((z'*G(:,:,k)*z+gamma*eye(1))\eye(1))*conj(e(k));
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

    meanCount = mean(count,2);

    w3{delay} = mean(wIndex,3);

    e3{delay} = mean(e2,2);

end

save(['.' filesep 'results' filesep 'testSM_PAPA_DFE_LinEq.mat'],'w3','e3','meanCount');

rmpath(['..' filesep 'simParameters' filesep]);

