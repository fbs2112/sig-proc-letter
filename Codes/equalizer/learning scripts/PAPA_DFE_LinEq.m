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
    
    globalLength = maxRuns + adapFiltLength + delayVector(delay) - 1 + length(h(:,1));

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);

    for index = 1:maxIt
        index

        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);

        xFiltered = zeros(globalLength,1);
        x = zeros(feedforwardLength,1);
        yHat = zeros(feedbackLength,1);

        input = randi([0,numberOfSymbols-1],globalLength,1);

        pilot = pammod(input,pamOrder,0,'gray');

        pilot = pilot.*sqrt(signalPower/var(pilot));

        xAux2 = filter(h(:,1),1,pilot);

        n = randn(globalLength,1) + randn(globalLength,1)*1i;
        powerSignal = xAux2'*xAux2./(globalLength);
        powerNoiseAux = n'*n/(globalLength);
        powerNoise = (powerSignal/SNR);
        n = n.*sqrt(powerNoise/powerNoiseAux);

        w = zeros(adapFiltLength,maxRuns) + 1e-6;
        
        hoIndex = 1;
        for k = (adapFiltLength + delayVector(delay) + length(h)):globalLength
            
            if k >= changingIteration
                hoIndex = 1;
            end
            
            xFiltered(k) = pilot(k:-1:k-length(h(:,hoIndex))+1).'*h(:,hoIndex) + n(k);

            x(:,k) = xFiltered(k:-1:k-feedforwardLength+1);

            yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength - 1 + 1));
            
            z = [x(:,k);yHat(:,k)];

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

save(['.' filesep 'results' filesep 'testPAPA_DFE_LinEq.mat'],'w3','e3');

rmpath(['..' filesep 'simParameters' filesep]);

