%Volterra Data Reuse Equalyzer 

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;

adapFiltLength = N;

numberOfSymbols = 2^numberOfBits;

delayVector = 1:N+length(h);%adapFiltLength + 10;

e3 = cell(length(delayVector),1);
w3 = cell(length(delayVector),1);

for delay = 1:length(delayVector)
    
    
    globalLength = maxRuns + adapFiltLength + delayVector(delay) - 1;

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);
            
    for index = 1:maxIt
        index

        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        input = randi([0,numberOfSymbols-1],globalLength,1);

        pilot = pammod(input,pamOrder,0,'gray');

        pilot = pilot.*sqrt(signalPower/var(pilot));

        xAux2 = filter(h,1,pilot);

        n = randn(globalLength,1) + randn(globalLength,1)*1i;
        powerSignal = xAux2'*xAux2./(globalLength);
        powerNoiseAux = n'*n/(globalLength);
        powerNoise = (powerSignal/SNR);
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;

        xFlip = flipud(buffer(xAux,N,N-1));

        w = zeros(adapFiltLength,globalLength) + 1e-1;

       for k = (adapFiltLength + delayVector(delay)):globalLength

            xAP = xFlip(:,k);

            d(k) = (pilot(-delayVector(delay) + k + 1)); 

            e(k) = d(k) - w(:,k)'*xAP(:,1);

            G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
            w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(1))\eye(1))*conj(e(k));            

        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end


    w3{delay} = mean(wIndex,3);

    e3{delay} = mean(e2,2);

end


save(['.' filesep 'results' filesep 'testPAPA_LinEq.mat'],'w3','e3');

rmpath(['..' filesep 'simParameters' filesep]);


