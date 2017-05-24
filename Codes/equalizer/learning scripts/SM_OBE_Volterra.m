%Teste Volterra SM-NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;


numberOfSymbols = 2^numberOfBits;

delayVector = 1:N+length(h);%adapFiltLength + 10;

e3 = cell(length(delayVector),1);
w3 = cell(length(delayVector),1);


for delay = 1:length(delayVector)
    delay
    globalLength = maxRuns + adapFiltLength + delayVector(delay) - 1;

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);

    count = zeros(globalLength,maxIt);

    w2 = zeros(adapFiltLength,globalLength,maxIt);
    for index = 1:maxIt
        index
        
        d = zeros(globalLength,1);
        P = zeros(adapFiltLength,adapFiltLength,maxIt);
        P(:,:,adapFiltLength + delayVector(delay)) = eye(adapFiltLength)*1e-6;
        sigma = zeros(globalLength,1);
        sigma(adapFiltLength + delayVector(delay)) = 1;
        delta = zeros(globalLength,1);
        lambda = zeros(globalLength,1);
        G = zeros(globalLength,1);
        

        input = randi([0,numberOfSymbols-1],globalLength,1);
        
        pilot = qammod(input,numberOfSymbols,0,'gray');
        
        pilot = pilot.*sqrt(signalPower/var(pilot));

        xAux2 = filter(h,1,pilot);
        xAux2 = xAux2 + 0.2*(xAux2.^2);
   
        n = randn(globalLength,1) + randn(globalLength,1)*1i;
        powerSignal = xAux2'*xAux2./(globalLength);
        powerNoiseAux = n'*n/(globalLength);
        powerNoise = (powerSignal/SNR);
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;
        
        xFlip = flipud(buffer(xAux,N,N-1));

        theta = zeros(adapFiltLength,globalLength);

        for k = (adapFiltLength + delayVector(delay)):globalLength

            xTDLAux = zeros(adapFiltLength - N,1);

            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = xFlip(l1(lIndex),k)*(xFlip(l2(lIndex),k));
            end

            xAP = [xFlip(:,k);xTDLAux];
        
            d(k) = (pilot(-delayVector(delay) + k + 1)); 

            delta(k) = d(k) - theta(:,k).'*xAP(:,1);
            G(k) = xAP.'*P(:,:,k)*conj(xAP);
            
            if abs(delta(k)) > barGamma
                lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);

                P(:,:,k+1) = P(:,:,k) - (lambda(k)*P(:,:,k)*conj(xAP)*xAP.'*P(:,:,k))/(1+lambda(k)*G(k));

                theta(:,k+1) = theta(:,k) + lambda(k)*P(:,:,k+1)*conj(xAP)*delta(k);

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

    meanCount = mean(count,2);

    w3{delay} = mean(wIndex,3);

    e3{delay} = mean(e2,2);
end

save(['.' filesep 'results' filesep 'testOBE_Eq.mat'],'e3','w3','meanCount');

rmpath(['..' filesep 'simParameters' filesep]);
