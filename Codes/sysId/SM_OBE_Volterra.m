%Volterra SM-OBE

clear;
clc;
close all;

addpath(['.' filesep 'simParameters' filesep]);

load param01.mat;
inputType = 'white';


globalLength = maxRuns + N - 1;
count = zeros(globalLength,maxIt);

misalignmentAux = zeros(globalLength,maxIt);
wIndex = zeros(adapFiltLength,globalLength,maxIt);
e2 = zeros(globalLength,maxIt);

for index = 1:maxIt
    index

    d = zeros(globalLength,1);
    P = zeros(adapFiltLength,adapFiltLength,maxIt);
    P(:,:,N) = eye(adapFiltLength)*1e-6;
    sigma = zeros(globalLength,1);
    sigma(N) = 1;
    delta = zeros(globalLength,1);
    lambda = zeros(globalLength,1);
    G = zeros(globalLength,1);
    
    input = pammod(randi([0 pamOrder-1],globalLength,1),pamOrder,0,'gray');
    
    if strcmp(inputType,'colored')
       input = filter([1 0],[1 -0.9],input);
    end
    
    input = input.*sqrt(signalPower/var(input));

    n = randn(globalLength,1);
    n = n.*sqrt(noisePower/var(n));

    theta = zeros((N^2+N)/2 + N,globalLength);    
    
    xFlip = flipud(buffer(input,N,N-1));
    woIndex = 1;
    
    for k = N:globalLength
        
        if k >= changingIteration
            woIndex = 2;
        end

        
        xTDLAux = zeros(adapFiltLength - N,1);
        
        for lIndex = 1:length(l1)
            xTDLAux(lIndex,1) = xFlip(l1(lIndex),k)*xFlip(l2(lIndex),k);
        end
        
        xAP = [xFlip(:,k);xTDLAux];
        
        d(k) = ((wo(:,woIndex).'*xAP)) + n(k);

        delta(k) = d(k) - theta(:,k)'*xAP; %error
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
        
        
        misalignmentAux(k,index) = norm(theta(:,k+1) - wo(:,woIndex)).^2/(norm(wo(:,woIndex)).^2);
        
    end
    wIndex(:,:,index) = (theta(:,1:globalLength));
    e2(:,index) = abs(delta).^2;
end

meanCount = mean(count,2);

w3 = mean(wIndex,3);

e3 = mean(e2,2);

misalignment = mean(misalignmentAux,2);

save(['.' filesep 'results' filesep 'resultsSM_0BE_Test.mat'],'misalignment','e3','meanCount','w3');


