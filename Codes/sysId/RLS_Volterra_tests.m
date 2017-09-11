%Volterra SM-OBE

clear;
clc;
close all;

addpath(['.' filesep 'simParameters' filesep]);

load param01.mat;
inputType = 'white';

lambda = 0.95;
maxIt = 20;
globalLength = maxRuns + N - 1;

misalignmentAux = zeros(globalLength,maxIt);
wIndex = zeros(adapFiltLength,globalLength,maxIt);
e2 = zeros(globalLength,maxIt);

for index = 1:maxIt
    index

    d = zeros(globalLength,1);
    Sd = zeros(adapFiltLength,adapFiltLength,globalLength);
    Sd(:,:,N-1) = eye(adapFiltLength)*(1-lambda)*signalPower;
    psi = zeros(adapFiltLength,globalLength);
    e = zeros(globalLength,1);
    
    input = pammod(randi([0 pamOrder-1],globalLength,1),pamOrder,0,'gray');
    
    if strcmp(inputType,'colored')
       input = filter([1 0],[1 -0.9],input);
    end
    
    input = input.*sqrt(signalPower/var(input));

    n = randn(globalLength,1);
    n = n.*sqrt(noisePower/var(n));

    w = zeros((N^2+N)/2 + N,globalLength);    
    
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

        e(k) = d(k) - w(:,k-1).'*xAP;
        
        psi(:,k) = Sd(:,:,k-1)*xAP;
       
        Sd(:,:,k) = (1/lambda)*(Sd(:,:,k-1)-((psi(:,k)*psi(:,k).')/(lambda+psi(:,k).'*xAP)));
        
        w(:,k) = w(:,k-1) + conj(e(k))*Sd(:,:,k)*xAP;
       
       
        
        misalignmentAux(k,index) = norm(w(:,k) - wo(:,woIndex)).^2/(norm(wo(:,woIndex)).^2);
        
    end
    wIndex(:,:,index) = (w(:,1:globalLength));
    e2(:,index) = abs(e).^2;
end

w3 = mean(wIndex,3);

e3 = mean(e2,2);

misalignment = mean(misalignmentAux,2);

save(['.' filesep 'results' filesep 'RLS.mat'],'w3','e3','misalignment');

rmpath(['.' filesep 'simParameters' filesep]);

