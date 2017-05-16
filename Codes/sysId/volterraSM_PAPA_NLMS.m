%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);

load param01.mat;
inputType = 'white';

L = 0;%PAPA NLMS


e3 = cell(length(L));
misalignment = cell(length(L));
meanCount = cell(length(L));


for LIndex = 1:length(L)
    
    globalLength = maxRuns + N + L(LIndex) - 1;
    count = zeros(maxIt,1);

    u = zeros(L(LIndex)+1,1);
    u(1) = 1;
    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    misalignmentAux = zeros(globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);
    
    for index = 1:maxIt
        index
        
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        mu = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        input = randn(globalLength,1);
        
        if strcmp(inputType,'colored')
           input = filter([1 0],[1 -0.9],input);
        end
        
        input = input.*sqrt(signalPower/var(input));
        
        n = randn(globalLength,1);
        n = n.*sqrt(noisePower/var(n));

        w = zeros(adapFiltLength,globalLength) + 0.1;

        xFlip = flipud(buffer(input,N,N-1));
        
        
        for m = 1:size(xFlip,2)
            for l3 = 1:L(LIndex)+1
                xTDLAux = zeros(adapFiltLength - N,size(xFlip,2));

                for lIndex = 1:length(l1)
                    xTDLAux(lIndex,m) = xFlip(l1(lIndex),m+l3-1-L(LIndex))*xFlip(l2(lIndex),m+l3-1-L(LIndex));
                end

            end

        end
        
        xFlipConc = [xFlip;xTDLAux];
        
        for k = N + L(LIndex):globalLength

            xAP = xFlipConc(:,k:-1:k-L(LIndex));
           
            d(k) = ((wo(:,1)'*xAP))  + n(k);

            e(k) = d(k) - w(:,k)'*xAP(:,1);
            absoluteValueError = abs(e(k));
            
            if absoluteValueError > barGamma
                mu(k) = 1 - barGamma/absoluteValueError;
                G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;
                count(k,index) = 1;
            else
                w(:,k+1) = w(:,k);
            end

            misalignmentAux(k,index) = norm(w(:,k+1) - wo(:,1)).^2/(norm(wo(:,1)).^2);
            
        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end

    meanCount{LIndex} = mean(count,2);

    w3 = mean(wIndex,3);
    misalignment{LIndex} = mean(misalignmentAux,2);
    
    e3{LIndex} = mean(e2,2);

end
save(['.' filesep 'results' filesep 'resultsPAPATest.mat'],'misalignment','e3','meanCount','w3');


