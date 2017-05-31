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

for LIndex = 1:length(L)
    
    globalLength = maxRuns + N + L(LIndex) - 1;

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    misalignmentAux = zeros(globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);
    
    for index = 1:maxIt
        index
        
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        
        input = pammod(randi([0 pamOrder-1],globalLength,1),pamOrder,0,'gray');
        
        if strcmp(inputType,'colored')
           input = filter([1 0],[1 -0.9],input);
        end
        
        input = input.*sqrt(signalPower/var(input));
        
        n = randn(globalLength,1);
        n = n.*sqrt(noisePower/var(n));

        w = zeros(adapFiltLength,globalLength) + 0.1;

        xFlip = flipud(buffer(input,N,N-1));        
        
        
        
        xFlipConc = zeros(adapFiltLength,globalLength);
        
        woIndex = 1;
        for k = N + L(LIndex):globalLength
            
            if k >= changingIteration
                woIndex = 2;
            end
            
            xTDLAux = zeros(length(l1),1);
            
            for l3 = 1:L(LIndex)+1
                for lIndex = 1:length(l1)
                    xTDLAux(lIndex) = xFlip(l1(lIndex),k+l3-1-L(LIndex))*xFlip(l2(lIndex),k+l3-1-L(LIndex));
                end
                
                xFlipConc(:,l3+k-L-1) = [xFlip(:,k+l3-1-L);xTDLAux];

            end
                        
            xAP = xFlipConc(:,k:-1:k-L(LIndex));
           
            d(k) = ((wo(:,woIndex)'*xAP))  + n(k);

            e(k) = d(k) - w(:,k)'*xAP(:,1);
            
            G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
            w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k:-1:k-L(LIndex)));


            misalignmentAux(k,index) = norm(w(:,k+1) - wo(:,woIndex)).^2/(norm(wo(:,woIndex)).^2);
            
        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end

    w3 = mean(wIndex,3);
    misalignment{LIndex} = mean(misalignmentAux,2);
    
    e3{LIndex} = mean(e2,2);

end
save(['.' filesep 'results' filesep 'results01.mat'],'w3','e3','meanCount');


