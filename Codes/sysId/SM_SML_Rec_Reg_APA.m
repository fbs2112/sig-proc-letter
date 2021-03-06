%This script evaluates the MSE for the case of Simple Multilinear
%Functionals for a given SML model


clear;
clc;
close all;

addpath(['.' filesep 'simParameters' filesep]);

load param01.mat;
inputType = 'white';

L = 0:0;

e3 = cell(length(L));
misalignment = cell(length(L));
meanCount = zeros(length(L),1);

for LIndex = 1:length(L)
    globalLength = maxRuns + M + L(LIndex) - 1;
    e2 = zeros(globalLength,maxIt);
    wIndex = zeros(M^2,globalLength,maxIt);
    
    misalignmentAux = zeros(globalLength,maxIt);
    count = zeros(maxIt,1);
    
    for index = 1:maxIt
        
        xFlip = zeros(M,L(LIndex)+1);
        yAux = zeros(globalLength,L(LIndex)+1,K);
        yAux2 = zeros(K,globalLength,L(LIndex)+1);
        y = zeros(globalLength,L(LIndex)+1);
        Y = zeros(L(LIndex)+1,K);
        T = zeros(L(LIndex)+1,M*K);
        
        d = zeros(L(LIndex)+1,globalLength);
        e = zeros(L(LIndex)+1,globalLength);
        
        xp = zeros(globalLength,1);
        xp(M + L(LIndex),1) = 2;
        delta = zeros(globalLength,1);
        delta(M + L(LIndex),1) = 1e-3;
        
        
        invQ = zeros(L(LIndex)+1,L(LIndex)+1,globalLength);
        
        invQ(:,:,M + L(LIndex)) = 1e-6*eye(L(LIndex)+1);
        
        w = zeros(M^2,globalLength);
        
        wC = zeros(M*K,globalLength);
        
        mu = zeros(globalLength,1);
        
        index
        input = randn(globalLength,1);
        
        if strcmp(inputType,'colored')
           input = filter([1 0],[1 -0.9],input);
        end
        
        input = input.*sqrt(signalPower/var(input));

        n = randn(globalLength,1);
        n = n.*sqrt(noisePower/var(n));
               

        for j = 1:K - 1
            w1 = zeros(M,1);
            w1(1) = 2^(1-j);
        end

        w2 = zeros(M,1);

        wC(:,M+L(LIndex)) = [w1;w2];


        for k = M + L(LIndex):globalLength

            wAux = reshape(wC(:,k),[],K);

            for l = 0:L(LIndex)
                xFlip(:,l+1) = input(k-l:-1:k-M+1-l);

                for s = 1:K

                    yAux(k,l+1,s) = xFlip(:,l+1).'*wAux(:,s);

                end
            end


            for l = 0:L(LIndex)

                for s = 1:K
                    aux = yAux(k,l+1,:);
                    aux(s) = 1;

                   yAux2(s,k,l+1) = prod(aux);
                end

                y(k,l+1) = yAux2(K,k,l+1)*yAux(k,l+1,K);

                Y(l+1,:) = yAux2(:,k,l+1);

                T(l+1,:,:) = kron(Y(l+1,:).',xFlip(:,l+1));
            end


            y2 = y(k,:);

            for l = L(LIndex):-1:0
                d(l+1,k) = kron(xFlip(:,l+1),xFlip(:,l+1)).'*ho + n(k-(L(LIndex) - l));
            end
            
            

            e(:,k) = d(:,k) - y2.';
            
            absoluteValueError = abs(e(1,k));
            
            
            
            Q = (Y*Y.') .* (xFlip.'*xFlip);

            auxMatrix = Q*invQ(:,:,k);

            invQ(:,:,k+1) = 1/(1-alpha) * (invQ(:,:,k) - alpha*invQ(:,:,k)*(((alpha/(1-alpha)) *auxMatrix + eye(L(LIndex)+1))\eye(L(LIndex)+1))*1/((1-alpha))*auxMatrix);
            
            if absoluteValueError > barGamma
                                

                mu(k) = 1 - barGamma/absoluteValueError;
                
                wC(:,k+1) = wC(:,k) + mu(k)*T.'*invQ(:,:,k+1)*e(:,k);

                w(:,k+1) = kron(wC(1:end/2,k+1),wC(end/2+1:end,k+1));
                misalignmentAux(k,index) =  norm(w(:,k+1)- ho).^2/(norm(ho).^2);
                count(index) = count(index)+1; 
            else
                wC(:,k+1) = wC(:,k);
                w(:,k+1) =  w(:,k);
                misalignmentAux(k,index) = misalignmentAux(k-1,index);
            end
            
        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e(1,:)).^2;
    end
    misalignment{LIndex} = mean(misalignmentAux,2);
    w3 = mean(wIndex,3);
    meanCount(LIndex) = mean(count);

    e3{LIndex} = mean(e2,2);
end

save(['.' filesep 'results' filesep 'testSML_Rec_Reg2.mat'],'e3','w3','misalignment','meanCount');
