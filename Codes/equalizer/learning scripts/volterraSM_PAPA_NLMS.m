%%
%Teste Volterra SM-NLMS

clear;
clc;
close all;

maxRuns = 1000;
maxIt = 1000;

N = 4;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N)/2 + N;


signalPower = 1;
noisePower = 0.01;
barGamma = sqrt(5*noisePower);
gamma = 1e-12;
lambdaW = 0.98;
varW = 0.0005;

kappa = 0.5;

% wo = zeros((N^2+N)/2 + N,1);
% wo([4 5 9]) = 1;

for L = 0:0
    count = zeros(maxIt,1);

    u = zeros(L+1,1);
    u(1) = 1;
    w2 = zeros(adapFiltLength,maxRuns,maxIt);

    for j = 1:maxIt
        j
        wo(:,1) = zeros((N^2+N)/2 + N,1);
        wo([1 4],1) = 1;
%         wo([1 4],2) = [1 -1];
%         wo([1 4],3) = [-1 -1];
        
        
        input = randn(maxRuns*2,1);
%         input = filter([1 0],[1 -0.9],input);
        input = input.*sqrt(signalPower/var(input));
        
        n = randn(maxRuns*2,1);
        n = n.*sqrt(noisePower/var(n));


        xTDL = flipud(buffer(input,N,N-1,'nodelay'));

        w = zeros((N^2+N)/2 + N,maxRuns) + 0.1;
        aux = 1;

        for k = L+1:maxRuns

            for l3 = 1:L+1
                counterAux = 1;
                xTDLAux = zeros((N*N+N)/2,1);
                
                for lIndex = 1:length(l1)
                    xTDLAux(lIndex,1) = xTDL(l1(lIndex),k+l3-1-L)*xTDL(l2(lIndex),k+l3-1-L);

%                         xTDLAux(counterAux,1) = xTDL(l1,k+l3-1)*xTDL(l2,k+l3-1);

%                     counterAux = counterAux + 1;
                end
                
                
                xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1-L);xTDLAux];
%                 xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1);xTDLAux];
            end

            xAP = xTDLConc(:,k:-1:k-L);
            
            

%             d(k) = xTDLConc(1,k).^2 + xTDLConc(2,k)*xTDLConc(3,k) +
%             xTDLConc(4,k) + n(k); nonlinear system

            d(k) = ((wo(:,aux)'*xAP))^2  + n(k);
%             wo(:,k) = zeros(adapFiltLength,1);
%             wo([5 9],k) = 1;
%             wo([1 4],k) = lambdaW + n(k);
            
            e(k) = d(k) - w(:,k)'*xAP(:,1);
            absoluteValueError = abs(e(k));
            if absoluteValueError > barGamma
                mu(k) = 1 - barGamma/absoluteValueError;
                G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;
            else
                mu(k) = 0;
                count(j) = count(j)+1;
                w(:,k+1) = w(:,k);
            end

            misalignmentAux(k,j) = norm(w(:,k+1) - wo(:,aux)).^2/(norm(wo(:,aux)).^2);
            
%             if ~mod(k,N)
%                 nW = randn(N,1);
%                 nW = nW.*sqrt(varW/var(nW));
%                 wo(1:N,k+1) = lambdaW*wo(1:N,k) + nW;
%             else
%                 wo(:,k+1) = wo(:,k);
%             end
%             wo(1:N,k)
%             if ~mod(k,600)
%                 aux = aux + 1;
%             end
%             
            
        end
        w2(:,:,j) = conj(w(:,1:maxRuns));
        e2(:,j) = abs(e).^2;
    end

    meanCount(L+1) = mean(count);

    w3 = mean(w2,3);
    misalignment(:,L+1) = mean(misalignmentAux,2);
    
    e3(:,L+1) = mean(e2,2);

end
save(['.' filesep 'results' filesep 'results43.mat'],'misalignment','e3','meanCount','w3');


