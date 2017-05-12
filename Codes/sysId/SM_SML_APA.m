%This script evaluates the MSE for the case of Simple Multilinear
%Functionals for a given SML model


clear;
clc;
close all;

maxRuns = 5000;
maxIt = 100;
signalPower = 1;
noisePower = 1e-3;
gamma = 1e-1;
barGamma = sqrt(5*noisePower);


K = 2;
M = 10;

mu = 0.05;

h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2);
% 
% h1 = randn(M,1);
% h2 = randn(M,1);
% ho = kron(h1,h2);

% L = 3;


L = 0:3;


for LIndex = 1:length(L)
    
    u = zeros(L(LIndex)+1,1);
    u(1) = 1;
    
    
    count = zeros(maxIt,1);

    for index = 1:maxIt
        
        xFlip = zeros(M,L(LIndex)+1);
        yAux = zeros(maxRuns,L(LIndex)+1,K);
        yAux2 = zeros(K,maxRuns,L(LIndex)+1);
        y = zeros(maxRuns,L(LIndex)+1);
        Y = zeros(L(LIndex)+1,K);
        T = zeros(L(LIndex)+1,M*K);
        
        d = zeros(L(LIndex)+1,maxRuns);
        e = zeros(L(LIndex)+1,maxRuns);
%         wC = randn(M*K,maxRuns);
        wC = zeros(M*K,maxRuns);
%         WC = ho;
        index
        input = randn(maxRuns*2,1);
        input = filter([1 0],[1 -0.9],input);
        input = input.*sqrt(signalPower/var(input));

        n = randn(maxRuns*2,1);
        n = n.*sqrt(noisePower/var(n));

        for j = 1:K - 1
            w1 = zeros(M,1);
            w1(1) = 2^(1-j);
        end

        w2 = zeros(M,1);

        wC(:,M+L(LIndex)) = [w1;w2];


        for k = M + L(LIndex):maxRuns

            w = reshape(wC(:,k),[],K);

            for l = 0:L(LIndex)
                xFlip(:,l+1) = input(k-l:-1:k-M+1-l);

                for s = 1:K

                    yAux(k,l+1,s) = xFlip(:,l+1).'*w(:,s);

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
                d(l+1,k) = kron(xFlip(:,l+1),xFlip(:,l+1)).'*ho + n(k);
            end

            e(:,k) = d(:,k) - y2.';
            
%             if abs(e(1,k))^2 > 10
%                 abs(e(1,k))^2
%                 condQ(index,k-1,LIndex)
%             end
            Q = gamma*eye(L(LIndex)+1) + (Y*Y.') .* (xFlip.'*xFlip);
%             condQ(index,k,LIndex) = cond(Q);
            absoluteValueError = abs(e(1,k));
            
            if absoluteValueError >= barGamma
                 mu(k) = 1 - barGamma/absoluteValueError;
                
                 wC(:,k+1) = wC(:,k) + mu(k)*T.'*(Q\eye(L(LIndex)+1))*e(1,k)*u;
            else
                wC(:,k+1) = wC(:,k);
                count(index) = count(index)+1; 
            end

        end
    %         w2(:,:,index) = conj(w(:,1:maxRuns));
        e2(:,index) = abs(e(1,:)).^2;
    end

     meanCount(LIndex) = mean(count);
%     w3 = mean(w2,3);
%     w4(:,1) = w3(:,end);
    e3(:,LIndex) = mean(e2,2);
 end

% % hold on
% % plot(10*log10((e3(:,2))))
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
save(['.' filesep 'resultsTest' filesep 'results07.mat'],'e3','meanCount');