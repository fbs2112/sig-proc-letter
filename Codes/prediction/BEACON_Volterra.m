%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);
addpath(['.' filesep 'data' filesep]);


load param01.mat;
inputType = 'white';

predictingOrder = 1;

dataAux = xlsread('daily-maximum-temperatures-in-me.xls');

data = dataAux(15:end,2);


figure

plot(data);

globalLength = 400;

data = buffer(data,globalLength);

N = 5;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N)/2 + N;

L = 0;%PAPA NLMS

e3 = cell(length(L));
powerNoise = 0.3;


monteCarloLoops = size(data,2);
globalLength = globalLength + predictingOrder + N - 2;


barGamma = 5;


for LIndex = 1:length(L)
    
    wIndex = zeros(adapFiltLength,globalLength,monteCarloLoops);
    misalignmentAux = zeros(globalLength,monteCarloLoops);
    e2 = zeros(globalLength,monteCarloLoops);
    
    for index = 1:monteCarloLoops-1
        index
        
        P = zeros(adapFiltLength,adapFiltLength,maxIt);
        P(:,:,predictingOrder + N + L(LIndex)) = eye(adapFiltLength)*1e-3;
        sigma = zeros(globalLength,1);
        sigma(predictingOrder + N + L(LIndex)) = 1;
        delta = zeros(globalLength,1);
        lambda = zeros(globalLength,1);
        G = zeros(globalLength,1);
        %             x = zeros(globalLength,1);
        
        
        
        n = randn(globalLength,1);
        n = n.*sqrt(powerNoise/var(n));
        
        theta = zeros(adapFiltLength,globalLength) + 0.1;
        
        x = data(:,index);
        
        x = [zeros(N-1,1);x];
        
        for k = predictingOrder + N + L(LIndex):globalLength
            %             x(k) = -0.85*x(k-1) + n(k);
            xAP = x(k-predictingOrder:-1:k-N+1-predictingOrder);
            
            xTDLAux = zeros(length(l1),1);
            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = xAP(l1(lIndex),1)*(xAP(l2(lIndex),1));
            end
            
            xAP = [xAP;xTDLAux];
            
            y(k,index) = theta(:,k)'*xAP;
            
            e(k) = x(k) - y(k,index);
            
            delta(k) = x(k) - y(k,index); %error
            
            ePlot(k) = delta(k)/(abs(x(k)) + 1e-6);
            
            if abs(delta(k)) > barGamma
                G(k) = xAP.'*P(:,:,k)*conj(xAP);
                lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                lambda2 = 1/lambda(k);
                
                P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(xAP)*xAP.'*P(:,:,k))/(lambda2+G(k)));
                
                theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(xAP)*delta(k);
                
                sigma(k+1) = sigma(k) - (lambda(k)*delta(k)^2)/(1+lambda(k)*G(k)) + lambda(k)*delta(k)^2;
                
                count(k,index) = 1;
            else
                lambda(k) = 0;
                P(:,:,k+1) = P(:,:,k);
                theta(:,k+1) = theta(:,k);
                sigma(k+1) = sigma(k);
            end
        end
        
        wIndex(:,:,index) = conj(theta(:,1:globalLength));
        e2(:,index) = abs(e).^2;
        e2Plot(:,index) = abs(ePlot);
        
    end
    w3 = mean(wIndex,3);
    
    e3{LIndex} = mean(e2(predictingOrder + N + L(LIndex):end,:),2);
    
    e3Plot{LIndex} = mean(e2Plot(predictingOrder + N + L(LIndex):end,:),2);
    
    meanCount{LIndex} = mean(count(predictingOrder + N + L(LIndex):end),2);
    
end
save(['.' filesep 'results' filesep 'resultsBEACON.mat'],'e3','e3Plot','meanCount','y');

% figure
% plot((10*log10(e3{1,1})));
% %     ylim([-30 60]);
% figure
% plot(10*log10(e3Plot{1,1}));
% 
% figure
% plot(200:400,data(200:400,index-1),'*-');
% hold on
% plot(200:400,y(200:400,index-1),'r*-');
% %
% %     ylim([-50 40]);
% 
% mean(meanCount{1}(predictingOrder + N + L(LIndex):end))*100
