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


% adapFiltLength = N;
mu = 0.1;
gamma = 1e-3;
lambda = 0.95;

monteCarloLoops = size(data,2);

globalLength = globalLength + predictingOrder + N - 2;

%  adapFiltLength = N;
%  L = 0;

% for i = 1:2

%     if i == 2
%         N = 65;
%         adapFiltLength = N;
%         mu = 0.1;
%     end

for LIndex = 1:length(L)
    
    wIndex = zeros(adapFiltLength,globalLength,monteCarloLoops);
    misalignmentAux = zeros(globalLength,monteCarloLoops);
    e2 = zeros(globalLength,monteCarloLoops);
    
    for index = 1:monteCarloLoops-1
        index
        
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        %             x = zeros(globalLength,1);
        Sd = zeros(adapFiltLength,adapFiltLength,globalLength);
        Sd(:,:,predictingOrder + N + L(LIndex)-1) = eye(adapFiltLength)*(1-lambda)*signalPower;
        psi = zeros(adapFiltLength,globalLength);
        
        
        
        n = randn(globalLength,1);
        n = n.*sqrt(powerNoise/var(n));
        
        w = zeros(adapFiltLength,globalLength) + 0.1;
        
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
            
            y(k,index) = w(:,k-1)'*xAP;
            
            e(k) = x(k) - y(k,index);
            
            
            ePlot(k) = e(k)/(abs(x(k)) + 1e-6);
            
            
            psi(:,k) = Sd(:,:,k-1)*xAP;
            
            Sd(:,:,k) = (1/lambda)*(Sd(:,:,k-1)-((psi(:,k)*psi(:,k).')/(lambda+psi(:,k).'*xAP)));
            
            w(:,k) = w(:,k-1) + conj(e(k))*Sd(:,:,k)*xAP;
        end
        
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
        e2Plot(:,index) = abs(ePlot);
        
    end
    w3 = mean(wIndex,3);
    
    e3{LIndex} = mean(e2(predictingOrder + N + L(LIndex):end,:),2);
    
    e3Plot{LIndex} = mean(e2Plot(predictingOrder + N + L(LIndex):end,:),2);
    
end
save(['.' filesep 'results' filesep 'resultsRLS.mat'],'e3','e3Plot','y');

%     figure
%     plot((10*log10(e3{1,1})));
% %     ylim([-30 60]);
%     figure
%     plot(10*log10(e3Plot{1,1}));
%
%     figure
%     plot(200:400,data(200:400,index-1),'*-');
%     hold on
%     plot(200:400,y(200:400,index-1),'r*-');
% %
%     ylim([-50 40]);

% end
