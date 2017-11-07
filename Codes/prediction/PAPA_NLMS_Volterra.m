%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);
addpath(['.' filesep 'data' filesep]);


load param01.mat;
inputType = 'white';

predictingOrder = 1;


% data = csvread('BeijingPM20100101_20151231.csv',[);

data = textscan('BeijingPM20100101_20151231.csv','%s', 16);

% dataAux = xlsread('AirQualityUCI.xlsx');
% data = dataAux(:,6);

data = data(abs(data)~=200);

figure

plot(data);

globalLength = 1500;

data = buffer(data,globalLength);
% data = data./var(data);
% 
N = 10;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N)/2 + N;


% xFlip = flipud(buffer(data,N,N-1));

L = 0;%PAPA NLMS




e3 = cell(length(L));
% misalignment = cell(length(L));
powerNoise = 0.3;


% adapFiltLength = N;
mu = 0.01;
gamma = 1e-3;

monteCarloLoops = size(data,2);

% monteCarloLoops = 100;

%  adapFiltLength = N;
%  L = 0;
for LIndex = 1:length(L)
    
    %     globalLength = maxRuns + N + L(LIndex) - 1;
%     globalLength = 3000 + adapFiltLength + L(LIndex) - 1;
%     
%     
%     
%     
%     
%     
% %     data2 = (buffer(data,globalLength));
%     
% %     monteCarloLoops = size(data2,2);
%     
    wIndex = zeros(adapFiltLength,globalLength,monteCarloLoops);
    misalignmentAux = zeros(globalLength,monteCarloLoops);
    e2 = zeros(globalLength,monteCarloLoops);
    
    for index = 1:monteCarloLoops-1
        index
        
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
        x = zeros(globalLength,1);
        
        
        
        n = randn(globalLength,1);
        n = n.*sqrt(powerNoise/var(n));
        
        w = zeros(adapFiltLength,globalLength) + 0.1;
               
        x = data(:,index);
        
        for k = predictingOrder + N + L(LIndex):globalLength
%             x(k) = -0.85*x(k-1) + n(k);
            xAP = x(k-predictingOrder:-1:k-N+1-predictingOrder);
            
            xTDLAux = zeros(length(l1),1);
            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = xAP(l1(lIndex),1)*(xAP(l2(lIndex),1));
            end
            
            xAP = [xAP;xTDLAux];
            
            e(k) = x(k) - w(:,k)'*xAP;
            %
            G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
            w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k));
        end
        
       wIndex(:,:,index) = conj(w(:,1:globalLength));
       e2(:,index) = abs(e).^2;
%         for i = 1:1000
% %             n = randn(maxRuns,1);
% %             n = n*sqrt(powerNoise/var(n));
%             
%             x = zeros(maxRuns,1);
%         
% %             x(1) = 0;
% %             x(2) = 0;
%             w = zeros(N,maxRuns) + 1e-6;
% 
%             n = randn(maxRuns,1);
%             n = n*sqrt(powerNoise/var(n));
%             for k = 3:maxRuns
%                 x(k) = -0.85*x(k-1) + n(k);
%                 xAP = x(k-predictingOrder:-1:k-N+1-predictingOrder);
% %               
%                 
%                 e(k) = x(k) - w(:,k)'*xAP;
%                 %
%                 G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
%                 w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k));
% %    
%             end
%             w2(:,i) = w(:,end);
%             e2(:,i) = e;
%             
%         end
%         mse = mean(e2.^2,2);
    end
%             xTDLAux = zeros(length(l1),1);
            
            %             xInput = data2(k - predictingOrder:-1:k-N+1 - predictingOrder,index);
%             xAP = x(k-predictingOrder:-1:k-N+1-predictingOrder);
            
            %             for l3 = 1:L(LIndex)+1
            %                 for lIndex = 1:length(l1)
            %                     xTDLAux(lIndex,1) = xInput(l1(lIndex),index)*(xInput(l2(lIndex),index));
            %                 end
            %
            %                 xFlipConc = [xInput;xTDLAux];
            %
            %             end
            
%             xAP = xFlipConc;
            %
%             d(k) =   x(k+1);
%             
%             e(k) = d(k) - w(:,k)'*[x(k-1);x(k-2)];
%             
% %             G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
% %             w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k));
%             
%             w(:,k+1) = w(:,k) + 2*mu*e(k)*[x(k-1);x(k-2)];
            %             misalignmentAux(k,index) = norm(w(:,k+1) - wo(:,woIndex)).^2/(norm(wo(:,woIndex)).^2);
            
%         end
%         wIndex(:,:,index) = conj(w(:,1:globalLength));
%         e2(:,index) = abs(e).^2;
%     end
%     
    w3 = mean(wIndex,3);
%     misalignment{LIndex} = mean(misalignmentAux,2);
    
    e3{LIndex} = mean(e2(predictingOrder + N + L(LIndex):end,:),2);
    
end
% save(['.' filesep 'results' filesep 'results01.mat'],'w3','e3','misalignment');

figure
plot(10*log10(e3{1,1}));

% plot(10*log10(mean(e2,2)))
