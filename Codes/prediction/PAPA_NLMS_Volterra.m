%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);
addpath(['.' filesep 'data' filesep]);


load param01.mat;
inputType = 'white';

% data = csvread('AirQualityUCI.csv',2,10,[2 10 5 10]); 

dataAux = xlsread('AirQualityUCI.xlsx');
data = dataAux(:,5);
% data = data./var(data);

figure

plot(data);

N = 5;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N)/2 + N;


xFlip = flipud(buffer(data,N,N-1)); 

L = 0;%PAPA NLMS


predictingOrder = 1;


e3 = cell(length(L));
misalignment = cell(length(L));



% adapFiltLength = N;

for LIndex = 1:length(L)
    
%     globalLength = maxRuns + N + L(LIndex) - 1;
    globalLength = 300 + adapFiltLength + L(LIndex) - 1;

%     globalLength = length(data) - predictingOrder;

    

    
    
    data2 = (buffer(data,globalLength));

    monteCarloLoops = size(data2,2);
    
%         monteCarloLoops = 1;
    wIndex = zeros(adapFiltLength,globalLength,monteCarloLoops);
    misalignmentAux = zeros(globalLength,monteCarloLoops);
    e2 = zeros(globalLength,monteCarloLoops);
    
    for index = 1:1
        index
        
        d = zeros(globalLength,1);
        e = zeros(globalLength,1);
        G = zeros(adapFiltLength,adapFiltLength,globalLength);
%         
%         input = pammod(randi([0 pamOrder-1],globalLength,1),pamOrder,0,'gray');
%         
%         if strcmp(inputType,'colored')
%            input = filter([1 0],[1 -0.9],input);
%         end
%         
%         input = input.*sqrt(signalPower/var(input));
%         
%         n = randn(globalLength,1);
%         n = n.*sqrt(noisePower/var(n));

        w = zeros(adapFiltLength,globalLength) + 0.1;

%         xFlip = flipud(buffer(input,N,N-1));        
        
        
        
        xFlipConc = zeros(adapFiltLength,globalLength);
        
        woIndex = 1;
        for k = predictingOrder + N + L(LIndex):globalLength
            
            
            xTDLAux = zeros(length(l1),1);
            
            xInput = data2(k - predictingOrder:-1:k-N+1 - predictingOrder,index);
            xAP = xInput;
            
            for l3 = 1:L(LIndex)+1
                for lIndex = 1:length(l1)
                    xTDLAux(lIndex,1) = xInput(l1(lIndex),index)*(xInput(l2(lIndex),index));
                end
                
                xFlipConc = [xInput;xTDLAux];

            end
                        
            xAP = xFlipConc;
%            
            d(k) =   data(k+predictingOrder);

            e(k) = d(k) - w(:,k)'*xAP(:,1);
            
            G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
            w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k:-1:k-L(LIndex)));


%             misalignmentAux(k,index) = norm(w(:,k+1) - wo(:,woIndex)).^2/(norm(wo(:,woIndex)).^2);
            
        end
        wIndex(:,:,index) = conj(w(:,1:globalLength));
        e2(:,index) = abs(e).^2;
    end

    w3 = mean(wIndex,3);
    misalignment{LIndex} = mean(misalignmentAux,2);
    
    e3{LIndex} = mean(e2,2);

end
% save(['.' filesep 'results' filesep 'results01.mat'],'w3','e3','misalignment');

figure

plot(10*log10(e3{1,1}))
