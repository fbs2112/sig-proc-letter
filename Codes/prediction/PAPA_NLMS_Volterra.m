%Volterra SM-PAPA NLMS

clear;
clc;
close all;


addpath(['.' filesep 'simParameters' filesep]);
addpath(['.' filesep 'data' filesep]);


load param01.mat;
inputType = 'white';

predictingOrder = 5;




dataAux = xlsread('daily-maximum-temperatures-in-me.xls');

% dataAux = csvread('monthly-lake-erie-levels-1921-19.csv');
data = dataAux(15:end,2);

data = data(abs(data)~=200);

figure

plot(data);

globalLength = 400;

data = buffer(data,globalLength);
% data = data./var(data);
% 
N = 10;

auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N)/2 + N;

L = 0;%PAPA NLMS

e3 = cell(length(L));
powerNoise = 0.3;


% adapFiltLength = N;
mu = 0.1;
gamma = 1e-3;

monteCarloLoops = size(data,2);


%  adapFiltLength = N;
%  L = 0;

for i = 1:2
    
    if i == 2
        N = 65;
        adapFiltLength = N;
        mu = 0.1;
    end

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

                if i == 1
                    xAP = [xAP;xTDLAux];
                end

                e(k) = x(k) - w(:,k)'*xAP;

                ePlot(k) = (x(k) - w(:,k)'*xAP)/(abs(x(k)) + 1e-6);
                
                if isinf(ePlot(k))
                    pause;
                end
                %
                G(:,:,k) = diag(((1 - kappa*mu)/adapFiltLength) + (kappa*mu*abs(w(:,k))/norm(w(:,k),1)));
                
                G(:,:,k) = eye(adapFiltLength,adapFiltLength);
                w(:,k+1) = w(:,k) + mu*G(:,:,k)*xAP*((xAP'*G(:,:,k)*xAP+gamma*eye(L+1))\eye(L+1))*conj(e(k));
            end

           wIndex(:,:,index) = conj(w(:,1:globalLength));
           e2(:,index) = abs(e).^2;
           e2Plot(:,index) = abs(ePlot);
   
        end
         w3 = mean(wIndex,3);
    %     misalignment{LIndex} = mean(misalignmentAux,2);

        e3{LIndex} = mean(e2(predictingOrder + N + L(LIndex):end,:),2);

        e3Plot{LIndex} = mean(e2Plot(predictingOrder + N + L(LIndex):end,:),2);

    end
    % save(['.' filesep 'results' filesep 'results01.mat'],'w3','e3','misalignment');

    figure
    plot((10*log10(e3{1,1})));
%     ylim([-30 60]);
    figure
    plot(10*log10(e3Plot{1,1}));
    
%     ylim([-50 40]);

end

% plot(10*log10(mean(e2,2)))
