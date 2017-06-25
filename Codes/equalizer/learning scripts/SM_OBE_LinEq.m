%Teste Volterra SM-NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load paramEq.mat;

N = 1;

adapFiltLength = N;
maxIt = 20;
maxRuns = 1000;
numberOfSymbols = 2^numberOfBits;
c = 0.99;
d2 = 1e-12;
% signalPower = 0.01;
% barGamma = barGamma/100;

delayVector = 1:N + length(h);%adapFiltLength + 10;

e3 = cell(length(delayVector),1);
w3 = cell(length(delayVector),1);
h = [1 2 3];

% h = h./norm(h);

for delay = 1:length(delayVector)
    delay
    globalLength = maxRuns + adapFiltLength + delayVector(delay) - 1;

    wIndex = zeros(adapFiltLength,globalLength,maxIt);
    e2 = zeros(globalLength,maxIt);

    count = zeros(globalLength,maxIt);

    for index = 1:maxIt
        index
        
        d = zeros(globalLength,1);
        P = zeros(adapFiltLength,adapFiltLength,globalLength);
        P(:,:,adapFiltLength + delayVector(delay)) = eye(adapFiltLength);
        sigma = zeros(globalLength,1);
        sigma(adapFiltLength + delayVector(delay)) = 1;
        delta = zeros(globalLength,1);
        lambda = zeros(globalLength,1);
        G = zeros(globalLength,1);
        
        s = zeros(adapFiltLength,globalLength);
        
        sigmaMatrix = zeros(adapFiltLength,adapFiltLength,globalLength);
        sigmaMatrix(:,:,adapFiltLength + delayVector(delay)) = eye(adapFiltLength);
        input = randi([0,numberOfSymbols-1],globalLength,1);
        
        pilot = pammod(input,pamOrder,0,'gray');
        
        pilot = pilot.*sqrt(signalPower/var(pilot));

        xAux2 = filter(h,1,pilot);
   
        n = randn(globalLength,1);
        powerSignal = xAux2'*xAux2./(globalLength);
        powerNoiseAux = n'*n/(globalLength);
        powerNoise = (powerSignal/SNR);
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;
        
        xFlip = flipud(buffer(xAux,N,N-1));
        
%         xFlip = xAux;
        
        theta = zeros(adapFiltLength,globalLength);

        for k = (adapFiltLength + delayVector(delay)):globalLength

            xAP = xFlip(:,k);
%             xAP = sigmaMatrix(:,:,k)*xAP;
%             
            
            d(k) = (pilot(-delayVector(delay) + k + 1)); 

            delta(k) = d(k) - theta(:,k).'*xAP(:,1);
            s(:,k) = 1./(c*(xAP.*conj(xAP)) + d2);
            sigmaMatrix(:,:,k+1) = diag(s(:,k));

            
%             if(G(k) < eps)
%                 G(k);
%             end
%             
            if abs(delta(k)) > barGamma
                
                
                t(k) = abs(delta(k))/barGamma;
%                 s(:,k) = 1./(c*(xAP.*conj(xAP)) + (1-c)*s(:,k-1) + d2);
%                  sigmaMatrix(:,:,k+1) 
%                 x2 = sigmaMatrix(:,:,k)*conj(xAP);
                G(k) = xAP.'*P(:,:,k)*conj(xAP);
                lambda(k) = (1/G(k))*((abs(delta(k))/barGamma) - 1);
                
                lambda2(k) = 1/lambda(k);
                
%                 S(:,:,k+1) = c*S(:,:,k) + conj(xAP)*xAP.'*
%                 P(:,:,k+1) = (P(:,:,k) - (lambda(k)*P(:,:,k)*x2*xAP.'*P(:,:,k))/(1+lambda(k)*G(k)));


%                 P(:,:,k+1) = (P(:,:,k) - (lambda(k)*P(:,:,k)*conj(xAP)*xAP.'*P(:,:,k))/(1+lambda(k)*G(k)));
                
                P(:,:,k+1) = lambda(k)*(P(:,:,k) - (P(:,:,k)*conj(xAP)*xAP.'*P(:,:,k))/(lambda2(k)+G(k)));
                
%                 P(:,:,k+1) = triu(P(:,:,k+1)) + triu(P(:,:,k+1))' - diag(diag(triu(P(:,:,k+1))));
                
                symMatrix(index,k) = issymmetric(P(:,:,k+1));
                
%                 P(:,:,k+1) = 1/(1-c)*(P(:,:,k) - (lambda(k)*c*P(:,:,k)*conj(xAP)*xAP.'*P(:,:,k))/(1+lambda(k)*G(k)-c));

% %                 P(:,:,k+1) = P(:,:,k)  - (lambda(k)*(P(:,:,k))*conj(xAP)*xAP.'*P(:,:,k))/(t(k));
%                 a(:,:,k) = ((P(:,:,k))*conj(xAP)*xAP.'*P(:,:,k));
%                 
%                 b(k) = (1+1/(abs(delta(k))/barGamma - 1));
%                 
%                 if isnan(delta(k))
%                     delta(k);
%                 end
%                 
%                 P(:,:,k+1) = P(:,:,k)  - a(:,:,k)/(G(k)*b(k));
                

                theta(:,k+1) = theta(:,k) + P(:,:,k+1)*conj(xAP)*delta(k);
                
%                 theta(:,k+1) = theta(:,k) + lambda(k)*P(:,:,k+1)*conj(xAP)*delta(k);


                sigma(k+1) = sigma(k) - (lambda(k)*delta(k)^2)/(1+lambda(k)*G(k)) + lambda(k)*delta(k)^2;

                count(k,index) = 1;
            else
                lambda(k) = 0;
                P(:,:,k+1) = P(:,:,k);
                theta(:,k+1) = theta(:,k);
                sigma(k+1) = sigma(k);
%                 s(:,k) = s(k-1);
%                 sigmaMatrix(:,:,k+1) = sigmaMatrix(:,:,k);
            end

       end
       wIndex(:,:,index) = (theta(:,1:globalLength));
       e2(:,index) = abs(delta).^2;
    end

    meanCount{delay} = mean(count,2);

    w3{delay} = mean(wIndex,3);

    e3{delay} = mean(e2,2);
end


for i = 1:length(delayVector)
    figure
    plot(10*log10(e3{i}));
end


save(['.' filesep 'results' filesep 'testSMOBE_LinEq.mat'],'e3','w3','meanCount');

rmpath(['..' filesep 'simParameters' filesep]);

