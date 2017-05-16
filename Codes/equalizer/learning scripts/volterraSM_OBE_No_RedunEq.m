%Teste Volterra SM-NLMS

clear;
clc;
close all;


addpath(['..' filesep 'simParameters' filesep]);


load param.mat;


delayVector = 3;%adapFiltLength + 10;


for delay = 1:length(delayVector)
    delay

    count = zeros(maxIt,1);

    w2 = zeros(adapFiltLength,maxRuns,maxIt);
    for j = 1:maxIt
        j
        
        P = zeros(adapFiltLength,adapFiltLength,maxRuns*2);
        sigma = zeros(maxRuns*2,1);
        delta = zeros(maxRuns,1);
        G = zeros(maxRuns*2,1);
        lambda = zeros(maxRuns*2,1);

        P(:,:,adapFiltLength + delayVector(delay) + 10) = eye(adapFiltLength)*1e-6;
        sigma(adapFiltLength + delayVector(delay) + 10) = 1;


        input = randi([0,3],maxRuns*2,1);
        pilot = qammod(input,4);

        pilot = pilot.*sqrt(signalPower/var(pilot));

        xAux2 = filter(h,1,pilot);
        xAux2 = xAux2 + 0.2*(xAux2.^2) + 0.05*(xAux2.^3);
   
        n = randn(maxRuns*2,1) + randn(maxRuns*2,1)*1i;
        powerSignal = xAux2'*xAux2./(maxRuns*2);
        powerNoiseAux = n'*n/(maxRuns*2);
        powerNoise = (powerSignal/SNRAux);
        n = n.*sqrt(powerNoise/powerNoiseAux);

        xAux = xAux2 + n;


        xTDL = flipud(buffer(input,N,N-1,'nodelay'));

        theta = zeros((N^2+N)/2 + N,maxRuns);

    %             w = zeros(N,maxRuns);


        for k = (adapFiltLength + delayVector(delay) + 10):maxRuns
            xAP = xAux(k:-1:k-N+1);

            xTDLAux = zeros((N*N+N)/2,1);

            for lIndex = 1:length(l1)
                xTDLAux(lIndex,1) = xAP(l1(lIndex))*(xAP(l2(lIndex)));
            end

            xTDLConc = [xAP;xTDLAux];

            d(k) = (pilot(-delayVector(delay) + k + 1)); 
    %             d(k) = (pilot(-L2 + adapFiltLength+1 - L)); 

            delta(k) = d(k) - theta(:,k).'*xTDLConc(:,1);
            G(k) = xTDLConc.'*P(:,:,k)*conj(xTDLConc);
    %                if abs(delta(k)) <= gamma
            if abs(delta(k)) <= gamma

                lambda(k) = 0;
                count(j) = count(j) + 1;
                P(:,:,k+1) = P(:,:,k);
                theta(:,k+1) = theta(:,k);
                sigma(k+1) = sigma(k); 

            else
                lambda(k) = (1/G(k))*((abs(delta(k))/gamma) - 1);
                 P(:,:,k+1) = (P(:,:,k)) - (lambda(k)*(P(:,:,k))*conj(xTDLConc)*xTDLConc.'*(P(:,:,k)))/(1+lambda(k)*G(k));

                theta(:,k+1) = theta(:,k) + lambda(k)*P(:,:,k+1)*conj(xTDLConc)*delta(k);

                sigma(k+1) = sigma(k) - (lambda(k)*delta(k)^2)/(1+lambda(k)*G(k)) + lambda(k)*delta(k)^2;

            end



       end
       w2(:,:,j) = (theta(:,1:maxRuns));
       e2(:,j) = abs(delta).^2;
    end

    meanCount = mean(count);

    w3 = mean(w2,3);

    wFinal(delay,:,1) = w3(:,end);

    e3(delay,:,1) = mean(e2,2);


    % save(['.' filesep 'results' filesep 'results07.mat'],'wFinal','e3','meanCount');


    %     for i = 1:L+1
    %         plot(10*log10((e3(:,i))))
    %         xlabel('Iterations','interpreter','latex');
    %         ylabel('MSE (dB)','interpreter','latex');
    %         hold on;
%     end

end

save(['.' filesep 'resultsMSE' filesep 'results24.mat'],'e3','wFinal','meanCount');



