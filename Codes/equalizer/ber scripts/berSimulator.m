%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers



clear;
clc;
close all;


addpath(['.' filesep 'resultsMSE']);
addpath(['.' filesep 'resultsBER']);


load results59.mat;

maxRuns = 3000;

N = 8;

adapFiltLength = (N^2+N)/2 + N;

% adapFiltLength = N;


auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);

% h = [0.34-(0.27*1i) 0.87+(0.43*1i) 0.34-(0.21*1i)]; 

% L = round((adapFiltLength + length(h))/2);

% delayVector = round((length(h) + adapFiltLength)/2);
% delayVector = 4:8;


M = 4;
numberOfSymbols = 1000;
numberOfBits = log2(M);

blockLength = numberOfSymbols*numberOfBits;

monteCarloLoops = 1000;

% equalyzerFilter = squeeze(wFinal(5,:,1)).';

h = [1 0.2 -0.3];

%  equalyzerFilter = wFinal.';
% equalyzerFilter = wFinal{:};

berAux = zeros(monteCarloLoops,1);

SNR = 0:5:30;
% SNR = 30;

signalPower = 1;
% noisePower = signalPower/db2pow(30);
% SNRAux = db2pow(30);


ber = zeros(length(SNR),size(e3,1));

for thresholdIndex = 7:size(e3,1)
%     equalyzerFilter = squeeze(wFinal(thresholdIndex,1,:));
    equalyzerFilter = squeeze(wFinal(thresholdIndex,1,:));
    for SNRIndex = 1:length(SNR)

        for j = 1:monteCarloLoops
            j
            equalyzedSignal = zeros(numberOfSymbols,1);


            binaryInputData = randi([0,1],blockLength+100,1);
            binaryInputData = reshape(binaryInputData,[],2);
            deciInputData = bi2de(binaryInputData);    
            pilot = qammod(deciInputData,2^numberOfBits,0,'gray');

            pilot = pilot.*sqrt(signalPower/var(pilot));

        %     xAux2 = filter(h,1,pilot); % linear case


        %             xAux2 = filter(h,1,pilot); % linear case

        %             xAux2(1:2,1) = pilot(1:2).^2; % nonlinear case
        %             xAux2(3,1) = pilot(3).^2 + pilot(2)*pilot(1);
        % 
        %             for m = 4:length(pilot)
        %                xAux2(m,1) = pilot(m).^2 + pilot(m-1)*pilot(m-2) + pilot(m-3);
        %             end22

    %          xAux2 = exp(-pilot);

             xAux2 = filter(h,1,pilot);
%              xAux2 = xAux2 + 0.2*(xAux2.^2);

             xAux2 = xAux2 + 0.2*(xAux2.^2) + 0.05*(xAux2.^3);
%              xAux2 = xAux2 + 0.2*(xAux2.^2);
%              xAux2 = exp(-xAux2);
        %               xAux2(1) = pilot(1);
        %               for m = 2:length(pilot)
        %                   xAux2(m,1) = pilot(m) + pilot(m-1);
        %               end


        %             xAux2(1,1) = pilot(1); % nonlinear case
        % %             xAux2(3,1) = pilot(3).^2 + pilot(2)*pilot(1);
        % 
        %             for m = 2:length(pilot)
        %                xAux2(m,1) = pilot(m) + pilot(m-1);
        %             end

        %         xAux2 = exp(pilot.*conj(pilot));


            n = randn(length(binaryInputData),1) + randn(length(binaryInputData),1)*1i;
            powerSignal = xAux2'*xAux2./(length(binaryInputData));
            powerNoiseAux = n'*n/(length(binaryInputData));
            powerNoise = (powerSignal/db2pow(SNR(SNRIndex)));
            n = n.*sqrt(powerNoise/powerNoiseAux);

            xAux = xAux2 + n;
            xAuxCorr = xAux;

            xAux = [zeros(N-1,1);xAux];

            for k = N:length(pilot)

                xFlip = xAux(k:-1:k-N+1);

                counterAux = 1;
                xTDLAux = zeros((N*N+N)/2,1);

                for lIndex = 1:length(l1)
                    xTDLAux(lIndex,1) = xFlip(l1(lIndex),1)*xFlip(l2(lIndex),1);
                end
        %         xTDLAux = [];

                xTDLConc = [xFlip;xTDLAux];

                equalyzedSignal(k,1) =  (equalyzerFilter)'*xTDLConc;


            end

        %     equalyzedSignal = filter(equalyzerFilter,1,xAux);
            [corr,lags] = xcorr(equalyzedSignal,xAux(N:end));
            [~,idx] = max(abs(corr));
            delay = abs(lags(idx));

        %     equalyzedSignal = equalyzedSignal(1:1000,:);

            decDemodSignal = qamdemod(equalyzedSignal,2^numberOfBits,0,'gray');

            binaryOutputData = de2bi(decDemodSignal,numberOfBits);


            berAux(j) = sum(sum(abs(binaryOutputData(delay+1:end,:) - binaryInputData(1:end-delay,:))))./blockLength;
        end

        ber(SNRIndex,thresholdIndex) = mean(berAux);


    end

    % semilogy(SNR,ber);

        % save(['.' filesep 'results' filesep 'results07.mat'],'wFinal','e3','meanCount');


    %     for i = 1:L+1
    %         plot(10*log10((e3(:,i))))
    %         xlabel('Iterations','interpreter','latex');
    %         ylabel('MSE (dB)','interpreter','latex');
    %         hold on;
    %     end


end


save(['.' filesep 'resultsBER' filesep 'resultsBER48.mat'],'SNR','ber');


rmpath(['.' filesep 'resultsMSE']);
rmpath(['.' filesep 'resultsBER']);

