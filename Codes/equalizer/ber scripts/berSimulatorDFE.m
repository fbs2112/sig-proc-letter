%This scripts generates bit error rates using pre trained adaptive filters
%as equalyzers



clear;
clc;
close all;


addpath(['.' filesep 'resultsMSE']);
addpath(['.' filesep 'resultsBER']);


load results52.mat;

maxRuns = 3000;

volterraFFFlag = 0;
volterraFBFlag = 1;


barGammaVector = 1;
feedforwardLength = 4;
feedbackLength = 2;

adaptfiltFF = (feedforwardLength^2+feedforwardLength)/2 + feedforwardLength;
adaptfiltFB = (feedbackLength^2+feedbackLength)/2 + feedbackLength;

adaptfilt = adaptfiltFF + adaptfiltFB;

auxMatrix = triu(ones(feedforwardLength));
[l1FF,l2FF] = find(auxMatrix);

auxMatrix = triu(ones(feedbackLength));
[l1FB,l2FB] = find(auxMatrix);


h = [1 0.2 -0.3];

% delayVector = 3;
delayVector = round((length(h) + feedforwardLength)/2);

if ~volterraFFFlag
    adaptfiltFF = feedforwardLength;
end

if ~volterraFBFlag
    adaptfiltFB = feedbackLength;
end

    

adapFiltLength = adaptfiltFF + adaptfiltFB;



% N = 8;
% 
% adapFiltLength = (N^2+N)/2 + N;
% 
% % adapFiltLength = N;
% 
% 
% auxMatrix = triu(ones(N));
% [l1,l2] = find(auxMatrix);

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


%  equalyzerFilter = wFinal.';
% equalyzerFilter = wFinal{:};

berAux = zeros(monteCarloLoops,1);

SNR = 0:5:30;
% SNR = 30;

signalPower = 1;
% noisePower = signalPower/db2pow(30);
% SNRAux = db2pow(30);


ber = zeros(length(SNR),size(e3,1));

for thresholdIndex = 1:size(e3,1)
%     equalyzerFilter = squeeze(wFinal(thresholdIndex,1,:));
    equalyzerFilter = squeeze(wFinal(1,1,:));
    for SNRIndex = 1:length(SNR)

        for j = 1:monteCarloLoops
            j
            equalyzedSignal = zeros(numberOfSymbols,1);


            binaryInputData = randi([0,1],blockLength+100,1);
            binaryInputData = reshape(binaryInputData,[],2);
            deciInputData = bi2de(binaryInputData);    
            pilot = qammod(deciInputData,2^numberOfBits,0,'gray');
            pilotPower = pilot'*pilot/(length(pilot));

            pilot = pilot.*sqrt(signalPower/pilotPower);

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

            xAux = [zeros(feedforwardLength-1,1);xAux];
            outputFB = zeros(length(pilot),1);
            
%             yAux = zeros(adaptfiltFB,1);

            for k = feedforwardLength:length(pilot)

                xFlip = xAux(k:-1:k-feedforwardLength+1);

                 if volterraFFFlag
                    
                    aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                    for lIndex = 1:length(l1FF)
                        aux(lIndex,1) = xFlip(l1FF(lIndex),1)*(xFlip(l2FF(lIndex),1));
                    end
                    xConc = [xFlip(:,1);aux];
                else
                    xConc = xFlip(:,1);
                 end

                
                 if ~volterraFFFlag && ~volterraFBFlag 
                    xConc = x(:,k);
                 end

                
        %         xTDLAux = [];

%                 xTDLConc = [xFlip;xTDLAux];

                outputFF(k) =  (equalyzerFilter(1:adaptfiltFF))'*xConc;
                equalyzedSignal(k,1) = signal2Symb(outputFF(k) + outputFB(k-1),signalPower);
%                 outputFFSymb(k) = qammod(qamdemod(outputFF(k-1),2^numberOfBits,0,'gray'),2^numberOfBits,0,'gray');
%                 outputFFSymb(k) = signal2Symb(outputFF(k-1),signalPower);
                inputFB(:,k) = equalyzedSignal(k:-1:k-feedbackLength + 1); 
                
                
                
                if volterraFBFlag
                    aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                    for lIndex = 1:length(l1FB)
                        aux(lIndex,1) = inputFB(l1FB(lIndex),k)*(inputFB(l2FB(lIndex),k));
                    end

                    yHatConc = [inputFB(:,k);aux];
                else
                    yHatConc = inputFB(:,k);
                end

                if ~volterraFFFlag && ~volterraFBFlag 
                    yHatConc = inputFB(:,k);
                end
                
                outputFB(k) = (equalyzerFilter(adaptfiltFF+1:end))'*yHatConc;
                
%                 equalyzedSignal(k,1) = outputFF(k) + outputFB;
                

            end

        %     equalyzedSignal = filter(equalyzerFilter,1,xAux);
            [corr,lags] = xcorr(equalyzedSignal,xAux(feedforwardLength:end));
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


save(['.' filesep 'resultsBER' filesep 'resultsBER42.mat'],'SNR','ber');


rmpath(['.' filesep 'resultsMSE']);
rmpath(['.' filesep 'resultsBER']);

