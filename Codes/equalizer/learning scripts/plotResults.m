clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep 'simParameters' filesep]);
addpath(['..' filesep '..' filesep]);
addpath(['..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);



linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

fileVector = 1:12;

meanCountVector = 5:12;

auxIndex = 1;

for i = 1:length(fileVector)
    
    if i > 9
        load(['results' num2str(fileVector(i)) '.mat']);
    else
        load(['results0' num2str(fileVector(i)) '.mat']);
    end
    
%     load teste.mat
%     for delay = 1:size(e3,1)
%         figure
%         for l = 1:size(e3,2)
%             aux = find(e3{delay},1);
% 
%             plot(10*log10((e3{delay}(aux:end))))
%             hold on;
%         end
% 
%         xlabel('Iterations','interpreter','latex');
%         ylabel('MSE (dB)','interpreter','latex');
% 
%         H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%         set(H,'interpreter','latex')
%     %     xlim([0 5000])
% 
%     end
%     close all

    bestDelay = 6;
    if i == 3 || i == 6
        bestDelay = 1;
    end

    figure
    for l = 1:size(e3,2)
        aux = find(e3{bestDelay},1);

        plot(10*log10((e3{bestDelay}(aux:end))))
        hold on;
    end

    xlabel('Iterations','interpreter','latex');
    ylabel('MSE (dB)','interpreter','latex');
    
    ylim([-25 5]);
    
    xlim([0 8000]);

%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
    %     xlim([0 5000])
    
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-01' filesep 'mse' num2str(i)],'en' , figProp );
    if ismember(i,meanCountVector)
       
%         ct{auxIndex} = meanCount;
        
        if i < 9
        
        
            meanCountTran(auxIndex) = mean(meanCount{bestDelay}(find(meanCount{bestDelay},1):2000));
            meanCountSS(auxIndex) = mean(meanCount{bestDelay}(2001:4000));


            meanCountTran_2(auxIndex) = mean(meanCount{bestDelay}(4001:4001+2000));
            meanCountSS_2(auxIndex) = mean(meanCount{bestDelay}(4001+2000+1:end));

        else
            meanCountTran(auxIndex) = mean(meanCount{bestDelay}(find(meanCount{bestDelay},1):200));
            meanCountSS(auxIndex) = mean(meanCount{bestDelay}(201:4000));


            meanCountTran_2(auxIndex) = mean(meanCount{bestDelay}(4001:4001+200));
            meanCountSS_2(auxIndex) = mean(meanCount{bestDelay}(4001+200+1:end));
        end

        
         auxIndex = auxIndex + 1;
        
    end

end









rmpath(['.' filesep 'results']);
rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep '..' filesep]);
rmpath(['..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);

