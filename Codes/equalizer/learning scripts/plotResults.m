clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep 'simParameters' filesep]);
addpath(['..' filesep '..' filesep]);


load paramEq.mat;


numberOfSymbols = 2^numberOfBits;

delayVector = 1:N+length(h);


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  


load testPAPA_DFE_Volterra.mat;

for delay = 1:size(e3,1)
    figure
    for l = 1:size(e3,2)
        aux = find(e3{delay},1);
        
        plot(10*log10((e3{delay}(aux:end))))
        hold on;
    end

    xlabel('Iterations','interpreter','latex');
    ylabel('MSE (dB)','interpreter','latex');

    H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
    set(H,'interpreter','latex')
    xlim([0 5000])
    
end


% 
% close all;
% 
% M = 17;
% 
% load testOBE_DFE_FF_FB.mat;
% 
% for delay = 1:length(delayVector)
%     figure
%     for l = 1:size(e3,2)
%         plot(10*log10((e3{delay}(M:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
% 
%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
%     xlim([0 5000])
%     
% end
% 
% close all;
% 
% M = 17;
% 
% delayVector = 1;
% 
% load testOBE_DFE.mat;
% 
% for delay = 1:length(delayVector)
%     figure
%     for l = 1:size(e3,2)
%         plot(10*log10((e3{delay}(M:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
% 
%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
%     xlim([0 5000])
%     
% end
% 
% 
% 
% 
% close all;
% 
% M = 17;
% 
% 
% load testSM_PAPA_LinEq.mat;
% 
% for delay = 1:size(e3,1)
%     figure
%     for l = 1:size(e3,2)
%         plot(10*log10((e3{delay}(M:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
% 
%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
%     xlim([0 500])
%     
% end
% 
% close all
% 
% 
% load testSM_PAPA_DFE_LinEq.mat;
% 
% for delay = 1:size(e3,1)
%     figure
%     for l = 1:size(e3,2)
%         plot(10*log10((e3{delay}(M:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
% 
%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
%     xlim([0 500])
%     
% end
% 
% 
% close all;
% 
% load testSM_PAPA_DFE_Volterra.mat;
% 
% for delay = 1:size(e3,1)
%     figure
%     for l = 1:size(e3,2)
%         plot(10*log10((e3{delay}(M:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
% 
%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
%     xlim([0 500])
%     
% end
% 
% 
% 
% close all;
% 
% load testPAPA_DFE_Volterra.mat;
% 
% for delay = 1:size(e3,1)
%     figure
%     for l = 1:size(e3,2)
%         plot(10*log10((e3{delay}(M:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
% 
%     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
%     set(H,'interpreter','latex')
% %     xlim([0 500])
%     
% end
% 







rmpath(['.' filesep 'results']);
rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep '..' filesep]);

