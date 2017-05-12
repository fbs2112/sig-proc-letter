clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep]);


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  


M = 10;

% load LMS.mat;
% 
% plot(10*log10((e3(M:end,1))))
% hold on;
% 
% 
% load NLMS.mat;
% 
% plot(10*log10((e3(M:end,1))))
% 
% load APA.mat;
% 
% plot(10*log10((e3(M:end,1))))
% 
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('LMS','NLMS','APA');
% 
% set(H,'interpreter','latex')
% 
% % formatFig( gcf ,['.' filesep 'figs 2017-04-27' filesep 'mseCase1'],'en' , figProp );
% 
% 
% 
% figure
% load NLMS_2.mat;
% 
% plot(10*log10((e3(M:end,1))))
% hold on
% load APA_2.mat;
% 
% plot(10*log10((e3(M:end,1))))
% 
% 
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% H = legend('NLMS','APA');
% 
% set(H,'interpreter','latex')

% formatFig( gcf ,['.' filesep 'figs 2017-04-27' filesep 'mseCase2'],'en' , figProp );






% 
% load results01.mat;
% figure
% 
% for l = 1:size(e3,2)
% % figure
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSML1'],'en' , figProp );
% 
% 
% 
% % load testInit.mat;
% % figure
% % 
% % for l = 1:size(e3,2)
% % % figure
% %     plot(10*log10((e3(M:end,l))))
% %     hold on;
% % end
% % xlabel('Iterations','interpreter','latex');
% % ylabel('MSE (dB)','interpreter','latex');
% % 
% % H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% % set(H,'interpreter','latex')
% 
% % formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSML1'],'en' , figProp );
% 
% 
% 
% load results02.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% 
% % formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseVolterra1'],'en' , figProp );
% 
% 
% 
% load results03.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% 
% formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSML2'],'en' , figProp );
% 
% 
% 
% load results04.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% 
% load results05.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% 
% formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSMSML1'],'en' , figProp );
% 
% load results06.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSMVolterra1'],'en' , figProp );
% 
% 
% 
% load results07.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSMSML2'],'en' , figProp );
% 
% load results08.mat;
% figure
% 
% for l = 1:size(e3,2)
%     plot(10*log10((e3(M:end,l))))
%     hold on;
% end
% 
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
% 
% 
% formatFig( gcf ,['.' filesep 'figs 2017-05-04' filesep 'mseSMVolterra2'],'en' , figProp );








load testSML.mat;
figure

for l = 1:size(e3,2)
    plot(10*log10((e3{l}(M:end))))
    hold on;
end

xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
set(H,'interpreter','latex')
xlim([0 5000])

load testSML2.mat;
figure

for l = 1:size(e3,2)
    plot(10*log10((e3{l}(M:end))))
    hold on;
end

xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
set(H,'interpreter','latex')
xlim([0 5000])


rmpath(['.' filesep 'results']);


