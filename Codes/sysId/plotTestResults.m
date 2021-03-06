clear;
clc;
close all;



addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep]);
addpath(['.' filesep 'Utils' filesep]);


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  



load results01.mat;
figure

for l = 1:size(e3,2)
    aux = find(e3{l},1);
    plot(10*log10((e3{l}(aux:end))))
    hold on;
    meanMSE(1) = mean(10*log10((e3{l}(999 - 499:999))));
    stdMse(1) = std(10*log10((e3{l}(999 - 499:999))));
    
end
convergenceSample(1) = find(10*log10(e3{l}(aux:999)) < meanMSE(1) + 2*stdMse(1)*sign(2*stdMse(1)),1,'first');




load results02.mat;

for l = 1:size(e3,2)
    aux = find(e3{l},1);
    plot(10*log10((e3{l}(aux:end))))
    hold on;
    meanMSE(2) = mean(10*log10((e3{l}(999 - 499:999))));
    stdMse(2) = std(10*log10((e3{l}(999 - 499:999))));
end
convergenceSample(2) = find(10*log10(e3{l}(aux:999)) < meanMSE(2) + 2*stdMse(2)*sign(2*stdMse(2)),1,'first');


load results04.mat;

for l = 1:size(e3,2)
    aux = find(e3(:,l),1);
    plot(10*log10((e3(aux:end,l))))
    hold on;
    
    meanMSE(3) = mean(10*log10((e3(999 - 499:999))));
    stdMse(3) = std(10*log10((e3(999 - 499:999))));
end

convergenceSample(3) = find(10*log10(e3(aux:999)) < meanMSE(3) + 2*stdMse(3)*sign(2*stdMse(3)),1,'first');

load results03.mat;

for l = 1:size(e3,2)
    aux = find(e3(:,l),1);
    plot(10*log10((e3(aux:end,l))))
    hold on;
    
    meanMSE(4) = mean(10*log10((e3(999 - 499:999))));
    stdMse(4) = std(10*log10((e3(999 - 499:999))));
end

convergenceSample(4) = find(10*log10(e3(aux:999)) < meanMSE(4) + 2*stdMse(4)*sign(2*stdMse(4)),1,'first');



xlabel('Iterations [$k$]','interpreter','latex');
ylabel('MSE [dB]','interpreter','latex');

H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
set(H,'interpreter','latex')
xlim([0 2000])
ylim([-35 15])


formatFig( gcf ,['.' filesep 'figs' filesep 'mseEq'],'en' , figProp );


figure

load results01.mat;

for l = 1:size(misalignment,2)
    aux = find(misalignment{l},1);
    plot(10*log10((misalignment{l}(aux:end))))
    hold on;
end

load results02.mat;

for l = 1:size(misalignment,2)
    aux = find(misalignment{l},1);
    plot(10*log10((misalignment{l}(aux:end))))
    hold on;
end


load results04.mat;

for l = 1:size(misalignment,2)
    aux = find(misalignment,1);
    plot(10*log10((misalignment(aux:end))))
    hold on;
end


load results03.mat;

for l = 1:size(misalignment,2)
    aux = find(misalignment,1);
    plot(10*log10((misalignment(aux:end))))
    hold on;
end

xlabel('Iterations [$k$]','interpreter','latex');
ylabel('Misalignment [dB]','interpreter','latex');

H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
set(H,'interpreter','latex')
xlim([0 2000])
ylim([-60 10])

formatFig( gcf ,['.' filesep 'figs' filesep 'misEq'],'en' , figProp );


load results02.mat;
aux = find(misalignment{l},1);

meanCountSMPAPA_NLMS = mean(meanCount{1}(aux:999))*100;

load results03.mat;

aux = find(misalignment,1);
meanCountBEACON = mean(meanCount(aux:999))*100;


meanCountSMPAPA_NLMS_Tran = mean(meanCount{1}(find(meanCount{1},1):200));
meanCountSMPAPA_NLMS_SS = mean(meanCount{1}(201:1000));


meanCountSMPAPA_NLMS_Tran_2 = mean(meanCount{1}(1001:1001+200));
meanCountSMPAPA_NLMS_SS_2 = mean(meanCount{1}(1001+200+1:end));



figure

for l = 1:size(misalignment,2)
    aux = find(misalignment{l},1);
    plot(10*log10((misalignment{l}(aux:end))))
    hold on;
end

xlabel('Iterations','interpreter','latex');
ylabel('Misalignment (dB)','interpreter','latex');

% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
xlim([0 2000])
ylim([-60 10])

% formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-01' filesep 'misSMPAPA_NLMS'],'en' , figProp );




xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
xlim([0 2000])
ylim([-35 10])
% formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-01' filesep 'mseSM_OBE'],'en' , figProp );


meanCountSM_OBE_Tran = mean(meanCount(find(meanCount,1):200));
meanCountSM_OBE_SS = mean(meanCount(201:1000));

meanCountSM_OBE_Tran_2 = mean(meanCount(1001:1001+200));
meanCountSM_OBE_SS_2 = mean(meanCount(1001+200+1:end));


figure

for l = 1:size(misalignment,2)
    aux = find(misalignment(:,l),1);
    plot(10*log10((misalignment(aux:end,l))))
    hold on;
end

xlabel('Iterations','interpreter','latex');
ylabel('Misalignment (dB)','interpreter','latex');

% H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% set(H,'interpreter','latex')
xlim([0 2000])
ylim([-60 15])


% formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-01' filesep 'misSM_OBE'],'en' , figProp );


rmpath(['.' filesep 'results']);
rmpath(['..' filesep '..' filesep]);
rmpath(['.' filesep 'Utils' filesep]);


