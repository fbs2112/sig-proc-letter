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


load testSML_MD.mat;
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


load testSML_Rec_Reg.mat;
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

load testSML_Rec_Reg2.mat;
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


