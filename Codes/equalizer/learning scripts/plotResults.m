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

fileVector = [2 4 5 6 8 9 12];

for i = 1:length(fileVector)
    load(['results0' num2str(fileVector(i)) '.mat']);

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
    %     xlim([0 5000])

    end
    close all

    bestDelay = 6;

    figure
    for l = 1:size(e3,2)
        aux = find(e3{bestDelay},1);

        plot(10*log10((e3{bestDelay}(aux:end))))
        hold on;
    end

    xlabel('Iterations','interpreter','latex');
    ylabel('MSE (dB)','interpreter','latex');

    H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
    set(H,'interpreter','latex')
    %     xlim([0 5000])


end

close all

load results05.mat;

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
%     xlim([0 5000])
    
end

bestDelay = 6;


figure
for l = 1:size(e3,2)
    aux = find(e3{bestDelay},1);

    plot(10*log10((e3{bestDelay}(aux:end))))
    hold on;
end

xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
set(H,'interpreter','latex')



close all

load results06.mat;

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
%     xlim([0 5000])
    
end

bestDelay = 6;


figure
for l = 1:size(e3,2)
    aux = find(e3{bestDelay},1);

    plot(10*log10((e3{bestDelay}(aux:end))))
    hold on;
end

xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
set(H,'interpreter','latex')












rmpath(['.' filesep 'results']);
rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep '..' filesep]);

