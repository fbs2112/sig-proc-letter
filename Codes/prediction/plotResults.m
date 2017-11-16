clear;
clc;
close all;




addpath(['.' filesep 'results']);
addpath(['.' filesep 'data' filesep]);

dataAux = xlsread('daily-maximum-temperatures-in-me.xls');
data = dataAux(15:end,2);
globalLength = 400;

data = buffer(data,globalLength);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  


figPath = ['.' filesep 'figs' filesep];

fileNames = {'resultsPAPA','resultsSMPAPA','resultsRLS','resultsBEACON'};


for i = 1:length(fileNames)
    
    load([fileNames{i} '.mat']);
    plot(10*log10(e3{1,1}));
    hold on;
end


H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
set(H,'interpreter','latex')
xlim([1 400]);

ylabel('MSE [dB]','interpreter','latex');
xlabel('Iterations','interpreter','latex');
formatFig( gcf ,[figPath 'mse'],'en' , figProp );

figure

for i = 1:length(fileNames)
    
    load([fileNames{i} '.mat']);
    
    plot(10*log10(e3Plot{1,1}));
    hold on;
end


H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
set(H,'interpreter','latex')
xlim([1 400]);

ylabel('MRE [dB]','interpreter','latex');
xlabel('Iterations','interpreter','latex');
formatFig( gcf ,[figPath 'mseN'],'en' , figProp );



index = 8;
for i = 1:length(fileNames)
    figure
    load([fileNames{i} '.mat']);
    plot(200:400,data(200:400,index-1),'*-');
    hold on
    plot(200:400,y(200:400,index-1),'r*-');
    
    
    H = legend('Real signal','Predicted signal');
    set(H,'interpreter','latex')
    ylim([10 40]);
    
    ylabel('Temperature [$^\circ$C]','interpreter','latex');
    xlabel('Samples','interpreter','latex');
    
    if i > 2
        updates(i-2) = meanCount{1,1};
    end
    formatFig( gcf ,[figPath 'track' num2str(i)],'en' , figProp );

    
end



figure


plot(dataAux(15:end,2))

xlim([0 length(dataAux(15:end,2))]);

ylabel('Temperature [$^\circ$C]','interpreter','latex');
xlabel('Samples','interpreter','latex');

formatFig( gcf ,[figPath 'data'],'en' , figProp );


