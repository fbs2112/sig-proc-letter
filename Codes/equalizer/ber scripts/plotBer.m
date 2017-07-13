clear;
clc;
close all;


addpath(['.' filesep 'results']);
addpath(['..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , 24 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);


fileVector = [1 5 9];



for i = 3:5
    
    figure
    
    for j = 1:length(fileVector)

        load (['resultsBER0' num2str(fileVector(j)) '.mat']);

        semilogy(SNR,squeeze(ber(i,:)))
        hold on;
    end
    H = legend('V-PNLMS','VSM-PNLMS','VM-BEACON');
    set(H,'location','SouthWest');
    set(H,'interpreter','latex')
    
    xlabel('SNR [dB]','interpreter','latex');
    ylim([1e-6 1e0])
    xlim([0 30]);
    
    set(gca,'xtick',SNR);
  
    if i == 3
        ylabel('BER','interpreter','latex');
    else
        set(gca,'ytick',[]);
    end
        
    
    formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'berFF'  num2str(i)],'en' , figProp );
    
end

