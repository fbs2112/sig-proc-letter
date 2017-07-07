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

fileVector = [38 39 40 42 43 44 46 47 48];

% fileVector = [33];


meanCountVector = 5:12;

auxIndex = 1;



for l = 3:length(fileVector)
    
    load(['results' num2str(fileVector(l)) '.mat']);
    convergenceSample2 = zeros(size(e4,1),size(e4,2));
    for i = 1:size(e4,1)
        for j = 1:size(e4,2)
            
            if fileVector(l) == 37 || fileVector(l) == 41 || fileVector(l) == 45
                x = e4{i};
            else
                x = e4{i,j};
            end
%             for k = 1:size(x,1)
%                 figure
                aux = find(x{1},1);
                figure
                xAux = 10*log10(x{1}(aux:4999));
                xAux2 = flipud(10*log10(x{1}(aux:4999)));
                plot(10*log10((x{1}(aux:end))))
                
                y(i,j) = mean(10*log10((x{1}(4999 - 999:4999))));
                stdMse(i,j) = std(10*log10((x{1}(4999 - 999:4999))));
                convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
%                 convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)- 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); %nao funcionou por causa de undershoot
                
                convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)+ 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
                
                convergenceSample2(i,j) = find(xAux2 > y(i,j) +3.3*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
                
                
                
%                 [~,minIdx] = min(xAux);
%                 while minIdx < convergenceSample2(i,j)
%                     convergenceSample2(i,j) = find(10*log10(x{1}(minIdx:4999)) < y(i,j)- 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first');
%                 end
                    
                
                treta = diff(convergenceSample2Aux);
                    treta = treta(treta>0);
                
                idx = 0;
                index = 1;
                
                while ~idx % nao funciona
                    x = median(xAux((i:100+index-1)));
                    index = index+1;
                    if x < y(i,j) || index == 4999 - 1000;
                        idx = 1;
                    end
                end

               
                    
                

        end
    end
%     figure
%     colormap(jet)
%     imagesc(y)
%     c = colorbar;
%     ylabel(c,'[dB]','interpreter','latex')
%     set(c,'ylim',[-25 10]);
%     set(c,'ytick',-25:5:10);
% %     caxis manual
% %     caxis([bottom top]);
% 
%     xlabel('$N_{\mathrm{FB}}$','interpreter','latex');
%     ylabel('$N_{\mathrm{FF}}$','interpreter','latex');
%     set(gca,'ytick',1:5);
% 
%     ylim([1 5])
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-07' filesep 'mse'  num2str(fileVector(l))],'en' , figProp );

    minAxis(l) = min(min(y));
    
    maxAxis(l) = max(max(y(y~=0)));
    
    
    minAxisConv(l) = min(min(convergenceSample2));
    
    maxAxisConv(l) = max(max(convergenceSample2(convergenceSample2~=0)));
    
    
    
end







for l = 1:length(fileVector)
    
    load(['results' num2str(fileVector(l)) '.mat']);
    for i = 1:size(e4,1)
        for j = 1:size(e4,2)
            
            if fileVector(l) == 37 || fileVector(l) == 41 || fileVector(l) == 45
                x = e4{i};
            else
                x = e4{i,j};
            end
%             for k = 1:size(x,1)
%                 figure
                aux = find(x{1},1);
% 
%                 plot(10*log10((x{1}(aux:end))))
                
                y(i,j) = mean(10*log10((x{1}(4999 - 999:4999))));
                stdMse(i,j) = std(10*log10((x{1}(4999 - 999:4999))));
                convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
                convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j) + 2*stdMse(i,j),1,'first');
%                 treta(i,j,:) = diff(10*log10(x{1}(aux:4999)));
%                 convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) > y(i,j) + 0.2*abs(y(i,j)),1,'last');
%                  H = legend('$N = 1$','$N = 2$','$N = 3$','$N = 4$','$N = 5$');
%                 set(H,'interpreter','latex')
%                 ylim([-25 20]);
%                 hold on
    
%                 xlim([0 10000]);
% %                 formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-09' filesep 'msePAPA_Volterra'],'en' , figProp );
%                  title([num2str(l) 'N_{FF} = ' num2str(i) ', N_{FB} = ' num2str(j)])
%             end
%             close all;
        end
    end
%     figure
%     colormap(jet)
%     imagesc(y)
%     c = colorbar;
%     ylabel(c,'[dB]','interpreter','latex')
%     set(c,'ylim',[-25 10]);
%     set(c,'ytick',-25:5:10);
%     caxis manual
%     caxis([floor(min(minAxis)) ceil(max(maxAxis))]);
% 
%     xlabel('$N_{\mathrm{FB}}$','interpreter','latex');
%     ylabel('$N_{\mathrm{FF}}$','interpreter','latex');
%     set(gca,'ytick',1:5);
% 
%     ylim([1 5])
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-07' filesep 'mse'  num2str(fileVector(l))],'en' , figProp );



    figure
    colormap(jet)
    imagesc(convergenceSample2)
    c = colorbar;
    ylabel(c,'Iterations until Convergence','interpreter','latex')
%     set(c,'ylim',[-25 10]);
%     set(c,'ytick',-25:5:10);
    caxis manual
    caxis([floor(min(minAxisConv)) ceil(max(maxAxisConv))]);

    xlabel('$N_{\mathrm{FB}}$','interpreter','latex');
    ylabel('$N_{\mathrm{FF}}$','interpreter','latex');
    set(gca,'ytick',1:5);

    ylim([1 5])
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-07' filesep 'conv'  num2str(fileVector(l))],'en' , figProp );

    
end



close all
fileVector = [37 41 45];




for i = 1:size(e4,1)
    figure

    for l = 1:length(fileVector)

        load(['results' num2str(fileVector(l)) '.mat']);

        x = e4{i};
%             for k = 1:size(x,1)
        aux = find(x{1},1);

        plot(10*log10((x{1}(aux:end))))

        mse(i,l) = mean(10*log10((x{1}(4999 - 999:4999))));
        convergenceSample2(i,l) = find(10*log10(x{1}(aux:4999)) < mse(i,l),1,'first');
                 
        hold on

        xlim([0 10000]);
% %                 
%                  title([num2str(l) 'N_{FF} = ' num2str(i) ', N_{FB} = ' num2str(j)])
%             end
%             close all;
    end
    
    H = legend('PAPA Volterra','SM-PAPA Volterra','Modified BEACON');
    set(H,'interpreter','latex')
    ylim([-20 10]);
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-07' filesep num2str(i)],'en' , figProp );
    
end


close all;


fileVector = [41 42 43 44 45 46 47 48];


for l = 1:length(fileVector)
    
    load(['results' num2str(fileVector(l)) '.mat']);
    for i = 1:size(e4,1)
        for j = 1:size(e4,2)
            
            if fileVector(l) == 41 || fileVector(l) == 45
                x = meanCount2{i};
            else
                x = meanCount2{i,j};
            end
            aux = find(x{1},1);
                
            updates(i,j) = mean(x{1}(4999 - 999:4999))*100;
%             convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
%                 
        end
    end

    minAxis2(l) = min(min(updates));
    
    maxAxis2(l) = max(max(updates(updates~=0)));
    
 
    
    
    
end







fileVector = [41 42 43 44 45 46 47 48];


for l = 1:length(fileVector)
    
    load(['results' num2str(fileVector(l)) '.mat']);
    for i = 1:size(e4,1)
        for j = 1:size(e4,2)
            
            if fileVector(l) == 41 || fileVector(l) == 45
                x = meanCount2{i};
            else
                x = meanCount2{i,j};
            end
            aux = find(x{1},1);
                
            updates(i,j) = mean(x{1}(4999 - 999:4999))*100;
%             convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
%                 
        end
    end

    figure
    colormap(jet)
    imagesc(updates)
    c = colorbar;
    ylabel(c,'Rate of Updates','interpreter','latex')
    set(c,'ylim',[0 100]);
    set(c,'ytick',0:10:100);
    caxis manual
    caxis([min(minAxis2) max(maxAxis2)]);

    xlabel('$N_{\mathrm{FB}}$','interpreter','latex');
    ylabel('$N_{\mathrm{FF}}$','interpreter','latex');
    set(gca,'ytick',1:5);

    ylim([1 5])
    formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-07' filesep 'update'  num2str(fileVector(l))],'en' , figProp );
    
    
    
end



















% for l = 1:length(fileVector)
%     
%     load(['results' num2str(fileVector(l)) '.mat']);
%     for i = 1:size(e4,1)
%         for j = 1:size(e4,2)
%             x = e4{i,j};
% %             for k = 1:size(x,1)
% %                 figure
%                 aux = find(x{i+1},1);
% 
%                 plot(10*log10((x{i+1}(aux:end))))
%                 H = legend('$N = 1$','$N = 2$','$N = 3$','$N = 4$','$N = 5$');
%                 set(H,'interpreter','latex')
%                 xlabel('Iterations','interpreter','latex');
%                 ylabel('MSE (dB)','interpreter','latex');
%                 ylim([-25 20]);    
%                 xlim([0 10000]);
%                 hold on
% 
% %                  title([num2str(l) 'N_{FF} = ' num2str(i) ', N_{FB} = ' num2str(j)])
% %             end
% %             close all;
%         end
%     end
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-09' filesep 'mseSMPAPA_Volterra'],'en' , figProp );

% end












% load testLen2.mat
% 
% 
% 
% for i = 1:size(e4,1)
%     
%     x = e4{i,1};
%     for k = 1:size(x,1)
%         figure
%         aux = find(x{k},1);
% 
%         plot(10*log10((x{k}(aux:end))))
%         title(['N = ' num2str(i)])
% 
%     end
%     close all;
%    
% end
%             
% 
% 
% 
% 
% 
% 
% 
% for i = 1:length(fileVector)
%     
%     if i > 9
%         load(['results' num2str(fileVector(i)) '.mat']);
%     else
%         load(['results0' num2str(fileVector(i)) '.mat']);
%     end
%     
% %     load teste.mat
% %     for delay = 1:size(e3,1)
% %         figure
% %         for l = 1:size(e3,2)
% %             aux = find(e3{delay},1);
% % 
% %             plot(10*log10((e3{delay}(aux:end))))
% %             hold on;
% %         end
% % 
% %         xlabel('Iterations','interpreter','latex');
% %         ylabel('MSE (dB)','interpreter','latex');
% % 
% %         H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% %         set(H,'interpreter','latex')
% %     %     xlim([0 5000])
% % 
% %     end
% %     close all
% 
%     bestDelay = 6;
%     if i == 3 || i == 6
%         bestDelay = 1;
%     end
% 
%     figure
%     for l = 1:size(e3,2)
%         aux = find(e3{bestDelay},1);
% 
%         plot(10*log10((e3{bestDelay}(aux:end))))
%         hold on;
%     end
% 
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
%     
%     ylim([-25 5]);
%     
%     xlim([0 8000]);
% 
% %     H = legend('$L = 0$','$L = 1$','$L = 2$','$L = 3$');
% %     set(H,'interpreter','latex')
%     %     xlim([0 5000])
%     
% %     formatFig( gcf ,['.' filesep 'figs' filesep '2017-06-01' filesep 'mse' num2str(i)],'en' , figProp );
%     if ismember(i,meanCountVector)
%        
% %         ct{auxIndex} = meanCount;
%         
%         if i < 9
%         
%         
%             meanCountTran(auxIndex) = mean(meanCount{bestDelay}(find(meanCount{bestDelay},1):2000));
%             meanCountSS(auxIndex) = mean(meanCount{bestDelay}(2001:4000));
% 
% 
%             meanCountTran_2(auxIndex) = mean(meanCount{bestDelay}(4001:4001+2000));
%             meanCountSS_2(auxIndex) = mean(meanCount{bestDelay}(4001+2000+1:end));
% 
%         else
%             meanCountTran(auxIndex) = mean(meanCount{bestDelay}(find(meanCount{bestDelay},1):200));
%             meanCountSS(auxIndex) = mean(meanCount{bestDelay}(201:4000));
% 
% 
%             meanCountTran_2(auxIndex) = mean(meanCount{bestDelay}(4001:4001+200));
%             meanCountSS_2(auxIndex) = mean(meanCount{bestDelay}(4001+200+1:end));
%         end
% 
%         
%          auxIndex = auxIndex + 1;
%         
%     end
% 
% end
% 
% 







rmpath(['.' filesep 'results']);
rmpath(['..' filesep 'simParameters' filesep]);
rmpath(['..' filesep '..' filesep]);
rmpath(['..' filesep '..' filesep 'sysId' filesep 'Utils' filesep]);

