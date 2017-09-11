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

% fileVector = {'teste6' 'teste5' 'teste4' 'teste' 'teste2' 'teste3'};

% fileVector = {'teste6' 'teste5' 'teste4' 'teste' 'teste2'};

fileVector = {'testeWam' 'teste4'};

fileVector = {'testeWam'};


for l = 1:length(fileVector)
    
    load([fileVector{l} '.mat']);
%     convergenceSample2 = zeros(size(e4,1),size(e4,2));
    for i = 1:size(e3,1)
        for j = 1:1%size(e4,2)
            
                x = e3{i};
%                 count = meanCount2{i,j};
                
%                 count{1} = zeros(10000,1);
%             for k = 1:size(x,1)
%                 figure
                aux = find(x,1);
                xAux = 10*log10(x(aux:4999));
                plot(xAux)
                
                hold on
%                 xAux2 = flipud(10*log10(x{1}(aux:4999)));
%                 plot(10*log10((x{1}(aux:end))))
%                 hold on
%                 y(i,j) = mean(10*log10((x{1}(4999 - 999:4999))));
%                 stdMse(i,j,l) = std(10*log10((x{1}(4999 - 999:4999))));
%                 convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
% %                 convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)- 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); %nao funcionou por causa de undershoot
%                 
%                 convergenceSample2(i,j,l) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)+ 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
            
        end
    end
%   convergenceSample2(l) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)+ 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
  updatesLin(l) = mean(meanCountLin{1}(aux:5000))*100;
  updatesNonLin(l) = mean(meanCountNonLin{1}(aux:5000))*100;
  updatesTotal(l) = mean(meanCountTotal{1}(aux:5000))*100;
end

H = legend('$\bar{\gamma} = 4\gamma$','$\bar{\gamma} = 2.5\gamma$','$\bar{\gamma} = 2\gamma$','$\bar{\gamma} = 1.5\gamma$','$\bar{\gamma} = 1\gamma$','$\bar{\gamma} = 0.5\gamma$');
set(H,'interpreter','latex')
% ylim([-15 10]);

ylabel('MSE [dB]','interpreter','latex');
xlabel('Iterations [$k$]','interpreter','latex');


















% 
fileVector = [39 43 51 47];


meanCountVector = 5:12;

auxIndex = 1;



figure

for l = 1:length(fileVector)
    
    load(['results' num2str(fileVector(l)) '.mat']);
%     convergenceSample2 = zeros(size(e4,1),size(e4,2));
    for i = 5:size(e4,1)
        for j = 1:1%size(e4,2)
            
            if fileVector(l) == 37 || fileVector(l) == 41 || fileVector(l) == 45
                x = e4{i};
            elseif fileVector(l) == 43 || fileVector(l) == 47
                x = e4{i,j};
                count = meanCount2{i,j};
            else
                x = e4{i,j};
                count{1} = zeros(10000,1);
            end
%             for k = 1:size(x,1)
%                 figure
                aux = find(x{1},1);
                xAux = 10*log10(x{1}(aux:4999));
                xAux2 = flipud(10*log10(x{1}(aux:4999)));
                plot(10*log10((x{1}(aux:end))))
                hold on
                y(i,j) = mean(10*log10((x{1}(4999 - 999:4999))));
                stdMse(i,j,l) = std(10*log10((x{1}(4999 - 999:4999))));
                convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
%                 convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)- 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); %nao funcionou por causa de undershoot
                
%                 convergenceSample2(i,j,l) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)+ 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
            
        end
    end
  convergenceSample2(l) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)+ 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
  updates(l) = mean(count{1}(aux:4999))*100;
end

H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');
set(H,'interpreter','latex')
ylim([-15 10]);

ylabel('MSE [dB]','interpreter','latex');
xlabel('Iterations','interpreter','latex');
% % if i > 3
%     set(gca,'ytick',[]);
% % end

formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'mse_DFE'],'en' , figProp );




fileVector = [42 43 44 46 47 48];


figProp = struct( 'size' , 32 , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

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
%                 figure
                xAux = 10*log10(x{1}(aux:4999));
                xAux2 = flipud(10*log10(x{1}(aux:4999)));
%                 plot(10*log10((x{1}(aux:end))))
                
                y(i,j) = mean(10*log10((x{1}(4999 - 999:4999))));
                stdMse2(i,j) = std(10*log10((x{1}(4999 - 999:4999))));
%                 convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
%                 convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)- 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); %nao funcionou por causa de undershoot
                
                convergenceSample2(i,j) = find(10*log10(x{1}(aux+1:4999)) < y(i,j)+ 2*stdMse2(i,j)*sign(2*stdMse2(i,j)),1,'first'); % funciona com  algumas ressalvas
                
%                 convergenceSample2(i,j) = find(xAux2 > y(i,j) +3.3*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first'); % funciona com  algumas ressalvas
                
                
                
%                 [~,minIdx] = min(xAux);
%                 while minIdx < convergenceSample2(i,j)
%                     convergenceSample2(i,j) = find(10*log10(x{1}(minIdx:4999)) < y(i,j)- 2*stdMse(i,j)*sign(2*stdMse(i,j)),1,'first');
%                 end
                    
                
%                 treta = diff(convergenceSample2Aux);
%                     treta = treta(treta>0);
                
%                 idx = 0;
%                 index = 1;
%                 
%                 while ~idx % nao funciona
%                     x = median(xAux((i:100+index-1)));
%                     index = index+1;
%                     if x < y(i,j) || index == 4999 - 1000;
%                         idx = 1;
%                     end
%                 end

               
                    
                

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



% fileVector = 42;    



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
%                 figure
%                 plot(10*log10(x{1}(aux:end)));
                stdMse2(i,j) = std(10*log10((x{1}(4999 - 999:4999))));
                convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
                convergenceSample2(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j) + 2*stdMse2(i,j)*sign(2*stdMse2(i,j)),1,'first');
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
    figure
    colormap(jet)
    imagesc(0:4,0:4,(y))
    
    if ~mod(l,3)
        c = colorbar;
        ylabel(c,'[dB]','interpreter','latex')
        set(c,'ylim',[-15 10]);
        set(c,'ytick',-15:5:10);
        caxis manual
        caxis([floor(min(minAxis)) ceil(max(maxAxis))]);
    end

    if l == 1 || l == 4
        ylabel('$M_{\mathrm{FF}}$','interpreter','latex');
        set(gca,'ytick',0:4);
    else
        set(gca,'ytick',[]);


    end
    
    set(gca,'xtick',0:4);

    xlabel('$M_{\mathrm{FB}}$','interpreter','latex');
    set(gca,'YDir','normal')


    ylim([0 4])
    formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'mseColor'  num2str(fileVector(l))],'en' , figProp );
% 


%     figure
%     colormap(jet)
%     imagesc(0:4,0:4,convergenceSample2)
%     c = colorbar;
%     ylabel(c,'Iterations until Convergence','interpreter','latex')
% %     set(c,'ylim',[-25 10]);
% %     set(c,'ytick',-25:5:10);
%     caxis manual
%     caxis([floor(min(minAxisConv)) ceil(max(maxAxisConv))]);
% 
%     xlabel('$M_{\mathrm{FB}}$','interpreter','latex');
%     ylabel('$M_{\mathrm{FF}}$','interpreter','latex');
%     set(gca,'ytick',0:4);

%     ylim([0 4])
%     formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-07' filesep 'conv'  num2str(fileVector(l))],'en' , figProp );

    
end



close all
fileVector = [37 41 49 45];




for i = 4:5
    figure

    for l = 1:length(fileVector)

        load(['results' num2str(fileVector(l)) '.mat']);

        x = e4{i};
%             for k = 1:size(x,1)
        aux = find(x{1},1);

        plot(10*log10((x{1}(aux:end))))

        mse(i,l) = mean(10*log10((x{1}(4999 - 999:4999))));
        stdMse(i,l) = std(10*log10((x{1}(4999 - 999:4999))));
        convergenceSample2(i,l) = find(10*log10(x{1}(aux:4999)) < mse(i,l)+2*stdMse(i,l)*sign(2*stdMse(i,l)),1,'first');
                 
        hold on

        xlim([0 10000]);
% %                 
%                  title([num2str(l) 'N_{FF} = ' num2str(i) ', N_{FB} = ' num2str(j)])
%             end
%             close all;
    end
%     ylabel('MSE [dB]','interpreter','latex');
    xlabel('Iterations','interpreter','latex');
    H = legend('V-PNLMS','VSM-PNLMS','V-RLS','VM-BEACON');

    set(H,'interpreter','latex')
    ylim([-15 10]);
    
    if i > 4
        set(gca,'ytick',[]);
    else
       ylabel('MSE [dB]','interpreter','latex');
    end
    
    formatFig( gcf ,['.' filesep 'figs' filesep '2017-07-12' filesep 'mse_FF_Eq' num2str(i)],'en' , figProp );
end


close all;



fileVector = [41 45];


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
                
            updates(i,j,l) = mean(x{1}(aux:4999))*100;
%             convergenceSample(i,j) = find(10*log10(x{1}(aux:4999)) < y(i,j),1,'first');
%                 
        end
    end

%     minAxis2(l) = min(min(updates));
%     
%     maxAxis2(l) = max(max(updates(updates~=0)));
    
 
    
    
    
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

