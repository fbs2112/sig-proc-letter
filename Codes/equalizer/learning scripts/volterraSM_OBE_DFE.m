%%
%Volterra NLMS DFE



clear;
clc;
close all;

maxRuns = 3000;
maxIt = 1000;

% N = 8;
% 
% adapFiltLength = (N^2+N)/2 + N;
% 
% % adapFiltLength = N;
% 

% auxMatrix = triu(ones(N));
% [l1,l2] = find(auxMatrix);

% h = [0.34-(0.27*1i) 0.87+(0.43*1i) 0.34-(0.21*1i)]; 

% L = round((adapFiltLength + length(h))/2);

signalPower = 1;
noisePower = signalPower/db2pow(30);
gamma = 4*sqrt(5*noisePower);

% barGammaVector = 1:0.5:4;


SNRAux = db2pow(30);
% mu = 0.8;


volterraFFFlag = 0;
volterraFBFlag = 1;

% lambda = 1;
% varW = 0;

% h(:,adapFiltLength + round((N+4)/2) + 10) = [1 0 0 1].';
% delayVector = round((length(h) + adapFiltLength)/2);
% delayVector = 4:8;

h = [1 0.2 -0.3];




barGammaVector = 1;
feedforwardLength = 4;
feedbackLength = 2;

adaptfiltFF = (feedforwardLength^2+feedforwardLength)/2 + feedforwardLength;
adaptfiltFB = (feedbackLength^2+feedbackLength)/2 + feedbackLength;

adaptfilt = adaptfiltFF + adaptfiltFB;

auxMatrix = triu(ones(feedforwardLength));
[l1FF,l2FF] = find(auxMatrix);

auxMatrix = triu(ones(feedbackLength));
[l1FB,l2FB] = find(auxMatrix);




% delayVector = 3;
delayVector = round((length(h) + feedforwardLength)/2);

if ~volterraFFFlag
    adaptfiltFF = feedforwardLength;
end

if ~volterraFBFlag
    adaptfiltFB = feedbackLength;
end

    

adapFiltLength = adaptfiltFF + adaptfiltFB;

% delayVector = round((N+4)/2);

% for N = 8:8
%     adapFiltLength = (N^2+N)/2 + N;
% 
%     auxMatrix = triu(ones(N));
%     [l1,l2] = find(auxMatrix);




for barGammaIndex = 1:length(barGammaVector)
    count = zeros(maxIt,1);

%     count = zeros(maxIt,length(barGammaVector));
    for delay = 1:length(delayVector)

        for L = 0:0
            
            u = zeros(L+1,1);
            u(1) = 1;
            

            w2 = zeros(adapFiltLength,maxRuns,maxIt);
            for j = 1:maxIt
                j
                
                P = zeros(adapFiltLength,adapFiltLength,maxRuns*2);
                sigma = zeros(maxRuns*2,1);
                delta = zeros(maxRuns,1);
                G = zeros(maxRuns*2,1);
                lambda = zeros(maxRuns*2,1);

                P(:,:,adapFiltLength + delayVector(delay) + 10) = eye(adapFiltLength)*1e-6;
                sigma(adapFiltLength + delayVector(delay) + 10) = 1;
                

                input = randi([0,3],maxRuns*2,1);
                pilot = qammod(input,4);

                pilot = pilot.*sqrt(signalPower/var(pilot));

    %             xAux2 = filter(h,1,pilot); % linear case

    %             xAux2(1:2,1) = pilot(1:2).^2; % nonlinear case
    %             xAux2(3,1) = pilot(3).^2 + pilot(2)*pilot(1);
    % 
    %             for m = 4:length(pilot)
    %                xAux2(m,1) = pilot(m).^2 + pilot(m-1)*pilot(m-2) + pilot(m-3);
    %             end

    %              for m = 1:length(pilot)
    %                 xAux2(m,1) = pilot(m) + 0.2*(pilot(m)^2);
    %              end
    %                  xAux2 = filter(h,1,pilot);
    %             xAux2 = exp(-pilot);

                  xAux2 = filter(h,1,pilot);
                  
                  xAux2 = xAux2 + 0.2*(xAux2.^2) + 0.05*(xAux2.^3);
    %                  xAux2 = xAux2.^2;
    %               xAux2(1) = pilot(1);
    %               for m = 2:length(pilot)
    %                   xAux2(m,1) = pilot(m) + pilot(m-1);
    %               end


    %             xAux2(1,1) = pilot(1); % nonlinear case
    % %             xAux2(3,1) = pilot(3).^2 + pilot(2)*pilot(1);
    % 
    %             for m = 2:length(pilot)
    %                xAux2(m,1) = pilot(m) + pilot(m-1);
    %             end

        %         xAux2 = exp(pilot.*conj(pilot));


                n = randn(maxRuns*2,1) + randn(maxRuns*2,1)*1i;
                powerSignal = xAux2'*xAux2./(maxRuns*2);
                powerNoiseAux = n'*n/(maxRuns*2);
                powerNoise = (powerSignal/SNRAux);
                n = n.*sqrt(powerNoise/powerNoiseAux);

                xAux = xAux2 + n;
    % 
    % 
%                 xTDL = (buffer(pilot,N,N-1,'nodelay'));

                theta = zeros(adapFiltLength,maxRuns);

    %                 w = zeros(N,maxRuns);


                for k = (adapFiltLength + delayVector(delay) + L + 10):maxRuns
                    
                    x(:,k) = xAux(k:-1:k-feedforwardLength+1);
                    
                    yHat(:,k) = (pilot(-delayVector(delay) + k + 1 -1:-1:-delayVector(delay) + k + 1 - feedbackLength - 1 + 1));
                    
                    if volterraFFFlag
                    
                        aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                        for lIndex = 1:length(l1FF)
                            aux(lIndex,1) = x(l1FF(lIndex),k)*(x(l2FF(lIndex),k));
                        end
                        xConc = [x(:,k);aux];
                    else
                        xConc = x(:,k);
                    end
                    
                    
                    if volterraFBFlag
                        aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                        for lIndex = 1:length(l1FB)
                            aux(lIndex,1) = yHat(l1FB(lIndex),k)*(yHat(l2FB(lIndex),k));
                        end

                        yHatConc = [yHat(:,k);aux];
                    else
                        yHatConc = yHat(:,k);
                    end
                    
                    if ~volterraFFFlag && ~volterraFBFlag 
                        xConc = x(:,k);
                        yHatConc = yHat(:,k);
                    end
                    
                    z = [xConc;yHatConc];
                    
                    d(k) = (pilot(-delayVector(delay) + k + 1));
                    
                    delta(k) = d(k) - theta(:,k).'*z;
                    
                    G(k) = z.'*P(:,:,k)*conj(z);
                    
                    
                    
                    
                    if abs(delta(k)) <= gamma*barGammaVector(barGammaIndex)

                        lambda(k) = 0;
                        count(j) = count(j) + 1;
                        P(:,:,k+1) = P(:,:,k);
                        theta(:,k+1) = theta(:,k);
                        sigma(k+1) = sigma(k); 

                    else
                        lambda(k) = (1/G(k))*((abs(delta(k))/(gamma*barGammaVector(barGammaIndex))) - 1);
                         P(:,:,k+1) = (P(:,:,k)) - (lambda(k)*(P(:,:,k))*conj(z)*z.'*(P(:,:,k)))/(1+lambda(k)*G(k));

                        theta(:,k+1) = theta(:,k) + lambda(k)*P(:,:,k+1)*conj(z)*delta(k);

                        sigma(k+1) = sigma(k) - (lambda(k)*delta(k)^2)/(1+lambda(k)*G(k)) + lambda(k)*delta(k)^2;

                    end
                    
                    
                    
                    
                    
                    
%                     absoluteValueError = abs(e(k));
% % 
%                     if absoluteValueError > barGamma*barGammaVector(barGammaIndex)
%                         mu(k) = 1 - barGamma*barGammaVector(barGammaIndex)/absoluteValueError;
%                         w(:,k+1) = w(:,k) + mu(k)*z*((z'*z+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;
%                     else
%                         mu(k) = 0;
%                         count(j) = count(j)+1;
%                         w(:,k+1) = w(:,k);
% %                         if count(j) >= maxIt
% %                             pause;
% %                         end
%                         
%                     end
                    
                    
                    
                    
                    
                    
                    
%                     w(:,k+1) = w(:,k) + mu*z*((z'*z+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;

%                     for l = 0:L
%                         xAP(:,l+1) = xAux(k-l:-1:k-N+1-l);
%                     end

    %                     xAP = flipud(xAux);

%                     xTDLConc = zeros(adapFiltLength,L+1);
% 
%                     for l3 = 1:L+1
%                         xTDLAux = zeros((N*N+N)/2,1);
% 
% 
% 
%                         for lIndex = 1:length(l1)
%                             xTDLAux(lIndex,1) = xAP(l1(lIndex),l3)*(xAP(l2(lIndex),l3));
% 
%     %                                 xTDLAux(counterAux,1) = xTDL(l1,k+l3-1)*xTDL(l2,k+l3-1);
%                         end
% 
% 
%                         xTDLConc(:,l3) = [xAP(:,l3);xTDLAux];
%         %                 xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1);xTDLAux];
%                     end

    %                     xTDLConc = xAP;

    %                 xAP = xTDLConc(:,k:-1:k-L);


%                     d(k) = (pilot(-delayVector(delay) + k + 1)); 
        %             d(k) = (pilot(-L2 + adapFiltLength+1 - L)); 

%                     e(k) = d(k) - w(:,k)'*xTDLConc(:,1);
%     %                 if abs(e(k)) > 100
%     %                     abs(e(k));
%     %                 end
%                     absoluteValueError = abs(e(k));
% 
%                     if absoluteValueError > barGamma*barGammaVector(barGammaIndex)
%                         mu(k) = 1 - barGamma*barGammaVector(barGammaIndex)/absoluteValueError;
%     %                     mu(k) = 0.2;
%                         w(:,k+1) = w(:,k) + mu(k)*xTDLConc*((xTDLConc'*xTDLConc+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;
%     %                     b(k) = ((xTDLConc'*xTDLConc+gamma*eye(L+1)));
%                     else
%                         mu(k) = 0;
%                         count(j) = count(j)+1;
%                         w(:,k+1) = w(:,k);
%                         if count(j) >= maxIt
%                             pause;
%                         end
%                         
%                     end
    %                     nW = randn(length(h(k)),1);
    %                     nW = nW.*sqrt(varW/var(nW));
    %                     h(:,k+1) = lambda*h(:,k) + nW;

                end
                w2(:,:,j) = (theta(:,1:maxRuns));
                e2(:,j) = abs(delta).^2;
            end

            meanCount(barGammaIndex) = mean(count);
            
%             count = zeros(maxIt,1);

            w3 = mean(w2,3);
            wFinal(barGammaIndex,delay,:,L+1) = w3(:,end);

            e3(barGammaIndex,delay,:,L+1) = mean(e2,2);

        end
        % save(['.' filesep 'results' filesep 'results07.mat'],'wFinal','e3','meanCount');




    %     end

    end
    
end


% for i = 1:length(delayVector)
%     figure
% plot(10*log10((e3(i,:))))
% xlabel('Iterations','interpreter','latex');
% ylabel('MSE (dB)','interpreter','latex');
% 
% end


% for i = 1:adapFiltLength+10
%     plot(10*log10((e3(i,:,1))))
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
%     hold on;
% end


save(['.' filesep 'resultsMSE' filesep 'results49.mat'],'wFinal','e3','meanCount');




