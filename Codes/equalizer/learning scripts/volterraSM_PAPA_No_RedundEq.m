%%
%Volterra Data Reuse Equalyzer 



clear;
clc;
close all;

maxRuns = 3000;
maxIt = 1000;

N = 8;

adapFiltLength = (N^2+N)/2 + N;

% adapFiltLength = N;


auxMatrix = triu(ones(N));
[l1,l2] = find(auxMatrix);

% h = [0.34-(0.27*1i) 0.87+(0.43*1i) 0.34-(0.21*1i)]; 

% L = round((adapFiltLength + length(h))/2);

signalPower = 1;
noisePower = signalPower/db2pow(30);
barGamma = 4*sqrt(5*noisePower);
gamma = 1e-12;
SNRAux = db2pow(30);
% lambda = 1;
% varW = 0;

% h(:,adapFiltLength + round((N+4)/2) + 10) = [1 0 0 1].';
% delayVector = round((length(h) + adapFiltLength)/2);
% delayVector = 4:8;

h = [1 0.2 -0.3];

delayVector = 3;

% delayVector = round((N+4)/2);

% for N = 8:8
%     adapFiltLength = (N^2+N)/2 + N;
% 
%     auxMatrix = triu(ones(N));
%     [l1,l2] = find(auxMatrix);

for delay = 1:length(delayVector)

    for L = 0:0
        count = zeros(maxIt,1);

        u = zeros(L+1,1);
        u(1) = 1;
        
        
        w2 = zeros(adapFiltLength,maxRuns,maxIt);
        for j = 1:maxIt
            j

            
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
              xAux2 = filter(h,1,pilot);
%               xAux2 = exp(-xAux2);
            
%               xAux2 = filter(h,1,pilot);
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
            xTDL = (buffer(pilot,N,N-1,'nodelay'));

            w = zeros((N^2+N)/2 + N,maxRuns) + 1e-6;

%                 w = zeros(N,maxRuns);


            for k = (adapFiltLength + delayVector(delay) + L + 10):maxRuns
%                     xAux2 = filter(h(:,k),1,xTDL(:,k));
% %                     xAux2 = xAux2.^2;
%                     
%                     
%                     n = randn(N,1) + randn(N,1)*1i;
%                     powerSignal = xAux2'*xAux2./(N);
%                     powerNoiseAux = n'*n/(N);
%                     powerNoise = (powerSignal/SNRAux);
%                     n = n.*sqrt(powerNoise/powerNoiseAux);
% 
%                     xAux = xAux2 + n;





                xAP = zeros(N,L+1);

                for l = 0:L
                    xAP(:,l+1) = xAux(k-l:-1:k-N+1-l);
                end

%                     xAP = flipud(xAux);

                xTDLConc = zeros(adapFiltLength,L+1);

                for l3 = 1:L+1
                    xTDLAux = zeros((N*N+N)/2,1);



                    for lIndex = 1:length(l1)
                        xTDLAux(lIndex,1) = xAP(l1(lIndex),l3)*(xAP(l2(lIndex),l3));

%                                 xTDLAux(counterAux,1) = xTDL(l1,k+l3-1)*xTDL(l2,k+l3-1);
                    end


                    xTDLConc(:,l3) = [xAP(:,l3);xTDLAux];
    %                 xTDLConc(:,l3+k-L-1) = [xTDL(:,k+l3-1);xTDLAux];
                end

%                     xTDLConc = xAP;

%                 xAP = xTDLConc(:,k:-1:k-L);


                d(k) = (pilot(-delayVector(delay) + k + 1)); 
    %             d(k) = (pilot(-L2 + adapFiltLength+1 - L)); 

                e(k) = d(k) - w(:,k)'*xTDLConc(:,1);
%                 if abs(e(k)) > 100
%                     abs(e(k));
%                 end
                absoluteValueError = abs(e(k));

                if absoluteValueError > barGamma
                    mu(k) = 1 - barGamma/absoluteValueError;
%                     mu(k) = 0.2;
                    G(:,:,k) = diag(((1 - kappa*mu(k))/adapFiltLength) + (kappa*mu(k)*abs(w(:,k))/norm(w(:,k),1)));
                    w(:,k+1) = w(:,k) + mu(k)*G(:,:,k)*xTDLConc*((xTDLConc'*G(:,:,k)*xTDLConc+gamma*eye(L+1))\eye(L+1))*conj(e(k))*u;
%                     b(k) = ((xTDLConc'*xTDLConc+gamma*eye(L+1)));
                else
                    mu(k) = 0;
                    count(j) = count(j)+1;
                    w(:,k+1) = w(:,k);
                end
%                     nW = randn(length(h(k)),1);
%                     nW = nW.*sqrt(varW/var(nW));
%                     h(:,k+1) = lambda*h(:,k) + nW;

            end
            w2(:,:,j) = conj(w(:,1:maxRuns));
            e2(:,j) = abs(e).^2;
        end

        meanCount(delay,L+1) = mean(count);

        w3 = mean(w2,3);
        wFinal(delay,:,L+1) = w3(:,end);

        e3(delay,:,L+1) = mean(e2,2);

    end
    % save(['.' filesep 'results' filesep 'results07.mat'],'wFinal','e3','meanCount');




%     end

end


for i = 1:length(delayVector)
    figure
plot(10*log10((e3(i,:))))
xlabel('Iterations','interpreter','latex');
ylabel('MSE (dB)','interpreter','latex');

end


% for i = 1:adapFiltLength+10
%     plot(10*log10((e3(i,:,1))))
%     xlabel('Iterations','interpreter','latex');
%     ylabel('MSE (dB)','interpreter','latex');
%     hold on;
% end


save(['.' filesep 'resultsMSE' filesep 'results26.mat'],'wFinal','e3','meanCount');




