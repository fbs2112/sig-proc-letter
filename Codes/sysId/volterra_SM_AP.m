%Volterra Affine Projection algorithm with kernel redundancy

clear;
clc;
close all;

maxRuns = 10000;
maxIt = 100;

N = 10;

auxMatrix = (ones(N));
[l1,l2] = find(auxMatrix);
adapFiltLength = (N^2+N);


h1 = [0.544 -0.252 0.593 0.236 -0.077 0.156 -0.5 0.025 -0.023 0.099].';
h2 = [-0.204 0.274 0.023 0.024 0.022 -0.274 -0.321 -0.070 0.712 0.433].';

ho = kron(h1,h2);

ho = [zeros(N,1);ho];


mu = 0.1;
signalPower = 1;
noisePower = 1e-3;
barGamma = sqrt(5*noisePower);

gamma = 1e-12;

L = 0:3;


wo(:,1) = zeros((N^2+N)/2 + N,1);
wo([4 11 15],1) = 1;


for LIndex = 1:length(L)
    count = zeros(maxIt,1);


    u = zeros(L(LIndex)+1,1);
    u(1) = 1;
    w2 = zeros(adapFiltLength,maxRuns,maxIt);

    for j = 1:maxIt
        
        j
     
        input = randn(maxRuns*2,1);
%         input = filter([1 0],[1 -0.9],input);
        input = input.*sqrt(signalPower/var(input));
        
        n = randn(maxRuns*2,1);
        n = n.*sqrt(noisePower/var(n));


        xTDL = flipud(buffer(input,N,N-1,'nodelay'));

        w = zeros(adapFiltLength,maxRuns);
        aux = 1;

        for k = L(LIndex)+1:maxRuns

            for l3 = 1:L(LIndex)+1
                counterAux = 1;
                xTDLAux = zeros((N^2+N)/2,1);
                
                for lIndex = 1:length(l1)
                    xTDLAux(lIndex,1) = xTDL(l1(lIndex),k+l3-1-L(LIndex))*xTDL(l2(lIndex),k+l3-1-L(LIndex));
                end
                
                xTDLConc(:,l3+k-L(LIndex)-1) = [xTDL(:,k+l3-1-L(LIndex));xTDLAux];
            end

            xAP = xTDLConc(:,k:-1:k-L(LIndex));
            

            d(k) = xAP(:,1).'*ho  + n(k);
            
            e(k) = d(k) - w(:,k)'*xAP(:,1);
            absoluteValueError = abs(e(k));
           
        
             if absoluteValueError >= barGamma
                 mu(k) = 1 - barGamma/absoluteValueError;
                 w(:,k+1) = w(:,k) + mu(k)*xAP*((xAP'*xAP+gamma*eye(L(LIndex)+1))\eye(L(LIndex)+1))*conj(e(k))*u;
             else
                 w(:,k+1) = w(:,k);
                count(j) = count(j)+1;
                 
             end
                        
        end
        w2(:,:,j) = conj(w(:,1:maxRuns));
        e2(:,j) = abs(e).^2;
    end

    meanCount(LIndex) = mean(count);

    w3 = mean(w2,3);
%     misalignment(:,L+1) = mean(misalignmentAux,2);
    
    e3(:,L(LIndex)+1) = mean(e2,2);

end
save(['.' filesep 'resultsTest' filesep 'results06.mat'],'e3','meanCount');


