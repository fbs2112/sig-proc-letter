function [L,DHat] = MattDavies(hessianMatrix)

n = size(hessianMatrix,1);

L = zeros(n);
DHat = zeros(n);

if hessianMatrix(1,1) > 0
    h00 = hessianMatrix(1,1);
else
    h00 = 1;
end

for k = 2:n
    m = k-1;
    L(m,m) = 1;
    
    if hessianMatrix(m,m) <= 0
        hessianMatrix(m,m) = h00;
    end
    
    for i = k:n
        L(i,m) = -hessianMatrix(i,m)/hessianMatrix(m,m);
        hessianMatrix(i,m) = 0;
     
        for j = k:n
            hessianMatrix(i,j) = hessianMatrix(i,j)+L(i,m)*hessianMatrix(m,j);
        end
    end
    
    if hessianMatrix(k,k) > 0 && hessianMatrix(k,k) < h00
        h00 = hessianMatrix(k,k);
    end
end

L(n,n) = 1;

if hessianMatrix(n,n) <=0
    hessianMatrix(n,n)= h00;
end

for i = 1:n
    DHat(i,i) = hessianMatrix(i,i);
end
    