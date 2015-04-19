function [Btot, Bvar, pp, Ypred, pVarExp] = eNetJack(X,Y)

parfor ii = 1:size(X,1)
    ii
    
    X2 = X;
    Y2 = Y;
    
    X2(ii,:) = [];
    Y2(ii) = [];
    
    Bvar(:,ii) = eNetOne(X2,Y2);
end

for jj = 1:size(Bvar,1)
    pp(jj) = signtest(Bvar(jj,:));
end

Btot = mean(Bvar,2);
Bvar = std(Bvar,0,2) / sqrt(size(X,1));

Ypred = Btot'*X';
pVarExp = (1 - var(zscore(Ypred) - zscore(Y'))) / var(zscore(Y));

end


function outter = eNetOne(X,Y)

[B, FitInfo] = lasso(X,Y,'CV',13);
outter = B(:,FitInfo.Index1SE);

end