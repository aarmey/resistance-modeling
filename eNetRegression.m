clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,end) = (double(data(IDX,end)));
end

cLine = double(data(:,1));
data(:,1) = [];

labb = {'Offset','Akt','Erk','GSK','cJun','JNK','P38','Akt-Erk','Akt-GSK',...
    'Akt-cJun','Akt-JNK','Akt-P38','Erk-GSK','Erk-cJun','Erk-JNK','Erk-P38',...
    'GSK-cJun','GSK-JNK','GSK-P38','cJun-JNK','cJun-P38','JNK-P38'};





%% Model including all cell lines

X = x2fx(double(data(:,1:(end-1))),'interactions');
Y = double(data(:,end));

% [B, FitInfo] = lasso(X,zscore(Y),'CV',26,'Options',statset('UseParallel',true));
% lassoPlot(B, FitInfo);


 [Btot, bVarTot, ppTot, Ypred, pVarExp] = eNetJack(X,Y);

%%

IDX = ppTot < 0.05;
errorbar(Btot(IDX),bVarTot(IDX),'.');
text(1:sum(IDX),Btot(IDX),labb(IDX));
line([0 sum(IDX)+1],[0 0]);

%%
for ii = 1:4
    curC = ii;
    X = x2fx(double(data(cLine == curC,1:(end-1))),'interactions');
    Y = double(data(cLine == curC,end));

    [B{ii}, bVar{ii}, pp{ii}, Ypred{ii}, pVarExp{ii}] = eNetJack(X,Y);
end

%%

names = {'SKBR3','BT474','PC9','HCC827'};

for ii = 1:4
    X = x2fx(double(data(cLine == curC,1:(end-1))),'interactions');
    
    xx = X(8,:);
    yy = B{ii};
    
    subplot(2,2,ii)
    title(names(ii));
    plot(xx(pp{ii} < 0.05),yy(pp{ii} < 0.05),'o');
end




%%

X = x2fx(double(data(cLine < 3,1:(end-1))),'interactions');
Y = double(data(cLine < 3,end));

[Btot, bVarTot, ppTot, Ypred, pVarExp] = eNetJack(X,Y);


%%



X = x2fx(double(data(cLine > 2,1:(end-1))),'interactions');
Y = double(data(cLine > 2,end));

[Btot, bVarTot, ppTot, Ypred, pVarExp] = eNetJack(X,Y);






%%

for ii = 1:2
    subplot(2,1,ii);
    
    IDX = pp{ii} < 0.05;

    errorbar(B{ii}(IDX),bVar{ii}(IDX),'.');
    text(1:sum(IDX),B{ii}(IDX),labb(IDX));
    line([0 sum(IDX)+1],[0 0]);
    axis([0 sum(IDX)+1 -1 0.5]);
end
%%

for ii = 1:4
    subplot(2,2,ii);
    
    IDX = pp{ii} < 0.05;
    
    X = x2fx(double(data(cLine == ii,1:(end-1))),'interactions');
    Y = double(data(cLine == ii,end));
    
    X = X(:,IDX);
    XX = B{ii}(IDX);
    
    predOut = zscore(XX'*X');
    
    
    bar([predOut; zscore(Y')]);
    
end

%%

for ii = 1:4
    subplot(2,2,ii);
    
    IDX = pp{ii} < 0.05;

    errorbar(B{ii}(IDX),bVar{ii}(IDX),'.');
    text(1:sum(IDX),B{ii}(IDX),labb(IDX));
    line([0 sum(IDX)+1],[0 0]);
    axis([0 sum(IDX)+1 -1 0.5]);
end


%%
ii = 4;
IDX = pp{ii} < 0.05;
B{ii}(IDX)
bVar{ii}(IDX)
labb(IDX)'
