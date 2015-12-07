
function eNetRegression
    init;
    
    genIndivData;
    
    plotIndivData;
end

function init()
    global data labb cLine names;

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
    names = {'SKBR3','BT474','PC9','HCC827'};
end

function plotIndivData() 
    global names;
    
    load jacData;
    
    for ii = 1:length(names)
        
        for jj = 1:length(B{ii}) %#ok<USENS>
            if B{ii}(jj) > 0
                ccc = 'r';
            else
                ccc = 'k';
            end
            
            if pp{ii}(jj) < 0.05
                mmm = 'o';
            else
                mmm = '^';
            end
            
            scatter(jj,ii, 100*abs(B{ii}(jj))^2 + 0.00001, [ccc mmm], 'filled');
            hold on;
        end
    end

    axis([0 length(B{ii})+1 0 length(names)+1]);
    axis equal;
end


function genIndivData()
    global cLine data;
    

    for ii = 1:4
        curC = ii;
        X = x2fx(double(data(cLine == curC,1:(end-1))),'interactions');
        Y = double(data(cLine == curC,end));

        [B{ii}, bVar{ii}, pp{ii}, Ypred{ii}, pVarExp{ii}] = eNetJack(X,Y);
        
    end
    
    save('jacData','B','bVar','pp','Ypred','pVarExp');
end

function allLines
    % Model including all cell lines
    clc;
    
    labb = {'Offset','Akt','Erk','GSK','cJun','JNK','P38','cLine','Akt-Erk','Akt-GSK',...
        'Akt-cJun','Akt-JNK','Akt-P38','Akt-cLine','Erk-GSK','Erk-cJun','Erk-JNK','Erk-P38',...
        'Erk-cLine','GSK-cJun','GSK-JNK','GSK-P38','GSK-cLine','cJun-JNK','cJun-P38',...
        'cJun-cLine','JNK-P38','JNK-cLine','P38-cLine'};
    
    global data cLine;
    
    dataIn = [double(data(:,1:(end-1))) (cLine > 2)];

    X = x2fx(dataIn,'interactions');
    Y = zscore(double(data(:,end)));
    
    [beta, sigma, E, covB, logL] = mvregress(X, Y);
    
    nlogL = mvregresslike(X,Y,beta,sigma,'mvn');
    
    mvL = @(xIn) -mvregresslike(X,Y,xIn,sigma,'mvn');
    
    rnd = slicesample(beta, 3000, 'logpdf', mvL);
    
    plot(rnd)
    
    
    
    %[B, bVar, pp, Ypred, pVarExp] = eNetJack(X,Y);
    
    %pVarExp
    


%      [Btot, bVarTot, ppTot, Ypred, pVarExp] = eNetJack(X,Y);
% 
% 
%     IDX = ppTot < 0.05;
%     errorbar(Btot(IDX),bVarTot(IDX),'.');
%     text(1:sum(IDX),Btot(IDX),labb(IDX));
%     line([0 sum(IDX)+1],[0 0]);
end



% %%
% 
% 
% 
% for ii = 1:4
%     X = x2fx(double(data(cLine == curC,1:(end-1))),'interactions');
%     
%     xx = X(8,:);
%     yy = B{ii};
%     
%     subplot(2,2,ii)
%     title(names(ii));
%     plot(xx(pp{ii} < 0.05),yy(pp{ii} < 0.05),'o');
% end
% 
% 
% 
% 
% %%
% 
% X = x2fx(double(data(cLine < 3,1:(end-1))),'interactions');
% Y = double(data(cLine < 3,end));
% 
% [Btot, bVarTot, ppTot, Ypred, pVarExp] = eNetJack(X,Y);
% 
% 
% %%
% 
% function allLines
% 
%     X = x2fx(double(data(cLine > 2,1:(end-1))),'interactions');
%     Y = double(data(cLine > 2,end));
% 
%     [Btot, bVarTot, ppTot, Ypred, pVarExp] = eNetJack(X,Y);
%     
% end



% 
% 
% 
% %%
% 
% for ii = 1:2
%     subplot(2,1,ii);
%     
%     IDX = pp{ii} < 0.05;
% 
%     errorbar(B{ii}(IDX),bVar{ii}(IDX),'.');
%     text(1:sum(IDX),B{ii}(IDX),labb(IDX));
%     line([0 sum(IDX)+1],[0 0]);
%     axis([0 sum(IDX)+1 -1 0.5]);
% end
% %%
% 
% for ii = 1:4
%     subplot(2,2,ii);
%     
%     IDX = pp{ii} < 0.05;
%     
%     X = x2fx(double(data(cLine == ii,1:(end-1))),'interactions');
%     Y = double(data(cLine == ii,end));
%     
%     X = X(:,IDX);
%     XX = B{ii}(IDX);
%     
%     predOut = zscore(XX'*X');
%     
%     
%     bar([predOut; zscore(Y')]);
%     
% end
% 
% %%
% 
% for ii = 1:4
%     subplot(2,2,ii);
%     
%     IDX = pp{ii} < 0.05;
% 
%     errorbar(B{ii}(IDX),bVar{ii}(IDX),'.');
%     text(1:sum(IDX),B{ii}(IDX),labb(IDX));
%     line([0 sum(IDX)+1],[0 0]);
%     axis([0 sum(IDX)+1 -1 0.5]);
% end
% 
% 
% %%
% ii = 4;
% IDX = pp{ii} < 0.05;
% B{ii}(IDX)
% bVar{ii}(IDX)
% labb(IDX)'



function [Btot, Bvar, pp, Ypred, pVarExp] = eNetJack(X,Y)

parfor_progress(size(X,1));

parfor ii = 1:size(X,1)
    ii
    
    X2 = X;
    Y2 = Y;
    
    X2(ii,:) = [];
    Y2(ii) = [];
    
    Bvar(:,ii) = eNetOne(X2,Y2);
    
    parfor_progress;
end

parfor_progress(0);

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
    outter = B(:,FitInfo.IndexMinMSE);
end









