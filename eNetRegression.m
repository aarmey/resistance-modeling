function eNetRegression
    clc; clear;
    init;
    genIndivData;
    %plotIndivData;
end

function init()
    global data labb cLine names;

    import bioma.data.*;

    data = DataMatrix('File', 'fullResistanceSet.xls');

    for cellLine = 1:4
        IDX = data(:,1) == cellLine;

        data(IDX,end) = zscore(double(data(IDX,end)));
    end

    cLine = double(data(:,1));
    data(:,1) = [];

    labb = {'Offset','Akt','Erk','GSK','cJun','JNK','P38','Akt-Erk','Akt-GSK',...
        'Akt-cJun','Akt-JNK','Akt-P38','Erk-GSK','Erk-cJun','Erk-JNK','Erk-P38',...
        'GSK-cJun','GSK-JNK','GSK-P38','cJun-JNK','cJun-P38','JNK-P38'};
    names = {'SKBR3','BT474','PC9','HCC827'};
end

function plotIndivData() 
    global names B bVar Ypred pVarExp; %#ok<NUSED>
    
    for ii = 1:2% length(names)
        
        for jj = 1:length(B{ii})
            if B{ii}(jj) > 0
                ccc = 'r';
            else
                ccc = 'k';
            end
            
            scatter(jj,ii, 100*abs(B{ii}(jj))^2 + 1.0E-7, [ccc 'o'], 'filled');
            hold on;
        end
    end

    axis([0 length(B{ii})+1 0 length(names)+1]);
    axis equal;
end

function genIndivData()
    global cLine data;
  
    labb = {'Offset','Akt','Erk','GSK','cJun','JNK','P38','cLine','Akt-Erk','Akt-GSK',...
        'Akt-cJun','Akt-JNK','Akt-P38','Akt-cLine','Erk-GSK','Erk-cJun','Erk-JNK','Erk-P38',...
        'Erk-cLine','GSK-cJun','GSK-JNK','GSK-P38','GSK-cLine','cJun-JNK','cJun-P38',...
        'cJun-cLine','JNK-P38','JNK-cLine','P38-cLine'};
    
    VarNames = {'Akt', 'Erk', 'GSK','cJun','JNK','P38', 'Viability'};
    
    dataIn = double(data(cLine < 3,1:(end-1)));

    X = x2fx(dataIn,'interactions');
    %X = dataIn;
    Y = zscore(double(data(cLine < 3,end)));
    
    [beta, sigma, E, covB, logL] = mvregress(X, Y);
    
    nlogL = mvregresslike(X,Y,beta,sigma,'mvn');
    
    mvL = @(xIn) -mvregresslike(X,Y,xIn,sigma,'mvn');
    
    rnd = slicesample(beta, 3000, 'logpdf', mvL);
    
    hist(rnd)
    
    
end




% function genIndivData()
%     global cLine data B bVar Ypred pVarExp;
%     
%     ii = 1;
%     X = x2fx((double(data(cLine > 2,1:(end-1)))),'interactions'); %  == curC
%     Y = zscore(double(data(cLine > 2,end)));
% 
%     [B{ii}, bVar{ii}, Ypred{ii}, pVarExp{ii}] = eNetJack(X,Y);
% 
%     ii = 2;
%     X = x2fx((double(data(cLine < 3,1:(end-1)))),'interactions'); %  == curC
%     Y = zscore(double(data(cLine < 3,end)));
% 
%     [B{ii}, bVar{ii}, Ypred{ii}, pVarExp{ii}] = eNetJack(X,Y);
% end

function [Btot, Bvar, Ypred, pVarExp] = eNetJack(X,Y)
    [Bvar, FitInfo] = lasso(X,Y,'CV',13);
    FitInfo
    Btot = Bvar(:,FitInfo.IndexMinMSE);

    Ypred = Btot'*X';
    pVarExp = (1 - var(zscore(Ypred) - zscore(Y'))) / var(zscore(Y));
end