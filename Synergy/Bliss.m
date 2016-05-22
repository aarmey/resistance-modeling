% For calculating Lowe synergy
% Data processing
function Bliss
    clear; clc;
    
    files = dir('*-avg.csv');
    
    for ii = 1:length(files)
        data = csvread(files(ii).name);
        [rsq(ii), varData(ii), residData(ii)] = plotData(data, 'JNK-IN-7 (uM)', ...
            'U0126 (uM)', 6, ii, files(ii).name);
    end
    
    subplot(6,6,31);
    bar(rsq)
    
    subplot(6,6,32);
    bar(residData./varData);
end

function [rsq, varData, varResid] = plotData(data, xlab, ylab, N, Nhere, name)
    Y = data(2:end,1);
    X = data(1,2:end);
    data(1,:) = []; data(:,1) = [];

    data = data(1:8,:) + data(9:end,:);
    Y = Y(1:8);

    data = 1 - (data / data(1,1));

    effX = data(1,:);
    effY = data(:,1);
    effXm = ones(length(effY),1)*effX;
    effYm = effY*ones(1,length(effX));



    blissPred = effXm + effYm - effXm.*effYm;
    
    resid = data - blissPred;

    crange = [-1 2];
    
    
    subplot(N,6,(Nhere-1)*3 + 1);
    imagesc(X, Y, -data, crange);
    colormap parula; %
    xlabel(xlab);
    ylabel(ylab);

    subplot(N,6,(Nhere-1)*3 + 2);
    imagesc(X, Y, -blissPred, crange);
    colormap parula;
    xlabel(xlab);
    ylabel(ylab);
    
    subplot(N,6,(Nhere-1)*3 + 3);
    imagesc(X, Y, -resid, crange);
    colormap parula;
    xlabel(xlab);
    ylabel(ylab);
    title(name);
    
    
    rsq = corrcoef(reshape(data,1,[]), reshape(blissPred,1,[]));
    rsq = rsq(1,2);
    
    varData = var(reshape(data,1,[]));
    varResid = var(reshape(resid,1,[]));
end