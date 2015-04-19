clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,end) = zscore(double(data(IDX,end)));
end


% clustergram(data(data(:,1)==1,2:(end-1)),'Standardize',1,'Colormap',redbluecmap)

idx = [];

for ii = 1:3
    idxx = nchoosek(2:7,ii);
    for jj = 1:size(idxx,1)
        idx(end+1,1:7) = false(1,7); %#ok<SAGROW>
        idx(end,idxx(jj,:)) = true; %#ok<SAGROW>
    end
end

vars = zeros(size(idx,1),7);
perf = zeros(4,size(idx,1));
ssize = zeros(1,size(idx,1));
sel = zeros(1,size(idx,1));
err = zeros(1,size(idx,1));

ddata = double(data);
ddata(isnan(ddata)) = 0;
data(:,:) = ddata;


parfor_progress(size(idx,1));

for ii = 1:(size(idx,1))
    idxes = find(idx(ii,:));
    
    for cellLine = 1:4
        IDX = ddata(:,1) == cellLine;

        mdl = fitlm(ddata(IDX,idxes),ddata(IDX,end),'interactions');

        perf(cellLine,ii) = mdl.ModelCriterion.AICc;
        err(cellLine,ii) = mdl.MSE;

        vars(ii,:) = idx(ii,:);
    end
    
    parfor_progress;
    
    ssize(ii) = sum(idx(ii,:));
end

parfor_progress(0);

%%

perf3 = sum(perf);
min(perf3)
perf3 = perf3 - min(perf3);
err1 = sum(err);

[~, sIDX1] = sort(perf3);
vars = vars(sIDX1,:);
err1 = err1(sIDX1);
ssize = ssize(sIDX1);
perf3 = perf3(sIDX1);
sel = sel(sIDX1);

subplot(3,1,1);
[hAx, hL1, hL2] = plotyy(1:length(err1), perf3, 1:length(err1), err1);
hL1.Marker = '.';
line([1 length(perf3)],[6 6]);
ylabel(hAx(1), 'AIC - AIC_{min}');
ylabel(hAx(2), 'Mean Squared Error');
axis(hAx(1), [1 length(err1) 0 max(perf3)]);
axis(hAx(2), [1 length(err1) 0 max(err1)]);
subplot(3,1,2);
imagesc(~vars(:,2:end)');
colormap('gray');

%%

figure(2);

for ii = 1:max(ssize)
    
    sBest(ii) = min(perf3(ssize == ii));
    bar(sBest);
end


