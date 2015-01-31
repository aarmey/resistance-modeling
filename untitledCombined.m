clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,2:end) = zscore(double(log(data(IDX,2:end))));
end

idx = [];

for ii = 1:4
    idxx = nchoosek(2:10,ii);
    for jj = 1:size(idxx,1)
        idx(end+1,1:10) = false(1,10);
        idx(end,idxx(jj,:)) = true;
    end
end

vars = zeros(size(idx,1),10);
perf = zeros(4,size(idx,1));

ddata = double(data);
ddata(isnan(ddata)) = 0;
data(:,:) = ddata;

parfor_progress(size(idx,1));


ddata(find(ddata == 0)) = rand(90,1)*0.01;

parfor ii = 1:size(idx,1)
    idxes = find(idx(ii,:));
    
    for cellLine = 1:4
        IDX = ddata(:,1) == cellLine; %#ok<PFBNS>

        try
            mdl = LinearModel.fit(ddata(IDX,idxes),ddata(IDX,end),'linear');

            perf(cellLine,ii) = mdl.ModelCriterion.AIC;
            err(cellLine,ii) = mdl.MSE;
            sel(ii) = 1;
            
            tmp = LinearModel.fit(ddata(IDX,idxes),ddata(IDX,end),'purequadratic');

            if (tmp.ModelCriterion.AIC < perf(cellLine,ii))
                perf(cellLine,ii) = tmp.ModelCriterion.AIC;
                err(cellLine,ii) = tmp.MSE;
                sel(ii) = 3;
            end
%             
%             tmp = LinearModel.fit(ddata(IDX,idxes),ddata(IDX,end),'interactions');
% 
%             if (tmp.ModelCriterion.AIC < perf(cellLine,ii))
%                 perf(cellLine,ii) = tmp.ModelCriterion.AIC;
%                 err(cellLine,ii) = tmp.MSE;
%                 sel(ii) = 2;
%             end
%             
%             tmp = LinearModel.fit(ddata(IDX,idxes),ddata(IDX,end),'quadratic');
% 
%             if (tmp.ModelCriterion.AIC < perf(cellLine,ii))
%                 perf(cellLine,ii) = tmp.ModelCriterion.AIC;
%                 err(cellLine,ii) = tmp.MSE;
%                 sel{ii} = 'quadratic';
%             end
            
        catch
            perf(cellLine,ii) = nan;
            err(cellLine,ii) = nan;
        end
    end
    
    parfor_progress;
    
    ssize(ii) = sum(idx(ii,:));
end

parfor_progress(0);

%%

perf3 = sum(perf);
perf3 = perf3 - min(perf3);
err1 = sum(err);

[~, sIDX1] = sort(perf3);
idx = idx(sIDX1,:);
err1 = err1(sIDX1);
ssize = ssize(sIDX1);
perf3 = perf3(sIDX1);


subplot(3,1,1);
hAx = plotyy(0:length(err1)-1, perf3, 0:length(err1)-1, err1);
line([1 length(perf3)],[6 6]);
ylabel(hAx(1), 'AIC - AIC_{min}');
ylabel(hAx(2), 'Mean Squared Error');
axis(hAx(1), [1 length(err1) 0 max(perf3)]);
axis(hAx(2), [1 length(err1) 0 max(err1)]);
subplot(3,1,2);
imagesc(idx(:,2:end)');
subplot(3,1,3);
hAx = plotyy(0:length(err1)-1, ssize, 0:length(err1)-1, sel);
axis(hAx(1), [1 length(err1) 0 max(ssize)+1]);
axis(hAx(2), [1 length(err1) 0 max(sel)+1]);