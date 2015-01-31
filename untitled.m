clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

idx = [];

for ii = 1:3
    idxx = nchoosek(1:9,ii);
    for jj = 1:size(idxx,1)
        idx(end+1,1:9) = false(1,9);
        idx(end,idxx(jj,:)) = 1;
    end
end

vars = zeros(size(idx,1),9);
perf1 = zeros(1,size(idx,1));
perf2 = perf1;

for ii = 1:size(idx,1)
    idxes = find(idx(ii,:));

    mdl1 = LinearModel.fit(double(log10(data(1:16,idxes))),double(data(1:16,end)),'quadratic','VarNames', data.ColNames([idxes end]));
    mdl2 = LinearModel.fit(double(log10(data(17:end,idxes))),double(data(17:end,end)),'quadratic','VarNames', data.ColNames([idxes end]));
    
    perf1(ii) = mdl1.ModelCriterion.AICc;
    perf2(ii) = mdl2.ModelCriterion.AICc;
    vars(ii,idxes) = 1;
end

perf1 = perf1 - min(perf1);
perf2 = perf2 - min(perf2);

[perf1, sIDX1] = sort(perf1);
vars1 = vars(sIDX1,:);

[perf2, sIDX2] = sort(perf2);
vars2 = vars(sIDX2,:);

%%

subplot(2,2,1);
plot(perf1);
axis([1 length(perf1) 0 max(perf1)]);
subplot(2,2,3);
imagesc(vars1');

subplot(2,2,2);
plot(perf2);
axis([1 length(perf2) 0 max(perf2)]);
subplot(2,2,4);
imagesc(vars2');

%%
clc;
for ii = 1:size(vars2,2)
    p(ii,1) = sum(vars2(:,ii)' .* [1:size(vars2,1)]);
    p(ii,2) = sum(vars(:,ii)' .* [1:size(vars,1)]);
end

for ii = 1:10000
    vec = [ones([28 1]); zeros([84-28 1])];
    vec = vec(randperm(84));
    
    pRand(ii) = sum(vec' .* [1:size(vec,1)]);
end

h = zeros(size(p));
h(p < quantile(pRand,0.05)) = 1;

% clustergram(data(1:16,:)','Standardize','row','RowLabels',txt,'ColumnLabels',...
%     conditions(1:16),'Colormap',redbluecmap);
% 
% clustergram(data(17:end,:)','Standardize','row','RowLabels',txt,'ColumnLabels',...
%     conditions(17:end),'Colormap',redbluecmap);