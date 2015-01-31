clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,2:end) = zscore(double(log(data(IDX,2:end))));
end

ddata = double(data);
ddata(isnan(ddata)) = 0;

%%

% cellLine = 1;
% 
% IDX = ddata(:,1) == cellLine;
% 
% h = glyphplot(ddata(IDX,2:end), 'glyph','star','varLabels',data.ColNames(2:end),'obslabels',data.RowNames(IDX));
% set(h(:,3),'FontSize',8);

names = {'EGF', 'HGF', 'PDGF', 'IGF', 'HRG', 'FGF', 'HBEGF', 'TGF', 'drug'};

selFun = @(x) ~cellfun(@isempty, strfind(data.RowNames,x));
drug = selFun('Erlot') | selFun('Vem') | selFun('Lapatinib');
inputt = [selFun(' EGF') selFun('HGF') selFun('PDGF') selFun('IGF') selFun('HRG') selFun('FGF') selFun('HBEGF') selFun('TGF') drug];

for cellLine = 1:4
    for ii = 2:10

        IDX = ddata(:,1) == cellLine;
        signal = ii;

        iddd = ones(size(ddata,1),1);
        facctor = @(f) inputt(IDX,:).*(iddd(IDX)*f);

        errFun = @(x) sum((sum(facctor(x),2) - ddata(IDX,signal)).^2);

        opts = optimoptions(@fmincon,'Display','off');
        x(cellLine,ii-1,:) = fmincon(errFun, zeros(size(names)), [], [], [], [], -10*ones(size(names)), 10*ones(size(names)), [], opts);
    end
    
%     subplot(2,2,cellLine);
%     plot(1:9,squeeze(x(cellLine,:,:)),'o');
%     legend(names);
end


for ii = 1:4
    subplot(2,2,ii)
    imagesc(squeeze(x(ii,:,:)))
end