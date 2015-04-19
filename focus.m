clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'fullResistanceSet.xls');

idxes = [2 5];

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,end) = zscore(double(log(data(IDX,end))));
end



ddata = double(data);
ddata(isnan(ddata)) = 0;
data(:,:) = ddata;
    
for cellLine = 1:4
    
    IDX = data(:,1) == cellLine;
    
    mdl{cellLine} = fitlm(ddata(IDX,idxes),ddata(IDX,end),...
                'interactions','VarNames', data.ColNames([idxes end]));
            
    mdl{cellLine}
    
    ppp(cellLine,:) = -log10(mdl{cellLine}.Coefficients.pValue);
end