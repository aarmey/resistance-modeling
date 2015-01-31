clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'fullResistanceSet.xls');

idxes = [2 3 6];

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,2:end) = zscore(double(log(data(IDX,2:end))));
end




%%
ddata = double(data);
ddata(isnan(ddata)) = 0;
data(:,:) = ddata;
    
for cellLine = 1:4
    
    IDX = data(:,1) == cellLine;
    
    mdl{cellLine} = LinearModel.fit(ddata(IDX,idxes),ddata(IDX,end),...
                'linear','VarNames', data.ColNames([idxes end]));
            
    mdl{cellLine}
end