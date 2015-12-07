clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,end) = zscore(double(data(IDX,end)));
    
    data(IDX,2:(end-1)) = zscore(log2(double(data(IDX,2:(end-1)))));
end

for cellLine = 1:4
    clustergram(data(data(:,1)==cellLine,2:(end-1)),'Standardize',0,'Colormap',redbluecmap)
end

