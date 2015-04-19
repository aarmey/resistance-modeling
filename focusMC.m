clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'fullResistanceSet.xls');

idxes = [2 3 5 6];

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,2:end) = zscore(double(log(data(IDX,2:end))));
end

data = data(data(:,1) == 4,[idxes end]);


for ii = 1:100
    a = rand;
    b = rand;
    c = rand;
    d = rand;
    
    
    
    
    
end
    
    
