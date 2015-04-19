clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    
    data(IDX,2:end) = double(log(data(IDX,2:end)));
end





for cellLine = 1:4
    IDX = data(:,1) == cellLine;
    X = double(data(IDX,2:end));
    
    [R,~,RLO,RUP] = corrcoef(X(1:(size(X,1)/2),:));
    
    Rwo(:,cellLine) = R(1:6,end); %#ok<*SAGROW>
    RLwo(:,cellLine) = RLO(1:6,end);
    RUwo(:,cellLine) = RUP(1:6,end);
    
    [R,~,RLO,RUP] = corrcoef(X((size(X,1)/2+1):end,:));
    
    Rw(:,cellLine) = R(1:6,end);
    RLw(:,cellLine) = RLO(1:6,end);
    RUw(:,cellLine) = RUP(1:6,end);
end


