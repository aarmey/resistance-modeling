clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'sigs-for-MC.xls');

IDX = [1, 9, 17, 25, 33, 41, 49, 56];

data = data(IDX,:);
ddata = double(data);

data.ColNames

ddata(ddata < 0) = 2;
data = log10(ddata);


colIDX = @(x,z) (x*2+z):(x*2+z);

sem = @(x) std(x) / sqrt(length(x));

for ii = 0:5
    for jj = 0:3
        
        rowIDX = (ii*3+1):(ii*3+3);
        
        dataOut(ii+1,jj+1) = nanmean(double(data(colIDX(jj,2),rowIDX))) - nanmean(double(data(colIDX(jj,1),rowIDX)));
        
        dataSEM(ii+1,jj+1) = sem(double(data(colIDX(jj,2),rowIDX))) + sem(double(data(colIDX(jj,1),rowIDX)));
        
    end
end
    



errorbar(dataOut,dataSEM)