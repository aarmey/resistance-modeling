function DM = getDataMat

import bioma.data.*;

[num, txt] = xlsread('FullSet.xls');

for ii = 2:size(txt,1)
    txtRow{ii-1} = [txt{ii,1} '-' txt{ii,2} '-' txt{ii,3}]; %#ok<AGROW>
end

DM = DataMatrix(num,txtRow,txt(1,4:end));

sites = {'AKT', 'MEK', 'JNK', 'ERK', 'cJUN', 'STAT3', 'GSK'};

for ii = 1:length(sites)
    neww = DataMatrix([nanmean(double(DM(:,sites{ii})),2), ...
        nanstd(double(DM(:,sites{ii})),0,2)/sqrt(3)], ...
        txtRow,{[sites{ii} 'm'],[sites{ii} 's']});
    DM(:,sites{ii}) = [];
    DM = horzcat(DM, neww); %#ok<AGROW>
end