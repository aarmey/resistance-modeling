%function outE = fitMC(inn)
clear global;
clc; clear;
% Inn should be length 30
inn = ones(30,1);

global data;

if isempty(data)
    import bioma.data.*;
    data = double(DataMatrix('File', 'fullResistanceSetRed2.xls'));
    data(:,1:2) = [];
end




sigMat = ones(3,6);

















% facctor = reshape(inn(1:18),[3 6]);
% exxpress = inn(19:end);
% 
% outE = 0;
% 
% % Loop through each block
% for ii = 0:5
%     curData = data((ii*4+1):(ii*4+4),3:end);
%     
%     facctor2 = facctor;
%     
%     for jj = 2:size(curData,1)
%         curData(jj,:) = (curData(jj,:) - curData(1,:));
%         
%         facctor2(jj-1,:) = facctor2(jj-1,:) * exxpress(floor(ii/2)*3 + jj - 1);
%     end
%     
%     curData(1,:) = [];
% 
%     outE = outE + sum(sum((curData - facctor2).^2));
% end
% 
% outE = outE + 10*norm(sum(facctor,2) - 1);
% 
% %end

