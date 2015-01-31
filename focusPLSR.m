clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'sigs-for-MC.xls');


IDX = 1:16;

data = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];

data(:,:) = zscore(double(data));
    
%sigs = (data(1:16,2:10));

jackMDS (data)


%plot(XL(:,1),XL(:,2),'o');
%hold on;
%plot(XS(:,1),XS(:,2),'go');



% for ii = 1:16
%     text(YS(ii,1),YS(ii,2),data.rownames{ii});
% end

%plot(YL(1),YL(2),'ro');
%axis([-3 3 -3 3]);




% ploterr(XL(:,1),XL(:,2),XLjk(:,1),XLjk(:,2),'o')
% hold on;
% for ii = 1:9
%     text(XL(ii,1),XL(ii,2),data.colnames{ii+1});
% end


% sigs = double(data(25:(25+7),2:end));
% [R,P,RLO,RUP] = corrcoef(sigs);
% 
% R = R(10,1:9);
% RLO = RLO(10,1:9);
% RUP = RUP(10,1:9);
% P = P(10,1:9);
