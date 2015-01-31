clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

skbr3sig = data(1:16,1:9);

bt474sig = data(17:end,1:9);

for ii = 1:size(skbr3sig,2)
    skbr3sig(:,ii) = double(skbr3sig(:,ii)) / double(skbr3sig(1,ii));
    
    bt474sig(:,ii) = double(bt474sig(:,ii)) / double(bt474sig(1,ii));
end



skbrViab = double(data(1:16,10));
bt474Viab = double(data(17:end,10));



f=fit(double(bt474sig(:,1)),bt474Viab,'poly2');
resid = f(double(bt474sig(:,1))) - bt474Viab;


plot(double(bt474sig(:,7)) - double(bt474sig(:,8)) - double(bt474sig(:,9)), resid,'o')

[r, p] = corr(double(bt474sig(:,7)) - double(bt474sig(:,8)) - double(bt474sig(:,9)), resid)
% 
% 
% subplot(1,2,1);
% scatter(,,25,,'filled')
% xlabel('Viability');
% ylabel('P38')
% 
% subplot(1,2,2);
% scatter(bt474Viab,double(bt474sig(:,9)),25,double(bt474sig(:,1)),'filled')
% xlabel('Viability');
% ylabel('NFKB')