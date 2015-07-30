clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'sigs-for-MC.xls');




IDX = 1:16;

data2 = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];
stat = jackPCA (data2);

subplot(1,4,1);

axis([-1 1 -1 1]);
axis square;
text(stat.c1m,stat.c2m,stat.ggnames);

hold on;
for ii = 1:length(stat.c1m)
    circle(stat.c1m(ii), stat.c2m(ii), stat.c1e(ii), stat.c2e(ii));
end



xlabel(['Principal Component 1 (' num2str(stat.exp(1),3) '%)']);
ylabel(['Principal Component 2 (' num2str(stat.exp(2),3) '%)']);
subplot(1,4,2);
hold on;


for ii = 1:length(stat.s1m)
    circle(stat.s1m(ii), stat.s2m(ii), stat.s1e(ii), stat.s2e(ii));
end


axis([-5 5 -4 4]);
axis square;
text(stat.s1m,stat.s2m,stat.stimuliS);
xlabel(['Principal Component 1 (' num2str(stat.exp(1),3) '%)']);
ylabel(['Principal Component 2 (' num2str(stat.exp(2),3) '%)']);
title('SKBR3');





IDX = 17:32;

data2 = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];
stat = jackPCA (data2);

subplot(1,4,3);
axis([-1 1 -1 1]);
axis square;
text(stat.c1m,stat.c2m,stat.ggnames);
hold on;
for ii = 1:length(stat.c1m)
    circle(stat.c1m(ii), stat.c2m(ii), stat.c1e(ii), stat.c2e(ii));
end
xlabel(['Principal Component 1 (' num2str(stat.exp(1),3) '%)']);
ylabel(['Principal Component 2 (' num2str(stat.exp(2),3) '%)']);
subplot(1,4,4);
hold on;

for ii = 1:length(stat.s1m)
    circle(stat.s1m(ii), stat.s2m(ii), stat.s1e(ii), stat.s2e(ii));
end


hold on;
axis([-5 5 -4 4]);
axis square;
text(stat.s1m,stat.s2m,stat.stimuliS);
xlabel(['Principal Component 1 (' num2str(stat.exp(1),3) '%)']);
ylabel(['Principal Component 2 (' num2str(stat.exp(2),3) '%)']);
title('PC9');



