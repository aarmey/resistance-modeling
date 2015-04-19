clc; clear;

import bioma.data.*;
data = DataMatrix('File', 'sigs-for-MC.xls');




IDX = [1 8:14];

data2 = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];
stat = jackPCA (data2);

subplot(4,2,1);
ploterr(stat.c1m,stat.c2m,stat.c1e,stat.c2e,'ro')
axis([-1 1 -1 1]);
axis square;
text(stat.c1m,stat.c2m,stat.ggnames);
xlabel('PC 1');
ylabel('PC 2');
subplot(4,2,2);
%hold on;
ploterr(stat.s1m,stat.s2m,stat.s1e,stat.s2e,'bo')
axis([-5 5 -4 4]);
axis square;
text(stat.s1m,stat.s2m,stat.stimuliS);
xlabel('PC 1');
ylabel('PC 2');
title('SKBR3');
















IDX = [15 22:28];

data2 = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];
stat = jackPCA (data2);

subplot(4,2,3);
ploterr(stat.c1m,stat.c2m,stat.c1e,stat.c2e,'ro')
axis([-1 1 -1 1]);
axis square;
text(stat.c1m,stat.c2m,stat.ggnames);
xlabel('PC 1');
ylabel('PC 2');
subplot(4,2,4);
%hold on;
ploterr(stat.s1m,stat.s2m,stat.s1e,stat.s2e,'bo')
axis([-5 5 -4 4]);
axis square;
text(stat.s1m,stat.s2m,stat.stimuliS);
xlabel('PC 1');
ylabel('PC 2');
title('BT474');




IDX = [29 36:42];

data2 = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];
stat = jackPCA (data2);

subplot(4,2,5);
ploterr(stat.c1m,stat.c2m,stat.c1e,stat.c2e,'ro')
axis([-1 1 -1 1]);
axis square;
text(stat.c1m,stat.c2m,stat.ggnames);
xlabel('PC 1');
ylabel('PC 2');
subplot(4,2,6);
%hold on;
ploterr(stat.s1m,stat.s2m,stat.s1e,stat.s2e,'bo')
axis([-5 5 -4 4]);
axis square;
text(stat.s1m,stat.s2m,stat.stimuliS);
xlabel('PC 1');
ylabel('PC 2');
title('P');







IDX = [43 49:54];

data2 = [data(IDX,1:3:end); data(IDX,2:3:end); data(IDX,3:3:end)];
stat = jackPCA (data2);

subplot(4,2,7);
ploterr(stat.c1m,stat.c2m,stat.c1e,stat.c2e,'ro')
axis([-1 1 -1 1]);
axis square;
text(stat.c1m,stat.c2m,stat.ggnames);
xlabel('PC 1');
ylabel('PC 2');
subplot(4,2,8);
hold on;
ploterr(stat.s1m,stat.s2m,stat.s1e,stat.s2e,'bo')
axis([-5 5 -4 4]);
axis square;
text(stat.s1m,stat.s2m,stat.stimuliS);
xlabel('PC 1');
ylabel('PC 2');
title('HCC');

