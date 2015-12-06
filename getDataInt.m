function [T, data, mmean, sem, auc] = getDataInt(GF, Plan, CellLine, Site, pplot, ccolor)

if nargin < 5
    pplot = 0;
end

if nargin < 6
    ccolor = 'k';
end

[num, txt] = xlsread('FullSet.xls');

colIDX = find(strcmp(txt(1,:), Site)) - 3;

GFIDX = find(strcmp(txt(:,3),GF));
PlanIDX = find(strcmp(txt(:,2),Plan));
CellIDX = find(strcmp(txt(:,1),CellLine));

rowIDX = intersect(intersect(GFIDX, PlanIDX), CellIDX)-1;


data = num(rowIDX, colIDX);


T = num(rowIDX, 1);
mmean = mean(data, 2);
sem = std(data, 0, 2) ./ sqrt(size(data, 2));

mmeanUp = mmean + sem;
mmeanDown = mmean - sem;


for ii = 1:size(data,2)
    auc(ii) = trapz(T, data(:,ii)-data(1,ii));
end

if pplot
    h = fill( [T; flipud(T)],  [mmeanUp; flipud(mmeanDown)], ccolor);
    h.LineWidth = 0.01;
    hold on;
    alpha(0.25);
    plot(T, mmean, ccolor, 'LineWidth', 2);
    axis([min(T) max(T) 0 max(mmeanUp)*1.1]);
    ylabel('FI')
    xlabel('Time (min)');
end