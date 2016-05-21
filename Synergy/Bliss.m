% For calculating Lowe synergy
% Data processing
clear; clc;
data = csvread('H522-U0126-IN7-avg.csv');
Y = data(2:end,1);
X = data(1,2:end);
data(1,:) = []; data(:,1) = [];

data = (data(1:8,:) + data(9:end,:))/2;
Y = Y(1:8);

data = data / data(1,1);
data = 1 - data;


effX = data(1,:);
effY = data(:,1);
effXm = ones(length(effY),1)*effX;
effYm = effY*ones(1,length(effX));



blissPred = effXm + effYm - effXm.*effYm;


cmin = min(min(data));
cmax = max(max(data));


subplot(1,3,1);
imagesc(data, [cmin cmax]);

subplot(1,3,2);
imagesc(blissPred, [cmin cmax]);
subplot(1,3,3);
imagesc(data - blissPred, [cmin cmax]);