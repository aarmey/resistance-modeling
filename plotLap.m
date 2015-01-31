clc; clear;

import bioma.data.*;

data = DataMatrix('File', 'fullResistanceSet.xls');

cellLine = 1;
IDX = data(:,1) == cellLine;

clustergram(data(IDX,2:(end-1)),'Standardize',1,'Cluster',2,'Colormap','parula');


cellLine = 2;
IDX = data(:,1) == cellLine;

clustergram(data(IDX,2:(end-1)),'Standardize',1,'Cluster',2,'Colormap','parula');