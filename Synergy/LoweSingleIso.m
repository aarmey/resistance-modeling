% For calculating Lowe synergy
% Data processing

function LoweSingleIso
    clear; clc;
    data = csvread('PC9-U0126-SP6-avg.csv');
    Y = data(:,1);
    X = data(1,:);
    data(1,:) = []; data(:,1) = [];
    
    
    scatter(X, Y, 5, data);
    
end

function fitAndPlot (DD, viab, color)
    E = @(Emax, IC50, m, D, vv) vv + Emax.*((D./IC50).^m)./(1 + (D./IC50).^m);


    options = gaoptimset('Display','none','PopulationSize',1000);


    xOut = ga(@(x) sum((E(x(:,1), 10.^x(:,2), x(:,3), DD, x(:,4)) - viab).^2), ...
        4, [], [], [], [], [0 -1 -40 0], [1 4 0 1], [], options);
    
    disp(xOut);

    DDx = unique(DD);
    for ii = 1:length(DDx)
        IDX = DD == DDx(ii);
        mmean(ii) = mean(viab(IDX));
        ssem(ii) = std(viab(IDX))/sqrt(sum(IDX)); 
    end
    

    fplot(@(xx) E(xOut(1), 10.^xOut(2), xOut(3), 10.^xx, xOut(4)), ...
        [log10(min(DD)) log10(max(DD))], color);
    hold on;
    errorbar(log10(DDx), mmean, ssem, ['.' color]);
    axis([-2.1, 2, 0, max(viab)+0.1]);
    ylabel('Viability');
    drawnow;
end