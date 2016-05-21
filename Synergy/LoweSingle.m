% For calculating Lowe synergy
% Data processing
function LoweSingle
    subplot(2,2,1);

    [DD2, viab2] = loadDataFilter ('HCC-U0126-IN7-avg.csv', 2);
    fitAndPlot (DD2, viab2, 'r');
    
    [DD2, viab2] = loadDataFilter ('PC9-U0126-IN7-avg.csv', 2);
    fitAndPlot (DD2, viab2, 'k');
    
    [DD, viab] = loadDataFilter ('HOP62-U0126-IN7-avg.csv', 2);
    fitAndPlot (DD, viab, 'm');
    
    [DD, viab] = loadDataFilter ('HOP92-U0126-IN7-avg.csv', 2);
    fitAndPlot (DD, viab, 'g');
    
    [DD, viab] = loadDataFilter ('H322M-U0126-IN7-avg.csv', 2);
    fitAndPlot (DD, viab, 'y');
    
    [DD, viab] = loadDataFilter ('H522-U0126-IN7-avg.csv', 2);
    fitAndPlot (DD, viab, 'b');
    
    axis([-2.1 1.5 0 1.6]);
    xlabel('U0126 Log(uM)');
    title('U0126');
    
    
    subplot(2,2,3);
    [DD, viab] = loadDataFilter ('HCC-U0126-IN7-avg.csv', 1);
    fitAndPlot (DD, viab, 'r');
    
    [DD, viab] = loadDataFilter ('PC9-U0126-IN7-avg.csv', 1);
    fitAndPlot (DD, viab, 'k');
    
    [DD, viab] = loadDataFilter ('HOP62-U0126-IN7-avg.csv', 1);
    fitAndPlot (DD, viab, 'm');
    
    [DD, viab] = loadDataFilter ('HOP92-U0126-IN7-avg.csv', 1);
    fitAndPlot (DD, viab, 'g');
    
    [DD, viab] = loadDataFilter ('H322M-U0126-IN7-avg.csv', 1);
    fitAndPlot (DD, viab, 'y');
    
    [DD, viab] = loadDataFilter ('H522-U0126-IN7-avg.csv', 1);
    fitAndPlot (DD, viab, 'b');
    
    title('JNK-IN-7');
    axis([-2.1 1 0 2]);
end

function [DD, viab] = loadDataFilter (filename, drug)
    [DD, viab] = loadData (filename);
    
    DD = DD + 0.01;

    IDX = DD(:,mod(drug,2) + 1) < 0.2;
    DD = DD(IDX,drug)';
    
    IDXtwo = DD < 0.02;
    
    viab = viab(IDX);
    viab = viab / mean(viab(IDXtwo));
end

function fitAndPlot (DD, viab, color)
    E = @(Emax, IC50, m, D) Emax.*((D./IC50).^m)./(1 + (D./IC50).^m);


    options = gaoptimset('Display','none','PopulationSize',1000);


    xOut = ga(@(x) sum((E(x(:,1), 10.^x(:,2), x(:,3), DD) - viab).^2), ...
        3, [], [], [], [], [-2 -1 -40], [2 4 40], [], options);
    
    disp(xOut);

    DDx = unique(DD);
    for ii = 1:length(DDx)
        IDX = DD == DDx(ii);
        mmean(ii) = mean(viab(IDX));
        ssem(ii) = std(viab(IDX))/sqrt(sum(IDX)); 
    end
    

    fplot(@(xx) E(xOut(1), 10.^xOut(2), xOut(3), 10.^xx), ...
        [log10(min(DD)) log10(max(DD))], color);
    hold on;
    errorbar(log10(DDx), mmean, ssem, ['.' color]);
    axis([-2.1, 2, 0, max(viab)+0.1]);
    ylabel('Viability');
    drawnow;
end