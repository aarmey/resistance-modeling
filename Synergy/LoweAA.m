% For calculating Lowe synergy
% Data processing
function LoweAA
    clc; clear;

    files = dir('*-avg.csv');
    
    for ii = 1:10000
        IDX = mod(ii,length(files)) + 1;
        
        LoweAAfile(files(IDX).name)
    end
end

function LoweAAfile(name)
    global DD viab;
    [DD, viab] = loadData (name);

    % Do the actual fitting
    maxP = [ 4,  4,  10, 10, 10,   1, min(viab)*2];
    minP = [-2, -2, -10,-10,-10, 0.2,           0];
    
    statP = gaoptimset('Display','iter','PlotFcns',@plotFun,'UseParallel',true); %
    
    [x{1}, f{1}] = ga(@(x) runGreco(DD, viab, x), length(minP), [],[],[],[], minP, maxP, [], statP);
    
    x{1}
    
    maxP(3) = 0; minP(3) = 0;
        
    [x{2}, f{2}] = ga(@(x) runGreco(DD, viab, x), length(minP), [],[],[],[], minP, maxP, [], statP);
    
    [~, resid1] = runGreco(DD, viab, x{1});
    [~, resid2] = runGreco(DD, viab, x{2});
    
    [~, pval] = vartest2(resid1 - viab, resid2 - viab);
    
    dlmwrite(['Nullout' name(1:(end-8)) '.csv'],...
        [f{1} x{1} f{2} x{2} pval],'-append');
end

function reform (DD, viab)
    X = unique(DD(:,1));
    Y = unique(DD(:,2));
    
    viabOut = zeros(length(X),length(Y));
    
    for ii = 1:length(X)
        for jj = 1:length(Y)
            IDXx = DD(:,1) == X(ii);
            IDXy = DD(:,2) == Y(jj);
            
            viabOut(ii,jj) = mean(viab(IDXx & IDXy));
        end
    end
    
    imagesc(X, Y, viabOut);
end

function state = plotFun(~, state, ~)
    global DD viab;
    
    [~, IDX] = min(state.Score);
    [~, dataFit] = runGreco(DD, viab, state.Population(IDX,:));
    
    state.Population(IDX,:)
    
    subplot(1,3,1);
    hold off;
    plot(viab,'ro');
    hold on;
    plot(dataFit,'b');
    drawnow;
    subplot(1,3,2);
    reform(DD, viab);
    subplot(1,3,3);
    reform(DD, dataFit);
end

function [fval, dataFit] = runGreco(DD, viab, xP)
    IC501 = 10.^xP(1);
    IC502 = 10.^xP(2);
    a = xP(3);
    m1 = 1/(xP(4));
    m2 = 1/(xP(5));
    D1 = DD(:,1)/IC501;
    D2 = DD(:,2)/IC502;
    Econ = xP(6);
    
    tempp = a*D1.*D2;
    mtemp = (m1+m2)/2;
    
    % Always scale E from 0 to 1
    ErrVal = @(EE) (D1./((EE./(Econ-EE)).^ m1)) ...
        + ( D2 ./ ((EE./(Econ-EE)) .^ m2)) ...
        + tempp./((EE./(Econ-EE)).^(mtemp))-1;
    
    try 
       fits = lsqnonlin(ErrVal,...
           0.5*ones(size(D1)), zeros(size(D1)), Econ*ones(size(D1)),...
           optimset('Display','off'));
    catch
       fval = 1E6;
       return;
    end
    
    dataFit = fits' + xP(7);
    
    fval = sum((dataFit - viab).^2);
end