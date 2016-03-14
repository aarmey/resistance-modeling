% For calculating Lowe synergy
% Data processing
function LoweAA
    files = dir('*-avg.csv');
    
    for ii = 1:10000
        IDX = mod(ii,length(files)) + 1;
        
        LoweAAfile(files(IDX).name)
    end
end

function LoweAAfile(name)
    [DD, viab] = loadData (name);

    % Do the actual fitting
    maxP = [ 2,  2,  0, 10, 10];
    minP = [ 0,  0,  0,-10,-10];
    
    optFun = @(x) runGreco(x,DD,viab);
    plotFun = @(x, flag) runGreco(x.x,DD,viab,1);
    
    statP = psoptimset('Display','iter', 'CompletePoll','on','CompleteSearch','on',...
        'PlotFcns',plotFun,'UseParallel',true);
    
    statG = gaoptimset('Display','iter','UseParallel',true,...
        'PlotFcns',@(~, state, ~) plotBestInd(DD, viab, state),...
        'PopulationSize',200); %#ok<NASGU>
    
    [x, f] = patternsearch(optFun, ...
        minP + (maxP - minP).*rand(size(minP)), [], [], [], [], minP, maxP, [], statP);
    %[x, f] = ga(optFun, length(minP), [], [], [], [], minP, maxP, [], statG);
    dlmwrite(['Nullout' name(1:(end-8)) '.csv'],[f x],'-append');
end

function fval = plotBestInd(DD, viab, state)
    [~, IDX] = min(state.Score);
    
    runGreco(state.Population(IDX,:),DD,viab,1);
    
    fval = state;
end

function fval = runGreco(xP,DD,viab, plott)
    if nargin < 4
        plott = 0;
    end

    IC501 = 10.^xP(1);
    IC502 = 10.^xP(2);
    a = xP(3);
    m1 = 1/(xP(4));
    m2 = 1/(xP(5));
    D1 = DD(:,1)/IC501;
    D2 = DD(:,2)/IC502;
    
    tempp = a*D1.*D2;
    mtemp = (m1+m2)/2;
    
    % Always scale E from 0 to 1
    ErrVal = @(EE) (D1./((EE ./ (1-EE)).^ m1)) ...
        + ( D2 ./ ((EE./(1-EE)) .^ m2)) ...
        + tempp./((EE./(1-EE)).^(mtemp))-1;
    
    fits = lsqnonlin(ErrVal,...
        0.5*ones(size(D1)), zeros(size(D1)), ones(size(D1)),...
        optimset('Display','off'));
    
    [x, fval] = fmincon(@(x) sum(((x(1)*fits' + x(2)) - viab).^2),...
        [mean(viab)/mean(fits) min(viab)], [], [], [], [], ...
        [0.5 0], [1.1 min(viab)*2],...
        [], optimset('Display','off'));
    
    if plott == 1
        hold off;
        plot(viab,'r.');
        hold on;
        plot(x(1)*fits' + x(2),'b');
        drawnow;
        
        fval = 0;
    end
end

function [DD, viab] = loadData (fileName)
    data = csvread(fileName);

    xAx = data(1,2:end);
    yAx = data(2:end,1);
    data(1,:) = []; data(:,1) = [];

    idx = 1;
    for ii = 1:length(xAx)
        for jj = 1:length(yAx)
            DD(idx,1) = xAx(ii); %#ok<AGROW>
            DD(idx,2) = yAx(jj); %#ok<AGROW>
            viab(idx) = data(jj,ii); %#ok<AGROW>
            idx = idx + 1;
        end
    end

    DD = DD + 0.001;
    viab = viab / max(viab);
end