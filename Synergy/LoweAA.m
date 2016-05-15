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
    maxP{1} = [ 2,  2,  100, 10, 10];
    minP{1} = [ 0,  0, -100,-10,-10];
    
    maxP{2} = maxP{1};
    maxP{2}(3) = 0;
    minP{2} = minP{1};
    minP{2}(3) = 0;
    
    optFun = @(x) runGreco(x,DD,viab);
    %plotFun = @(x, flag) runGreco(x.x,DD,viab,1);
    
    statP = psoptimset('Display','final', 'CompletePoll','on','CompleteSearch','on',...
        'UseParallel',true); % 'PlotFcns',plotFun
    
    for ii = 1:2
        [x{ii}, f{ii}] = patternsearch(optFun, ...
            minP{ii} + (maxP{ii} - minP{ii}).*rand(size(minP{ii})),...
            [], [], [], [], minP{ii}, maxP{ii}, [], statP);
    end
    
    [~, resid1] = optFun(x{1});
    resid1 = resid1 - viab;
    
    [~, resid2] = optFun(x{2});
    resid2 = resid2 - viab;
    
    [~, pval] = vartest2(resid1,resid2);
    
    dlmwrite(['Nullout' name(1:(end-8)) '.csv'],[f{1} x{1} f{2} x{2} pval],'-append');
end

function [fval, dataFit] = runGreco(xP,DD,viab, plott)
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
    
    dataFit = x(1)*fits' + x(2);
    
    if plott == 1
        hold off;
        plot(viab,'r.');
        hold on;
        plot(dataFit,'b');
        drawnow;
        
        fval = 0;
    end
end