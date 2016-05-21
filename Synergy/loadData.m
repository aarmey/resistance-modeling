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