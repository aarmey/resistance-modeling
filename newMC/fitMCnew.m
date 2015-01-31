function outE = fitMCnew(inn)

if size(inn,1) < size(inn,2)
    inn = inn';
end

global data;

if isempty(data)
    import bioma.data.*;
    data = double(DataMatrix('File', 'sigs-for-MC-red.xls'));
    data(isinf(data)) = nan;
    
    for cols = 0:3
        for rowBlocks = 0:3
            curData = data(rowBlocks*8+1:rowBlocks*8+8,cols*3+1:cols*3+3);
            curData = curData / mean(curData(1,:));
            
            difff = mean(curData(1,:)) - mean(curData(5,:));
            
            curData(5:end,:) = curData(5:end,:) + difff;
            
            
            data(rowBlocks*8+1:rowBlocks*8+8,cols*3+1:cols*3+3) = curData;
        end
    end

end


% IGFR, PDGFR, MET
% Akt, GSK, Erk, cJun
% This is the RTK to signal conversion
signalConv = reshape(inn(1:(3*4)),[3 4]);



facctor = @(x) [ones(1,4); ((x*ones(1,4) .* signalConv) + 1); ones(1,4); ((x*ones(1,4) .* signalConv) + 1)];

tripp = @(x) [x, x, x];



outt = [facctor(inn(13:15)); facctor(inn(16:18)); facctor(inn(19:21)); facctor(inn(22:24))];


outt = [tripp(outt(:,1)), tripp(outt(:,2)), tripp(outt(:,3)), tripp(outt(:,4))];


outE = sum(sum((outt - data).^2));

outE = outE + sum(sum(signalConv,2).^10 - 1);