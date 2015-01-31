clc; clear;
clear global;

fminLog = @(x) fitMCnew(10.^x);
fminPDF = @(x) -fitMCnew(10.^x);

opts = saoptimset( 'Display', 'iter', 'TimeLimit',60*10, 'AnnealingFcn', @annealingboltz,...
    'TemperatureFcn',@temperatureboltz,'DisplayInterval',1000);

SIZZ = 24;

inn = simulannealbnd(fminLog, rand(SIZZ,1), -3*ones(SIZZ,1), [ones(12,1) 3*ones(SIZZ-12,1)], opts);





rnd = slicesample(inn,1E4,'logpdf',fminPDF,'burnin',1000);

