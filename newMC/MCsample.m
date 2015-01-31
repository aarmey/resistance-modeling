clc; clear;
clear global;

fminPDF = @(x) -fitMC(10.^x);

SIZZ = 27;

inn = fmincon(@(x) -fminPDF(x), zeros(SIZZ,1),[],[],[],[], -2*ones(SIZZ,1), ones(SIZZ,1));

fitMC(10.^inn)


%%
rnd = slicesample(inn,1E4,'logpdf',fminPDF,'burnin',1000);


%%

MA = quantile(rnd, 0.95);
MI = quantile(rnd, 0.5);


bar(MI(19:end))