function outE = fitMC2(inn)

global data;

if isempty(data)
    import bioma.data.*;
    data = log(double(DataMatrix('File', 'fullResistanceSet.xls')));
    data(isinf(data)) = nan;
end

% EGFR, HER2, HER3, HER4, IGFR, MET, PDGFR, FGFR
% Akt, Erk, GSK, P70, cJun, STAT3, JNK, P38, NFkB
% This is the RTK to signal conversion
signalConv = reshape(inn(1:(9*8)),[8 9]);

outE = 1000*norm(sum(signalConv,2) - 1);

SKBR_TGF = inn(73:74);
SKBR_IGF = inn(75);
SKBR_LAP_HER2 = inn(76);
SKBR_LAP_HER3 = inn(77);
SKBR_PDGF = inn(78);
SKBR_MET = inn(79);
SKBR_HBEGF = inn(80:83);
SKBR_EGF = inn(84:85);
SKBR_HRG = inn(86:87);

SKBR3_LAP_Vec = SKBR_LAP_HER2*signalConv(2,:) + SKBR_LAP_HER3*signalConv(3,:);

SKBR3_data = zeros(16,9);
SKBR3_data(2,:) = SKBR_TGF(2)*signalConv(2,:) + SKBR_TGF(1)*signalConv(1,:);
SKBR3_data(3,:) = SKBR_HRG(1)*signalConv(2,:) + SKBR_HRG(2)*signalConv(3,:);
SKBR3_data(4,:) = SKBR_IGF*signalConv(5,:);
SKBR3_data(5,:) = SKBR_PDGF*signalConv(7,:);
SKBR3_data(6,:) = SKBR_EGF(1)*signalConv(1,:) + SKBR_EGF(2)*signalConv(2,:);
SKBR3_data(7,:) = SKBR_MET*signalConv(6,:);
SKBR3_data(8,:) = SKBR_HBEGF(1)*signalConv(1,:) + SKBR_HBEGF(2)*signalConv(2,:) + SKBR_HBEGF(3)*signalConv(3,:) + SKBR_HBEGF(4)*signalConv(4,:);
SKBR3_data(9,:) = -SKBR3_LAP_Vec;
SKBR3_data(10,:) = SKBR_TGF(1)*signalConv(1,:) - SKBR3_LAP_Vec;
SKBR3_data(11,:) = SKBR_HRG(2)*signalConv(3,:) - SKBR3_LAP_Vec;
SKBR3_data(12,:) = SKBR_IGF*signalConv(5,:) - SKBR3_LAP_Vec;
SKBR3_data(13,:) = SKBR_PDGF*signalConv(7,:) - SKBR3_LAP_Vec;
SKBR3_data(14,:) = SKBR_EGF(1)*signalConv(1,:) - SKBR3_LAP_Vec;
SKBR3_data(15,:) = SKBR_MET*signalConv(6,:) - SKBR3_LAP_Vec;
SKBR3_data(16,:) = SKBR_HBEGF(1)*signalConv(1,:) + SKBR_HBEGF(3)*signalConv(3,:) + SKBR_HBEGF(4)*signalConv(4,:) - SKBR3_LAP_Vec;

SKBR3_data(1:16,:) = SKBR3_data(1:16,:) - data(1:16,2:10);

for ii = 1:16
    outE = outE + norm(SKBR3_data(ii,:));
end


%%%%% Following section is actually BT474
SKBR_TGF = inn(88:89);
SKBR_IGF = inn(90);
SKBR_LAP_HER2 = inn(91);
SKBR_LAP_HER3 = inn(92);
SKBR_PDGF = inn(93);
SKBR_MET = inn(94);
SKBR_HBEGF = inn(95:98);
SKBR_EGF = inn(99:100);
SKBR_HRG = inn(101:102);

SKBR3_LAP_Vec = SKBR_LAP_HER2*signalConv(2,:) + SKBR_LAP_HER3*signalConv(3,:);

SKBR3_data = zeros(16,9);
SKBR3_data(2,:) = SKBR_TGF(2)*signalConv(2,:) + SKBR_TGF(1)*signalConv(1,:);
SKBR3_data(3,:) = SKBR_HRG(1)*signalConv(2,:) + SKBR_HRG(2)*signalConv(3,:);
SKBR3_data(4,:) = SKBR_IGF*signalConv(5,:);
SKBR3_data(5,:) = SKBR_PDGF*signalConv(7,:);
SKBR3_data(6,:) = SKBR_EGF(1)*signalConv(1,:) + SKBR_EGF(2)*signalConv(2,:);
SKBR3_data(7,:) = SKBR_MET*signalConv(6,:);
SKBR3_data(8,:) = SKBR_HBEGF(1)*signalConv(1,:) + SKBR_HBEGF(2)*signalConv(2,:) + SKBR_HBEGF(3)*signalConv(3,:) + SKBR_HBEGF(4)*signalConv(4,:);
SKBR3_data(9,:) = -SKBR3_LAP_Vec;
SKBR3_data(10,:) = SKBR_TGF(1)*signalConv(1,:) - SKBR3_LAP_Vec;
SKBR3_data(11,:) = SKBR_HRG(2)*signalConv(3,:) - SKBR3_LAP_Vec;
SKBR3_data(12,:) = SKBR_IGF*signalConv(5,:) - SKBR3_LAP_Vec;
SKBR3_data(13,:) = SKBR_PDGF*signalConv(7,:) - SKBR3_LAP_Vec;
SKBR3_data(14,:) = SKBR_EGF(1)*signalConv(1,:) - SKBR3_LAP_Vec;
SKBR3_data(15,:) = SKBR_MET*signalConv(6,:) - SKBR3_LAP_Vec;
SKBR3_data(16,:) = SKBR_HBEGF(1)*signalConv(1,:) + SKBR_HBEGF(3)*signalConv(3,:) + SKBR_HBEGF(4)*signalConv(4,:) - SKBR3_LAP_Vec;

SKBR3_data(1:16,:) = SKBR3_data(1:16,:) - data((1:16)+16,2:10);

for ii = 1:16
    outE = outE + norm(SKBR3_data(ii,:));
end
%%%%% End BT474



%%%% START PC9 Section %%%%%

PC9_EGF = inn(103:104);
PC9_HGF = inn(105);
PC9_IGF = inn(106);
PC9_PDGF = inn(107);
PC9_FGF = inn(108);
PC9_ERLOT = inn(109:111);
PC9_HRG = inn(112:115);

ERLOT_VEC = PC9_ERLOT(1)*signalConv(1,:) + PC9_ERLOT(2)*signalConv(2,:) + PC9_ERLOT(3)*signalConv(3,:);

PC9_data = zeros(16, 9);
PC9_data(2,:) = PC9_EGF(2)*signalConv(2,:) + PC9_EGF(1)*signalConv(1,:);
PC9_data(3,:) = PC9_HGF*signalConv(6,:);
PC9_data(4,:) = PC9_IGF*signalConv(5,:);
PC9_data(5,:) = PC9_PDGF*signalConv(7,:);
PC9_data(6,:) = PC9_HRG(1)*signalConv(1,:) + PC9_HRG(2)*signalConv(2,:) + PC9_HRG(3)*signalConv(3,:);
PC9_data(7,:) = PC9_FGF*signalConv(8,:);



PC9_data(9,:) = -ERLOT_VEC;
%PC9_data(10,:) = PC9_EGF(2)*signalConv(2,:) + PC9_EGF(1)*signalConv(1,:) - ERLOT_VEC;
PC9_data(11,:) = PC9_HGF*signalConv(6,:) - ERLOT_VEC;
PC9_data(12,:) = PC9_IGF*signalConv(5,:) - ERLOT_VEC;
PC9_data(13,:) = PC9_PDGF*signalConv(7,:) - ERLOT_VEC;
%PC9_data(14,:) = 
PC9_data(15,:) = PC9_FGF*signalConv(8,:) - ERLOT_VEC;
PC9_data(16,:) = -ERLOT_VEC;




PC9_data = PC9_data(1:16,[2 3 4 6 8 9]-1) - data((1:16)+46,[2 3 4 6 8 9]);

for ii = 1:16
    outE = outE + norm(PC9_data(ii,:));
end

%%%% End PC9 %%%%%





















%end

