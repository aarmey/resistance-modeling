%
clc; clear;
import bioma.data.*;

load DMobj;

RTKnames = {'MET', 'FGFR1', 'IGF1R', 'PDGFRB', 'PDGFRA', 'ERBB3', 'EGFR'};

CellLines = {'EFM19_BREAST','SKBR3_BREAST', 'BT474_BREAST', 'AU565_BREAST', 'HCC1954_BREAST', ...
    'CHL1_SKIN','FADU_UPPER_AERODIGESTIVE_TRACT','PC14_LUNG',...
    'HCC827_LUNG', 'HCC4006_LUNG', 'CALU3_LUNG','NCIH358_LUNG','NCIH1648_LUNG',...
    'NCIH1666_LUNG', 'EBC1_LUNG', 'NCIH1703_LUNG', 'A204_SOFT_TISSUE',...
    'RT112_URINARY_TRACT','RT4_URINARY_TRACT','KELLY_AUTONOMIC_GANGLIA',...
    'NCIH2228_LUNG','A375_SKIN','MALME3M_SKIN','COLO201_LARGE_INTESTINE'};

MET = [0 0 0 1 1 2 2 2 2 2 1 0 0 2 -1 2 0 2 0 0 0 1 1 1];
IGF1R = [0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

subplot(1,2,1);
scatter(MET(MET > -1), DMobj(CellLines(MET > -1),'MET'), 30, 'k','filled');

ylabel('MET probe value');
set(gca, 'FontName', 'Myriad Pro')
set(gca,'FontSize',16)
axis square;
axis([-0.5 2.5 4 12]);


subplot(1,2,2);
scatter(IGF1R(IGF1R > -1), DMobj(CellLines(IGF1R > -1),'IGF1R'), 30, 'k','filled');

ylabel('IGF1R probe value');
set(gca, 'FontName', 'Myriad Pro')
set(gca,'FontSize',16)
axis square;
axis([-0.5 2.5 4 12]);

print2eps('untitledRecp');