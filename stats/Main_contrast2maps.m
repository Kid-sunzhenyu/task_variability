clear,clc;
addpath('../variability')
% for contrast the task variability against the rest variability

mpath = '/zfs/musc/david/HCP4variability';
outpath = [mpath '/contrasts'];
mkdir(outpath)
proc = '12mr_gsr';
nLR = 59412;

load([tpath '/AllIndividuals_InterSubject_Variability_REST2_' proc '_LR.mat'])
rindimap1 = indimap;
load([tpath '/AllIndividuals_InterSubject_Variability_REST1_' proc '_LR.mat'])
rindimap2 = indimap;

[~,p,~,stat] = ttest(rindimap1, rindimap2);
tmap = stat.tstat;
filename = ['contrastmap_ttest_REST2_REST1_' proc];
Func_write_func_gifti_32k(filename, tmap, outpath, Lhdr, Rhdr)

tasks = {'LANGUAGE', 'WM', 'SOCIAL', 'MOTOR', 'GAMBLING', 'EMOTION', 'SOCIAL'};
for t = 1:7
    tpath = [mpath '/' tasks{t}];
    load([tpath '/AllIndividuals_InterSubject_Variability_' tasks{t} '_' proc '_LR.mat'])
    tindimap = indimap;
    [~,p,~,stat] = ttest(rindimap1, rindimap2);
    tmap = stat.tstat;
    filename = ['contrastmap_ttest_' tasks{t} '_REST2_' proc];
    Func_write_func_gifti_32k(filename, tmap, outpath, Lhdr, Rhdr)
end