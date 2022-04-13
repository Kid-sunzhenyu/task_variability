clear,clc
%addpath(genpath('/data/user/limeiling/SoftWare/WorkBench/CIFTIMatlabReaderWriter'));
%https://github.com/Washington-University/cifti-matlab.git
taskn = 'MOTOR';

OrigPath = ['/zfs/musc/Weigang/HCP_Data/HCP_3T_' taskn '/Preprocessing/'];

SIDs = textread('list_100Unrelated.txt','%s');

segNum = 2;

sess = {'LR', 'RL'};

for s = 1:1 %length(SIDs)
    sub = SIDs{s};
    for i = 1:segNum
        seg = sess{i};
        
        runname = ['tfMRI_' taskn '_' seg];
        
        InPath = [OrigPath sub '/MNINonLinear/Results/' runname];
        OutPath = ['/zfs/musc/david/HCP4variability/data/' sub '/' runname];
        mkdir(OutPath)
  
        cifti = open_wbfile([InPath '/' runname '_Atlas_hp200_s4.dtseries.nii']);        

       seg_data = cifti.cdata'; % TimePoints*nVertvex

        %% Filter
       seg_data_bpss = y_IdealFilter(seg_data,0.72,[0.01,0.08]);% 1200*91282
       GS = nanmean(seg_data_bpss,2);%1200*1
       GSdt = [0, GS(2:end)-GS(1:end-1)];

       %% Regression
       % Load Motion regressors
       Motion = importdata([InPath '/Movement_Regressors.txt']);% 1200*12

       PreData = [];
       regressors = [Motion ones(size(GS),1)];% 1200*12
       beta = pinv(regressors)*seg_data_bpss;% 13*91282
       seg_data_bpss_resid = seg_data_bpss-(regressors(:,1:end-1)*beta(1:end-1,:));% 1200*91282
       cifti.cdata = seg_data_bpss_resid';
       ciftisave(cifti,[OutPath '/' sub '/' runname '_Atlas_hp200_s4_bpss_12mr.dtseries.nii']);

       % Regression with GSR
       PreData_gsr = [];
       regressors = [Motion GS GSdt ones(size(GS,1))];% 1200*12
       beta = pinv(regressors)*seg_data_bpss;% 13*91282
       seg_data_bpss_resid_gsr = seg_data_bpss-(regressors(:,1:end-1)*beta(1:end-1,:));% 1200*91282
       cifti.cdata = seg_data_bpss_resid_gsr';
       ciftisave(cifti,[OutPath '/' sub '/' runname '__LR_Atlas_hp200_s4_bpss_12mr_gsr.dtseries.nii']);

    end
end
