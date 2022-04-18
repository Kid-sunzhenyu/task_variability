clear,clc
addpath(genpath('/mnt/Task_variability/cifti-matlab'));

taskn = 'REST1';

OrigPath = ['/mnt/HCP_DATASET/HCP_RESTING_STATE_FIX-Denoised(Compact)'];

SIDs = textscan('lists/list_82.txt','%s');
SIDs = importdata('lists/list_82.txt');
segNum = 2;

sess = {'LR', 'RL'};

for s = 1:length(SIDs)
    sub = SIDs(s);
    for i = 1:segNum
        seg = sess{i};
        disp([sub ':' seg])
        runname = ['rfMRI_' taskn '_' seg];
        name = [taskn '_' seg]
        InPath = [OrigPath '/' num2str(sub) '/' num2str(sub) '/MNINonLinear/Results/' runname];
        motion = ['/mnt/MotionInfo/' num2str(sub)]
        if exist([InPath '/' runname '_Atlas_hp2000_clean.dtseries.nii'])~=0
	       OutPath = ['/mnt/HCP_TASK_Output/Task_' taskn '/' num2str(sub) '/' runname];
           mkdir(OutPath)	
           cifti = cifti_read([InPath '/' runname '_Atlas_hp2000_clean.dtseries.nii']);        

           seg_data = cifti.cdata'; % TimePoints*nVertvex

            %% Filter
           seg_data_bpss = y_IdealFilter(seg_data,0.72,[0.01,0.08]);% 1200*91282
           GS = nanmean(seg_data_bpss,2);%1200*1
           GSdt = [0; GS(2:end)-GS(1:end-1)];
          
           %% Regression
           % Load Motion regressors
           Motion = importdata([motion '/' name '_Movement_Regressors.txt']);% 1200*12

           PreData = [];
           regressors = [Motion ones(size(GS,1),1)];% 1200*12
           beta = pinv(regressors)*seg_data_bpss;% 13*91282
           seg_data_bpss_resid = seg_data_bpss-(regressors(:,1:end-1)*beta(1:end-1,:));% 1200*91282
           cifti.cdata = seg_data_bpss_resid';
           ciftisave(cifti,[OutPath '/'  runname '_Atlas_hp2000_clean_bpss_12mr.dtseries.nii']);

           % Regression with GSR
           PreData_gsr = [];
           regressors = [Motion GS GSdt ones(size(GS,1),1)];% 1200*12
           beta = pinv(regressors)*seg_data_bpss;% 13*91282
           seg_data_bpss_resid_gsr = seg_data_bpss-(regressors(:,1:end-1)*beta(1:end-1,:));% 1200*91282
           cifti.cdata = seg_data_bpss_resid_gsr';
           cifti_write(cifti,[OutPath '/' runname '_Atlas_hp2000_clean_bpss_12mr_gsr.dtseries.nii']);
        end
    end
end