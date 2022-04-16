clear,clc;
addpath(genpath('/mnt/Task_variability/cifti-matlab')) % for cifti-matlab
addpath(genpath('lib')) % for my own codes

load fsLR_32k_config.mat
subs = importdata('lists/list_82.txt');
commons = importdata('lists/list_common_82.txt');
csubs = subs(commons==1);
nsubs = sum(commons);

taskn = 'REST1'; % 'REST2' 
mpath = ['/mnt/HCP_TASK_Output/Task_' taskn]; % change to ur path
rest_range = 1:300; % set REST frame range as 1:300 to mimic task length
sess = {'LR', 'RL'};
proc = '12mr';

Lhdr = gifti('rois/L.fslr_downsample_900mesh_parcellation_sm1.func.gii');
Lrois = Lhdr.cdata;
maxLrois = max(Lrois);
Rhdr = gifti('rois/R.fslr_downsample_900mesh_parcellation_sm1.func.gii');
Rrois = Rhdr.cdata;
maxRrois = max(Rrois);

OutPath = [mpath '/results/' taskn];
mkdir(OutPath)

nLR = 59412; % L = 29696; R = 29716

nsess = 2;
IntraVariance = zeros(nsubs, nLR);
for s = 1:nsubs %length(subs)
    sub = num2str(csubs(s));
    Rmat = zeros(nLR, 1483, nsess);
    for i = 1:nsess
	    i
        tasksess = ['rfMRI_' taskn '_' sess{i}];
        DataPath = [mpath '/' sub '/' tasksess];
        filename = [DataPath '/' tasksess '_Atlas_hp2000_clean_bpss_' proc '.dtseries.nii'];
        
        chdr = cifti_read(filename);
        alldata = single(chdr.cdata(:,rest_range));    % rest range    
        Lsurfdata = single(zeros(32492,size(alldata,2)));
        Lsurfdata(Lvertlist,:) = alldata([Lstart:Lcount],:);
        Rsurfdata = single(zeros(32492,size(alldata,2)));
        Rsurfdata(Rvertlist,:) = alldata([Rstart:Rstart+Rcount-1],:);
        
        Lsurfdata = single(zeros(32492,size(alldata,2)));
        Lsurfdata(Lvertlist,:) = alldata(1:29696,:);
        Rsurfdata = single(zeros(32492,size(alldata,2)));
        Rsurfdata(Rvertlist,:) = alldata([Rstart:Rstart+Rcount-1],:);
        
        data_rois_all = single(zeros(int32(maxRrois),size(alldata,2)));
        for n = 1:maxLrois
            data_rois_all(n,:) = nanmean(Lsurfdata(Lrois == n,:));
        end
        
        for n = maxLrois+1:maxRrois
            data_rois_all(n,:) = nanmean(Rsurfdata(Rrois == n,:));
        end
        
        Rmat(:,:,i) = single(corr(alldata(1:nLR,:)', data_rois_all'));
    end
    Rmat(isnan(Rmat)) = 0;
    count = 0;
    AveRmat = zeros(nLR, 1);
    clear alldata
    for m = 1:nsess
        for n = m+1:nsess
            count = count + 1;
            tmp = my_matcorr(squeeze(Rmat(:,:,m))', squeeze(Rmat(:,:,n))');
            tmp(isnan(tmp)) = 0;
            AveRmat = AveRmat + tmp';%diag(tmp);
        end
    end
    count;
    IntraVariance(s, :) = 1 - AveRmat/count;
    %save_mgh(IntraVariance(s,:), [OutPath '/lh.' sub '_intravariance_fs4.mgh'],eye(4))
end
meanIntraVariance = mean(IntraVariance);
save([OutPath '/IntraVariance_' taskn '_' proc '_LR.mat'],'meanIntraVariance')

filename = ['IntraVariance_' taskn '_' proc];
Func_write_func_gifti_32k(filename, meanIntraVariance, OutPath, Lhdr, Rhdr)
