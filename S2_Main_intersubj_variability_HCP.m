clear,clc;
addpath(genpath('/zfs/musc/david/codes/tools/'))
addpath(genpath('/zfs/musc/david/codes/matlab'))
load fsLR_32k_config.mat

subs = importdata('list_100Unrelated.txt');
commons = importdata('list_common_subs.txt');
csubs = subs(commons==1);
nsubs = sum(commons);

mpath = '/zfs/musc/david/HCP4variability'; % change to ur path
taskn = 'WM'; % 'LANGUAGE' % 'REST1'
sess = {'LR', 'RL'};
proc = '12mr_gsr';

Lhdr = gifti('../rois/L.fslr_downsample_900mesh_parcellation_sm1.func.gii');
Lrois = Lhdr.cdata;
maxLrois = max(Lrois);
Rhdr = gifti('../rois/R.fslr_downsample_900mesh_parcellation_sm1.func.gii');
Rrois = Rhdr.cdata;
maxRrois = max(Rrois);

OutPath = [mpath '/results/' taskn];
mkdir(OutPath)

nLR = 59412; % L = 29696; R = 29716
nsess = 2;
InterVariance = zeros(nsess, nLR);
for i = 1:nsess %length(subs)
    Rmat = zeros(nLR, 1483, nsubs);
    for s = 1:nsubs
	    s
        sub = num2str(csubs(s));
        tasksess = ['tfMRI_' taskn '_' sess{i}];
        DataPath = [mpath '/data/' sub '/' tasksess];
        filename = [DataPath '/' tasksess '_Atlas_hp200_s4_bpss_' proc '.dtseries.nii'];
        
        chdr = cifti_read(filename);
        alldata = single(chdr.cdata);       
        Lsurfdata = single(zeros(32492,size(alldata,2)));
        Lsurfdata(Lvertlist,:) = alldata([Lstart:Lcount],:);
        Rsurfdata = single(zeros(32492,size(alldata,2)));
        Rsurfdata(Rvertlist,:) = alldata([Rstart:Rstart+Rcount-1],:);
        
        data_rois_all = single(zeros(int32(maxRrois),size(alldata,2)));
        for n = 1:maxLrois
            data_rois_all(n,:) = nanmean(Lsurfdata(Lrois == n,:));
        end
        
        for n = maxLrois+1:maxRrois
            data_rois_all(n,:) = nanmean(Rsurfdata(Rrois == n,:));
        end
        
        Rmat(:,:,s) = single(corr(alldata(1:nLR,:)', data_rois_all'));
    end
    Rmat(isnan(Rmat)) = 0;
    count = 0;
    AveRmat = zeros(nsubs,nsubs,nLR);
    clear alldata
    for m = 1:nsubs %length(subs)
        for n = 1:nsubs %length(subs)
            count = count + 1;
            tmp = my_matcorr(squeeze(Rmat(:,:,m))', squeeze(Rmat(:,:,n))');
            tmp(isnan(tmp)) = 0;
            AveRmat = AveRmat + tmp';
        end
    end
    count
    InterVariance(i, :) = 1-AveRmat/count;
    %save_mgh(IntraVariance(s,:), [OutPath '/lh.' sub '_intravariance_fs4.mgh'],eye(4))
end
meanInterVariance = mean(InterVariance);
save([OutPath '/InterVariance_' taskn '_' proc '_LR.mat'], 'InterVariance', 'meanInterVariance')

filename = ['InterVariance_' taskn '_' proc];
Func_write_func_gifti_32k(filename, meanInterVariance, OutPath, Lhdr, Rhdr)

% Regress Out
Variability_norm = zeros(nsess,nLR);
Variability = zeros(nsess,nLR);

load([OutPath '/IntraVariance_' taskn '_' proc '_LR.mat'])
for i = 1:nsess
    Intra_value = meanIntraVariance;
    Inter = InterVariance(i,:)';
    X = [Intra_value', ones(nLR,1)];
    beta = pinv(X)*Inter;
    tmp = Inter - X*beta;
    Variability_norm(i,:) = tmp;
    
    tmp = Inter - X(:,1:end-1)*beta(1:end-1);
    Variability(i,:) = tmp;
end
save([OutPath '/InterSubject_Variability_session_' taskn '_' proc '_LR.mat'], 'Variability_norm', 'Variability')

% normalized variability
meanVariability = mean(Variability_norm);
filename = ['InterSubject_Variability_norm_' taskn '_' proc];
Func_write_func_gifti_32k(filename, meanVariability, OutPath, Lhdr, Rhdr)


% variability
meanVariability = mean(Variability);
filename = ['InterSubject_Variability_norm_' taskn '_' proc];
Func_write_func_gifti_32k(filename, meanVariability, OutPath, Lhdr, Rhdr)

% meanVariability = mean(Variability);
% Lhdr.cdata(Lvertlist) = meanVariability(1:29696);
% save(Lhdr, [OutPath '/InterSubject_Variability_' taskn '_' proc '_L.func.gii'])
% Rhdr.cdata = 0*Rhdr.cdata;
% Rhdr.cdata(Rvertlist) = meanVariability(Rstart:end)';
% save(Rhdr, [OutPath '/InterSubject_Variability_' taskn '_' proc '_R.func.gii'])