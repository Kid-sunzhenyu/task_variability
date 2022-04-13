#! /bin/csh

set InPath =  #Change this
set OutPath = /autofs/space/bidlin9_001/users/HCP_FIX/DataProcessed_FIX_s4/Q3
set Downsamplefolder = /autofs/space/bidlin9_001/users/HCP_FIX/DataProcessed_FIX_s4/Q3
set TemplatePath = autofs/space/bidlin9_002/users/Jianxun/codes/tool/workbench/Pipelines/global/templates/standard_mesh_atlases/resample_fsaverage/
#set subjects_dir = /usr/local/freesurfer/subjects
set count = 56
set stop = 67 #change

set FWHM = 4
set sigma = 1.5 #FWHM = 4, sigma = FHWM/2.355

set att_file = list_100Unrelated.txt

while($count <= $stop)
    set sub = `head -n $count $att_file | tail -n 1 | awk '{print $1}'`
    echo "${count}:${sub}"

set cmd = "wb_command -cifti-smoothing ${InPath}/${sub}/rfMRI_REST_Atlas_hp2000_clean_bpss_mr_gsr.dtseries.nii ${sigma} ${sigma} COLUMN ${InPath}/${sub}/rfMRI_REST_Atlas_hp2000_clean_bpss_mr_gsr_s${FWHM}.dtseries.nii -left-surface $Downsamplefolder/${sub}/MNINonLinear/fsaverage_LR32k/${sub}.L.midthickness.32k_fs_LR.surf.gii -right-surface $Downsamplefolder/${sub}/MNINonLinear/fsaverage_LR32k/${sub}.R.midthickness.32k_fs_LR.surf.gii"

	pbsubmit -n 4  -c "$cmd"
  @ count = $count + 1
end



