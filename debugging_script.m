% SCRIPT FOR DEBUGGING

disp('===DEBUGGING===')

load('sub_acc_asl_eng.mat', 'prd_acc');  % [n_voxels x n_subjects]
acc_asl_eng = prd_acc;

load('sub_acc_eng_asl.mat', 'prd_acc');
acc_eng_asl = prd_acc;

load('asl_eng_conjunction.nii.gz');
conj_mask = niftiread('asl_eng_conjunction.nii.gz') > 0;
conj_voxels = find(conj_mask(:));

% Subset to supramodal regions only
acc_conj_asl_eng = acc_asl_eng(conj_voxels, :);
acc_conj_eng_asl = acc_eng_asl(conj_voxels, :);

% Average across directions (optional)
mean_acc_per_subject = (sum(acc_conj_asl_eng, 1) + sum(acc_conj_eng_asl, 1)) / 2;
