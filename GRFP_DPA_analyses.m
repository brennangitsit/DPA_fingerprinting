% DPA_fingerprint_analysis.m
% This script performs Distributed Pattern Activity (DPA) fingerprinting analysis
% to identify brain regions where activation patterns are:
%   1. Similar across language modalities (ASL and English)
%   2. Unique to each individual
%   3. Stable within individuals across modalities
% 
% Adapted from Liu et al. (2020) NeuroImage paper on supramodal language networks

clear;
close all;

%% ==================== CONFIGURATION SECTION ====================
% Path configuration 
projPath = "/Users/bterhunecotter/_GRFP_data_local"; % Base project directory
dataPath = fullfile(projPath, 'stats/0fwhm');
resultsDir = fullfile(projPath, 'results/dpa');
toolboxPath = "/Users/bterhunecotter/_Repositories/NIFTI_toolbox";

% Create results directory if it doesn't exist
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Add NIFTI toolbox to path
addpath(genpath(toolboxPath));

% Load data parameters
subjectIDs = {'01', '02', '03', '04', '05', '06', '07', '08', '09', ...
              '10', '11', '12', '14', '15', '16', '17', '18', '19'}; % skip 20 because of anat overlap
n_subjects = length(subjectIDs);

% Analysis parameters
searchlightRadius = 5; % Radius in voxels (~15mm)
targetVoxelCount = 125; % Number of voxels in each searchlight sphere
p_threshold = 0.005; % Statistical threshold (FDR corrected)
min_cluster_size = 30; % Minimum cluster size

% Define modalities
taskname = {'asl', 'eng'};
modalityPairs = {{'asl', 'eng'}, {'eng', 'asl'}}; % Analyze both directions

%% ==================== DATA PREPARATION ====================
% Load subject masks and create a common analysis mask
disp('Loading subject masks and creating common analysis mask...');

% Initialize variables
commonMask = [];
maskDims = [];

% Loop through subjects to get mask dimensions from the first subject
subID = subjectIDs{1};
subjMaskFile = fullfile(dataPath, ['sub-' subID], 'mask.nii');
subjMask = niftiread(subjMaskFile);
info = niftiinfo(subjMaskFile);
maskDims = size(subjMask);
commonMask = zeros(maskDims);

% Create a common mask by finding voxels present in a threshold percentage of subjects
for s = 1:n_subjects
    subID = subjectIDs{s};
    subjMaskFile = fullfile(dataPath, ['sub-' subID], 'mask.nii');
    
    if exist(subjMaskFile, 'file')
        subjMask = niftiread(subjMaskFile);
        commonMask = commonMask + double(subjMask > 0);
    else
        warning('Mask file not found for subject %s', subID);
    end
end

% Create a common analysis mask (voxels present in at least 80% of subjects)
subjectThreshold = ceil(0.8 * n_subjects);
finalMask = commonMask >= subjectThreshold;

% Save the common mask using NIfTI tools
info.Datatype = 'uint8';  % Ensure correct datatype
niftiwrite(uint8(finalMask), fullfile(resultsDir, 'common_analysis_mask.nii'), info);

disp(['Created common analysis mask with ' num2str(sum(finalMask(:))) ' voxels']);

%% ==================== LOAD CONTRAST MAPS ====================
% Load the t-maps (contrast maps) for ASL and English conditions
disp('Loading contrast maps...');

% Initialize data arrays
asl_data = nan([maskDims, n_subjects]);  % ASL contrast t-maps
eng_data = nan([maskDims, n_subjects]);  % English contrast t-maps

% Loop through subjects to load their data
for s = 1:n_subjects
    subID = subjectIDs{s};
    
    % Load ASL sentence > pseudosign contrast (contrast 14)
    asl_file = fullfile(dataPath, ['sub-' subID], 'spmT_0014.nii');
    if exist(asl_file, 'file')
        temp = niftiread(asl_file);
        asl_data(:,:,:,s) = temp;
    else
        warning('ASL contrast file not found for subject %s', subID);
    end

    % Load English sentence > pseudoword contrast (contrast 18)
    eng_file = fullfile(dataPath, ['sub-' subID], 'spmT_0018.nii');
    if exist(eng_file, 'file')
        temp = niftiread(eng_file);
        eng_data(:,:,:,s) = temp;
    else
        warning('English contrast file not found for subject %s', subID);
    end
end

% Save the data arrays for future use
save(fullfile(resultsDir, 'asl_eng_sent-pseudo_tmaps.mat'), 'asl_data', 'eng_data');
disp('T-maps saved.');

%% ==================== DEFINE SEARCHLIGHTS ====================
% Generate searchlights using the common mask
disp('Defining searchlights...');

commonMaskPath = fullfile(resultsDir, 'common_analysis_mask.nii');

% Create searchlight list and get voxel indices
L = DPA_define_searchlight_volume(commonMaskPath, ...
                                   [searchlightRadius, targetVoxelCount]);

ID_set = L.LI;
voxel = L.voxel;

% Save searchlight list for future use
save(fullfile(projPath, 'results', 'searchlight_list.mat'), 'L');

%% ==================== SEARCHLIGHT ANALYSIS ====================
disp('Starting DPA searchlight analysis...');

DPA_run_fingerprint_analysis_SL( ...
    resultsDir, ...
    fullfile(resultsDir, 'common_analysis_mask.nii'), ...
    fullfile(resultsDir, 'searchlight_list.mat'), ...
    fullfile(resultsDir, 'asl_eng_sent-pseudo_tmaps.mat'), ...
    {'asl', 'eng'}, ...
    modalityPairs ...
);

%% ==================== PERMUTATION TESTING ====================
disp('Performing permutation testing to assess statistical significance...');

% This section would implement permutation testing
% For each voxel, randomly reassign subject identities and repeat the identification
% process to create a null distribution
% Compare actual identification rates to this null distribution
% Implement if resources permit (this is computationally intensive)

%% ==================== CREATE CONJUNCTION MAP ====================
disp('Creating conjunction map of significant regions...');

% Load accuracy maps
asl_acc_path = fullfile(resultsDir, 'asl_eng_acc.nii.gz');
eng_acc_path = fullfile(resultsDir, 'eng_asl_acc.nii.gz');

asl_acc_img = niftiread(asl_acc_path);
eng_acc_img = niftiread(eng_acc_path);
mask_info = niftiinfo(asl_acc_path);  % Use header from one of the maps

% Apply threshold (default 50% accuracy, adjust as needed)
asl_thresholded = asl_acc_img > 0.5;
eng_thresholded = eng_acc_img > 0.5;

% Create conjunction map
conjunction_map = asl_thresholded & eng_thresholded;

% Compute mean accuracy within conjunction
mean_acc_map = (asl_acc_img + eng_acc_img) / 2;
mean_acc_map(~conjunction_map) = 0;  % Zero outside conjunction

% Save conjunction map
mask_info.Datatype = 'uint8';
niftiwrite(uint8(conjunction_map), ...
    fullfile(resultsDir, 'asl_eng_conjunction.nii.gz'), ...
    mask_info, 'Compressed', true);

% Save mean accuracy map
mask_info.Datatype = 'single';
niftiwrite(single(mean_acc_map), ...
    fullfile(resultsDir, 'asl_eng_mean_acc.nii.gz'), ...
    mask_info, 'Compressed', true);

disp('Conjunction and mean accuracy maps saved successfully.');


%% ==================== REPORTING ====================
disp(' ');
disp('==================== REPORTING ====================');
disp(['Analysis complete. Results saved to: ' resultsDir]);
disp(' ');

% Count significant voxels
n_sig_voxels = nnz(conjunction_map);
disp(['Identified ' num2str(n_sig_voxels) ' voxels showing supramodal language processing.']);

% Calculate connected clusters in the conjunction map
CC = bwconncomp(conjunction_map, 6);  % 6-connectivity (default for 3D)
cluster_sizes = cellfun(@numel, CC.PixelIdxList);
[sorted_sizes, sorted_idx] = sort(cluster_sizes, 'descend');

% Report the largest clusters
n_clusters_to_show = min(5, numel(sorted_sizes));
disp(['Found ' num2str(CC.NumObjects) ' clusters. Showing top ' num2str(n_clusters_to_show) ':']);
for i = 1:n_clusters_to_show
    cluster_voxels = CC.PixelIdxList{sorted_idx(i)};
    cluster_mean_acc = mean(mean_acc_map(cluster_voxels));
    fprintf('  Cluster %d: %d voxels, mean accuracy = %.3f\n', ...
        i, sorted_sizes(i), cluster_mean_acc);
end

disp('===================================================');

%% ==================== ROI ANALYSIS ===========================
disp('Starting DPA ROI analysis...');

roiDir = "/Users/bterhunecotter/_NEURO/ROIs";
roiMapPath = fullfile(roiDir, 'parcelsASL_homologous_r53x65x56.nii');

DPA_run_fingerprint_analysis_ROI( ...
    resultsDir, ...
    roiMapPath, ... % ROI image
    fullfile(resultsDir, 'asl_eng_sent-pseudo_tmaps.mat'), ...
    {'asl', 'eng'}, ...
    modalityPairs ...
);

disp('ROI analysis complete.');