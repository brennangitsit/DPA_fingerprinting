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
% Path configuration - Change these paths to match your directory structure
projPath = "/Volumes/BTCruiser/GRFP"; % Base project directory
dataPath = fullfile(projPath, 'data');
resultsDir = fullfile(projPath, 'results');
toolboxPath = "/Users/bterhunecotter/MyDrive/SDSU_LLCN/NEURO/Repositories/NIFTI_toolbox";

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

% MODIFICATION: Handle individual subject masks
% Initialize empty mask for counting overlaps
commonMask = [];
maskDims = [];

% Loop through subjects to get mask dimensions from first subject
subID = subjectIDs{1};
subjMaskFile = fullfile(dataPath, ['sub-' subID], 'anat', ['sub-' subID '_mask.nii']);
subjMask = load_nii(subjMaskFile);
maskDims = size(subjMask.img);
commonMask = zeros(maskDims);
hdr = subjMask.hdr; % Save header for later use

% Create a common mask by finding voxels present in a threshold percentage of subjects
for s = 1:n_subjects
    subID = subjectIDs{s};
    subjMaskFile = fullfile(dataPath, ['sub-' subID], 'anat', ['sub-' subID '_mask.nii']);
    
    if exist(subjMaskFile, 'file')
        subjMask = load_nii(subjMaskFile);
        commonMask = commonMask + double(subjMask.img > 0);
    else
        warning('Mask file not found for subject %s', subID);
    end
end

% Create a common analysis mask (voxels present in at least 80% of subjects)
subjectThreshold = ceil(0.8 * n_subjects);
finalMask = commonMask >= subjectThreshold;

% Save the common mask for reference
common_mask_nii = make_nii(finalMask);
common_mask_nii.hdr = hdr;
save_nii(common_mask_nii, fullfile(resultsDir, 'common_analysis_mask.nii'));

disp(['Created common analysis mask with ' num2str(sum(finalMask(:))) ' voxels']);

%% ==================== LOAD CONTRAST MAPS ====================
% Load the t-maps (contrast maps) for ASL and English conditions
disp('Loading contrast maps...');

% Initialize data arrays
asl_data = nan([maskDims, n_subjects]); % ASL contrast t-maps
eng_data = nan([maskDims, n_subjects]); % English contrast t-maps

% Loop through subjects to load their data
for s = 1:n_subjects
    subID = subjectIDs{s};
    
    % Load ASL sentence > pseudosign contrast (contrast 14 in current iteration)
    asl_file = fullfile(dataPath, 'stats', '4fwhm', ['sub-' subID], 'spmT_0014.nii');
    if exist(asl_file, 'file')
        temp = load_nii(asl_file);
        asl_data(:,:,:,s) = temp.img;
    else
        warning('ASL contrast file not found for subject %s', subID);
    end
    
    % Load English sentence > pseudoword contrast (contrast 18 in current iteration)
    eng_file = fullfile(dataPath, 'stats', '4fwhm', ['sub-' subID], 'spmT_0018.nii');
    if exist(eng_file, 'file')
        temp = load_nii(eng_file);
        eng_data(:,:,:,s) = temp.img;
    else
        warning('English contrast file not found for subject %s', subID);
    end
end

% Save these matrices for future use
save(fullfile(dataPath, 'asl_eng_sent-pseudo_tmaps.mat'), 'asl_data', 'eng_data');
disp('T-maps saved');

%% ==================== DEFINE SEARCHLIGHTS ====================
% Generate searchlights using the common mask
disp('Defining searchlights...');

% Create searchlight list using the common mask
L = step1_defineSearchlight_volume(fullfile(resultsDir, 'common_analysis_mask.nii'), ...
                                   fullfile(resultsDir, 'common_analysis_mask.nii'), ...
                                   [searchlightRadius, targetVoxelCount]);

% Load searchlight information
load(fullfile(projPath, 'script', 'searchlight_list.mat'));
ID_set = L.LI;
voxel = L.voxel;

%% ==================== FINGERPRINTING ANALYSIS ====================
disp('Starting DPA fingerprinting analysis...');

% Create source/target data pairs
source = {asl_data, eng_data};

% Initialize results arrays to store identification results
n_voxels = prod(maskDims);

% Loop through modality pairs
for pairIdx = 1:length(modalityPairs)
    % Get modality pair information
    originModality = modalityPairs{pairIdx}{1};
    targetModality = modalityPairs{pairIdx}{2};
    
    originIdx = find(strcmp(taskname, originModality));
    targetIdx = find(strcmp(taskname, targetModality));
    
    origin_data = source{originIdx};
    target_data = source{targetIdx};
    
    disp(['Analyzing ' originModality ' -> ' targetModality ' identification...']);
    
    % Initialize arrays for this modality pair
    prd_id = nan(n_voxels, n_subjects);    % Predicted subject IDs
    prd_SI = nan(n_voxels, n_subjects);    % Similarity values
    prd_acc = zeros(n_voxels, n_subjects); % Accuracy (1=correct, 0=incorrect)
    prd_acc_SI = zeros(n_voxels, n_subjects); % Similarity values for correct IDs
    
    % Loop through subjects
    for subj = 1:n_subjects
        disp(['  Processing subject ' num2str(subj) ' of ' num2str(n_subjects)]);
        
        % Get origin data for this subject
        sub_origin = origin_data(:,:,:,subj);
        
        % Loop through searchlights
        for searchlight = 1:length(ID_set)
            if mod(searchlight, 1000) == 0
                disp(['    Searchlight ' num2str(searchlight) ' of ' num2str(length(ID_set))]);
            end
            
            % Get voxel indices for this searchlight
            ID = ID_set{searchlight};
            v_id = voxel(searchlight);
            
            % Initialize arrays for t-values in this searchlight
            t_sub_origin = nan(targetVoxelCount, 1);
            t_target = nan(targetVoxelCount, n_subjects);
            
            % Extract t-values for voxels in this searchlight
            for n = 1:size(ID, 2)
                a = ID(:, n);
                if a(1) <= size(sub_origin, 1) && a(2) <= size(sub_origin, 2) && a(3) <= size(sub_origin, 3)
                    t_sub_origin(n) = sub_origin(a(1), a(2), a(3));
                    
                    % Extract t-values for all subjects in target modality
                    for otherSubj = 1:n_subjects
                        t_target(n, otherSubj) = target_data(a(1), a(2), a(3), otherSubj);
                    end
                end
            end
            
            % Skip searchlights with no valid data
            if ~all(isnan(t_sub_origin)) && ~all(t_sub_origin == 0)
                % Correlate this subject's origin pattern with all subjects' target patterns
                corr_values = corr(t_sub_origin, t_target);
                
                % Find the subject with highest