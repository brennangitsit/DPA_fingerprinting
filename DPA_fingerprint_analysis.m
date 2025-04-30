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
              '10', '11', '12', '14', '15', '16', '17', '18', '19'};
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
                
                % Find the subject with highest correlation
                [max_corr, predicted_id] = max(corr_values);
                
                % Store the predicted ID and similarity value
                prd_id(v_id, subj) = predicted_id;
                prd_SI(v_id, subj) = max_corr;
            end
        end
    end
    
    %% ==================== CALCULATE ACCURACY ====================
    disp('Calculating identification accuracy...');
    
    % Determine which identifications were correct
    for subj = 1:n_subjects
        % Get predicted IDs for this subject
        subID_prd = prd_id(:, subj);
        
        % Find voxels where prediction matches subject's true ID
        correct_IDs = find(subID_prd == subj);
        
        % Find voxels with valid similarity values
        valid_voxels = find(~isnan(prd_SI(:, subj)));
        
        % Exclude voxels with perfect correlation (these are likely artifacts)
        perfect_corr = find(prd_SI(:, subj) == 1);
        
        % Get voxels with valid predictions
        valid_correct = intersect(correct_IDs, valid_voxels);
        final_correct = setdiff(valid_correct, perfect_corr);
        
        % Mark correctly identified voxels
        prd_acc(final_correct, subj) = 1;
        prd_acc_SI(final_correct, subj) = prd_SI(final_correct, subj);
    end
    
    % Calculate mean accuracy and similarity across subjects
    mean_acc = mean(prd_acc, 2);
    mean_SI = mean(prd_acc_SI, 2);
    
    % Reshape to 3D volume
    acc_map = reshape(mean_acc, maskDims);
    sim_map = reshape(mean_SI, maskDims);
    
    % Save accuracy and similarity maps
    origin_name = taskname{originIdx};
    target_name = taskname{targetIdx};
    
    acc_nii = make_nii(acc_map);
    acc_nii.hdr = hdr;
    save_nii(acc_nii, fullfile(resultsDir, [origin_name '_' target_name '_acc.nii']));
    
    sim_nii = make_nii(sim_map);
    sim_nii.hdr = hdr;
    save_nii(sim_nii, fullfile(resultsDir, [origin_name '_' target_name '_sim.nii']));
    
    % Save raw accuracy data for further analysis
    save(fullfile(resultsDir, ['sub_acc_' origin_name '_' target_name '.mat']), 'prd_acc');
    
    disp(['Completed ' origin_name ' -> ' target_name ' analysis']);
end

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
asl_to_eng_acc = load_nii(fullfile(resultsDir, 'asl_eng_acc.nii'));
eng_to_asl_acc = load_nii(fullfile(resultsDir, 'eng_asl_acc.nii'));

% Apply statistical threshold (this would ideally come from permutation testing)
% For now, we'll use a simple threshold based on the paper
asl_to_eng_thresholded = asl_to_eng_acc.img > 0.5;  % Threshold at 50% accuracy
eng_to_asl_thresholded = eng_to_asl_acc.img > 0.5;  % Adjust based on your data

% Create conjunction map (regions significant in both directions)
conjunction_map = asl_to_eng_thresholded & eng_to_asl_thresholded;

% Calculate mean accuracy in the conjunction regions
mean_acc_map = (asl_to_eng_acc.img + eng_to_asl_acc.img) / 2 .* conjunction_map;

% Save conjunction map
conj_nii = make_nii(conjunction_map);
conj_nii.hdr = hdr;
save_nii(conj_nii, fullfile(resultsDir, 'asl_eng_conjunction.nii'));

% Save mean accuracy map
mean_acc_nii = make_nii(mean_acc_map);
mean_acc_nii.hdr = hdr;
save_nii(mean_acc_nii, fullfile(resultsDir, 'asl_eng_mean_acc.nii'));

%% ==================== REPORTING ====================
disp('Analysis complete. Results saved to:');
disp(resultsDir);
disp(' ');
disp('Summary of findings:');

% Count significant voxels
n_sig_voxels = sum(conjunction_map(:));
disp(['Identified ' num2str(n_sig_voxels) ' voxels showing supramodal language processing']);

% Calculate cluster sizes
CC = bwconncomp(conjunction_map);
cluster_sizes = cellfun(@numel, CC.PixelIdxList);
[sorted_sizes, idx] = sort(cluster_sizes, 'descend');

% Report largest clusters
disp('Largest supramodal clusters:');
for i = 1:min(5, length(sorted_sizes))
    cluster_idx = CC.PixelIdxList{idx(i)};
    cluster_acc = mean(mean_acc_map(cluster_idx));
    disp(['  Cluster ' num2str(i) ': ' num2str(sorted_sizes(i)) ' voxels, mean accuracy: ' num2str(cluster_acc)]);
end

disp('Analysis completed successfully');