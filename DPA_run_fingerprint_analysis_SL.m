function DPA_run_fingerprint_analysis_SL(resultsDir, commonMaskPath, searchlightListPath, tmapFile, taskNames, modalityPairs)
% RUN_DPA_FINGERPRINT_ANALYSIS
% Performs cross-modality individual identification (fingerprinting)
% using searchlight-based correlations between t-maps.

% Load subject t-maps (4D: x, y, z, subject)
fprintf('Loading t-maps...\n');
load(tmapFile, 'asl_data', 'eng_data');
source = {asl_data, eng_data};
n_subjects = size(asl_data, 4);

% Load the common mask metadata for output NIfTI files
mask_info = niftiinfo(commonMaskPath);
maskDims = mask_info.ImageSize;

% Load searchlight list
fprintf('Loading searchlight list...\n');
load(searchlightListPath, 'L');
ID_set = L.LI;
voxel = L.voxel;
n_voxels = prod(maskDims);

% Loop over modality pairs (e.g., ASL -> ENG)
for pairIdx = 1:length(modalityPairs)
    originModality = modalityPairs{pairIdx}{1};
    targetModality = modalityPairs{pairIdx}{2};

    origin_data = source{strcmp(taskNames, originModality)};
    target_data = source{strcmp(taskNames, targetModality)};
    target_data_2d = reshape(target_data, [], n_subjects);  % reshape once before loops

    fprintf('Running fingerprinting: %s -> %s\n', originModality, targetModality);

    % Initialize results
    prd_id = nan(n_voxels, n_subjects);
    prd_SI = nan(n_voxels, n_subjects);
    prd_acc = zeros(n_voxels, n_subjects);
    prd_acc_SI = zeros(n_voxels, n_subjects);

    % Loop through subjects
    for subj = 1:n_subjects
        sub_origin = origin_data(:,:,:,subj);

        % Loop through searchlights
        for sl = 1:length(ID_set)
            ID = ID_set{sl};
            v_id = voxel(sl);

            % Extract linear indices (1D)
            if isempty(ID)
                continue;
            end

            % Flatten indices
            t_sub_origin = sub_origin(ID);
            t_target = target_data_2d(ID, :);

            % Check for valid pattern
            if ~all(isnan(t_sub_origin)) && ~all(t_sub_origin == 0)
                corr_values = corr(t_sub_origin, t_target);
                [max_corr, predicted_id] = max(corr_values);
                prd_id(v_id, subj) = predicted_id;
                prd_SI(v_id, subj) = max_corr;
            end
        end
    end

    %% Compute accuracy and similarity
    fprintf('Calculating accuracy...\n');
    for subj = 1:n_subjects
        predicted = prd_id(:, subj);
        correct_idx = find(predicted == subj);
        valid_idx = find(~isnan(prd_SI(:, subj)));
        perfect_corr = find(prd_SI(:, subj) == 1);
        final_idx = setdiff(intersect(correct_idx, valid_idx), perfect_corr);
        prd_acc(final_idx, subj) = 1;
        prd_acc_SI(final_idx, subj) = prd_SI(final_idx, subj);
    end

    %% Save group maps
    mean_acc = mean(prd_acc, 2);
    mean_SI = mean(prd_acc_SI, 2);

    acc_nii = reshape(mean_acc, maskDims);
    sim_nii = reshape(mean_SI, maskDims);

    mask_info.Datatype = 'single';

    outPrefix = fullfile(resultsDir, [originModality '_' targetModality]);
    niftiwrite(single(acc_nii), outPrefix + "_acc.nii", mask_info, 'Compressed', true);
    niftiwrite(single(sim_nii), outPrefix + "_sim.nii", mask_info, 'Compressed', true);

    % Save raw subject accuracy
    save(resultsDir + "/sub_acc_" + originModality + "_" + targetModality + ".mat", 'prd_acc');

    fprintf('Finished: %s -> %s\n', originModality, targetModality);
end
end
