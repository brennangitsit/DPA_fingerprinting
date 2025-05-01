function DPA_run_fingerprint_analysis_ROI(resultsDir, roiMapPath, tmapFile, taskNames, modalityPairs)
% DPA_RUN_FINGERPRINT_ANALYSIS_ROI
% Performs cross-modality fingerprinting using predefined ROIs

fprintf('Loading t-maps...\n');
load(tmapFile, 'asl_data', 'eng_data');
source = {asl_data, eng_data};
n_subjects = size(asl_data, 4);

fprintf('Loading ROI map...\n');
roiMap = niftiread(roiMapPath);
roiInfo = niftiinfo(roiMapPath);
roiLabels = unique(roiMap(:));
roiLabels(roiLabels == 0) = [];  % remove background
n_rois = length(roiLabels);

% Collect voxel indices for each ROI
roi_voxel_indices = cell(n_rois, 1);
for r = 1:n_rois
    roi_voxel_indices{r} = find(roiMap == roiLabels(r));
end

% Store results for all pairs
all_results = struct();

% Loop over modality pairs
for pairIdx = 1:length(modalityPairs)
    originModality = modalityPairs{pairIdx}{1};
    targetModality = modalityPairs{pairIdx}{2};

    origin_data = source{strcmp(taskNames, originModality)};
    target_data = source{strcmp(taskNames, targetModality)};
    target_data_2d = reshape(target_data, [], n_subjects);

    fprintf('Running fingerprinting: %s -> %s\n', originModality, targetModality);

    % Initialize results
    prd_id = nan(n_rois, n_subjects);
    prd_SI = nan(n_rois, n_subjects);
    prd_acc = zeros(n_rois, n_subjects);
    prd_acc_SI = zeros(n_rois, n_subjects);

    for subj = 1:n_subjects
        sub_origin = origin_data(:,:,:,subj);
        sub_origin_flat = reshape(sub_origin, [], 1);

        for roi = 1:n_rois
            ID = roi_voxel_indices{roi};
            if isempty(ID)
                continue;
            end

            t_sub_origin = sub_origin_flat(ID);
            t_target = target_data_2d(ID, :);

            if ~all(isnan(t_sub_origin)) && ~all(t_sub_origin == 0)
                corr_values = corr(t_sub_origin, t_target);
                [max_corr, predicted_id] = max(corr_values);
                prd_id(roi, subj) = predicted_id;
                prd_SI(roi, subj) = max_corr;
            end
        end
    end

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

    %% Save results per ROI
    mean_acc = mean(prd_acc, 2);
    mean_SI = mean(prd_acc_SI, 2);

    output = table(roiLabels, mean_acc, mean_SI, 'VariableNames', {'ROI', 'Accuracy', 'Similarity'});
    outTablePath = fullfile(resultsDir, [originModality '_' targetModality '_ROI_results.csv']);
    writetable(output, outTablePath);

    save(fullfile(resultsDir, ['sub_acc_' originModality '_' targetModality '_ROIs.mat']), 'prd_acc');

    all_results.(originModality).(targetModality) = output;
    fprintf('Finished: %s -> %s\n', originModality, targetModality);
end

%% ==================== CONJUNCTION AND REPORT ====================
if isfield(all_results, 'asl') && isfield(all_results.asl, 'eng') && ...
   isfield(all_results, 'eng') && isfield(all_results.eng, 'asl')

    disp('Creating conjunction report across directions...');
    asl2eng = all_results.asl.eng;
    eng2asl = all_results.eng.asl;

    assert(isequal(asl2eng.ROI, eng2asl.ROI), 'ROI mismatch between directions');

    conj_mask = (asl2eng.Accuracy > 0.5) & (eng2asl.Accuracy > 0.5);
    mean_acc = (asl2eng.Accuracy + eng2asl.Accuracy) / 2;

    conjunction_table = table(asl2eng.ROI, conj_mask, mean_acc, ...
        'VariableNames', {'ROI', 'IsSupramodal', 'MeanAccuracy'});

    writetable(conjunction_table, fullfile(resultsDir, 'ROI_conjunction_summary.csv'));

    n_supramodal = sum(conj_mask);
    disp(['Identified ' num2str(n_supramodal) ' supramodal ROIs.']);
    disp('Top supramodal regions:');
    [sorted_acc, sort_idx] = sort(mean_acc(conj_mask), 'descend');
    top_rois = asl2eng.ROI(conj_mask);
    for i = 1:min(5, length(sorted_acc))
        fprintf('  ROI %d: Mean Accuracy = %.3f\n', top_rois(sort_idx(i)), sorted_acc(i));
    end
end
end
