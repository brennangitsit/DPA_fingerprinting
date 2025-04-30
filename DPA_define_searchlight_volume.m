function L = DPA_define_searchlight_volume(maskPath, sphere)
% DPA_define_searchlight_volume Generate a volumetric searchlight definition.
% 
% INPUTS:
%   maskPath - string path to a NIfTI binary mask file (.nii or .nii.gz)
%   sphere   - either [R] (radius in voxels) or [R, C] where C is target voxel count
%
% OUTPUT:
%   L - structure with fields:
%         LI:     n x 1 cell of linear indices for each searchlight
%         voxel:  n x 1 vector of center voxel indices
%         voxmin: n x 3 min voxel coordinate for each searchlight
%         voxmax: n x 3 max voxel coordinate for each searchlight

    % Read mask
    mask = niftiread(maskPath);
    mask = mask > 0;

    % Identify searchlight centers
    centerIdx = find(mask);
    [x, y, z] = ind2sub(size(mask), centerIdx);
    centers = [x'; y'; z'];

    % Voxels to consider
    [vx, vy, vz] = ind2sub(size(mask), find(mask));
    all_voxels = [vx'; vy'; vz'];

    % Searchlight configuration
    fixedRadius = isscalar(sphere);
    radius = sphere(1);
    if ~fixedRadius
        targetCount = sphere(2);
    end

    % Preallocate output
    nCenters = size(centers, 2);
    LI = cell(nCenters, 1);
    voxmin = zeros(nCenters, 3);
    voxmax = zeros(nCenters, 3);

    % Build searchlights
    for k = 1:nCenters
        dist = sqrt(sum((centers(:,k) - all_voxels).^2, 1));
        if fixedRadius
            voxel_subset = all_voxels(:, dist <= radius);
        else
            candidates = find(dist <= radius);
            [~, sorted_idx] = sort(dist(candidates));
            n_vox = min(targetCount, numel(sorted_idx));
            voxel_subset = all_voxels(:, candidates(sorted_idx(1:n_vox)));
        end

        coords = voxel_subset';  % Convert to N x 3
        LI{k} = sub2ind(size(mask), coords(:,1), coords(:,2), coords(:,3));
        voxmin(k,:) = min(voxel_subset, [], 2)';
        voxmax(k,:) = max(voxel_subset, [], 2)';
    end

    % Output structure
    L.LI = LI;
    L.voxel = centerIdx;
    L.voxmin = voxmin;
    L.voxmax = voxmax;
end
