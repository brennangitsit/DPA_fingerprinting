% Code to inspect the content of sub_tmaps.mat
function inspect_sub_tmaps()
    % Load the file
    load('/Users/bterhunecotter/MyDrive/SDSU_LLCN/NEURO/Repositories/DPA_fingerprinting/sub_tmaps.mat');
    
    % Get variable names in the file
    var_names = who;
    disp('Variables in sub_tmaps.mat:');
    disp(var_names);
    
    % For each variable, display size and data type
    for i = 1:length(var_names)
        var_name = var_names{i};
        var_value = eval(var_name);
        var_size = size(var_value);
        var_type = class(var_value);
        
        fprintf('Variable: %s\n', var_name);
        fprintf('  Size: [%s]\n', strjoin(string(var_size), ' '));
        fprintf('  Type: %s\n', var_type);
        
        % If it's a 4D matrix, show details about subjects
        if length(var_size) == 4
            fprintf('  Number of subjects: %d\n', var_size(4));
            
            % Check if there are any NaN values
            nan_count = sum(isnan(var_value(:)));
            fprintf('  NaN values: %d\n', nan_count);
            
            % Show some basic statistics on the first subject
            first_subj = var_value(:,:,:,1);
            fprintf('  Stats for subject 1: min=%f, max=%f, mean=%f\n', ...
                min(first_subj(:)), max(first_subj(:)), mean(first_subj(~isnan(first_subj(:)))));
        end
        
        fprintf('\n');
    end
end