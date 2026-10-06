% Create table with correct responses to each stimulus for each participant
% as well as the mapping of stimuli to learning contexts (which was both
% counterbalanced between participants) for later use when simulating data
% for parameter recovery

correct_response_table = stim_rew_prob;
correct_response_table.id = cellstr(correct_response_table.id); % convert id column to string
correct_response_table.id = extractBefore(correct_response_table.id, 7); % only keep first 6 characters
ids = unique(correct_response_table.id); % create list of ids for loop % create list of ids for use in loops

% 1. Stimulus to learning contexts mapping
learning_context_mapping = correct_response_table;
learning_context_mapping = learning_context_mapping(:, {'id', 'stim', 'feedback', 'modality'}); % reduce columns
[~, idx] = unique(learning_context_mapping(:, {'id','stim','feedback', 'modality'}), 'rows');
learning_context_mapping = learning_context_mapping(idx, :); % only keep unique rows
clear idx;
learning_context_mapping.id = categorical(learning_context_mapping.id); % ensure id variable can be recognized in loop

mapping = {}; % prepare empty object

for i = 1:numel(ids)
    
    % subset rows for current participant i
    id = ids(i);
    
    for j = 1:6
        
        Ti = learning_context_mapping(learning_context_mapping.id == id & learning_context_mapping.stim == j, :);
        
        r = 1; % check by default only first row
        
        while ismissing(Ti.feedback(r)) % but check if missing and in this case take the next row
            r=r+1
        end
        
        if Ti.feedback(r) == 1 & Ti.modality(r) == 1 || Ti.feedback(r) == 0 & Ti.modality(r) == 0
            
            context = 'getreward';
        
        elseif Ti.feedback(r) == 1 & Ti.modality(r) == 0 || Ti.feedback(r) == 0 & Ti.modality(r) == 1 
            
            context = 'avoidloss';
        end
        
        newRow = {id, j, context};
        mapping(end+1, :) = newRow;
        clear Ti r context;
    end
            
    clear id;
end

% convert cell to table and add column names
stim_context_mapping = cell2table(mapping, ...
        'VariableNames', {'id','stim','context'});

stim_context_mapping.id = categorical(stim_context_mapping.id);

% 2. Correct response mapping

correct_response_table = correct_response_table(correct_response_table.correct == 1, :); % reduce to rows with correct responses
correct_response_table = correct_response_table(:, {'id', 'stim', 'choice'}); % reduce columns
correct_response_table.Properties.VariableNames{'choice'} = 'correct_response'; %rename column

[~, idx] = unique(correct_response_table(:, {'id','stim'}), 'rows');
correct_response_table = correct_response_table(idx, :); % only keep unique rows


% To ensure that no row is missing (e.g. a specific stimulus for one
% participant):

pairs = [1 2; 3 4; 5 6]; % define stimulus pairs that have opposed correct responses
final_rows = {}; % prepare empty object

% ensure variable types are as intended
correct_response_table.id = categorical(correct_response_table.id); % categorical
correct_response_table.stim = double(correct_response_table.stim); % numeric
correct_response_table.correct_response = double(correct_response_table.correct_response); % numeric

for i = 1:numel(ids) % loop through participants
    
    % subset rows for current participant i
    id = ids(i);
    Ti = correct_response_table(correct_response_table.id == id, :);

    for p = 1:size(pairs,1) % loop through stimulus pairs
        s1 = pairs(p,1);
        s2 = pairs(p,2);

        has_s1 = any(Ti.stim == s1);
        has_s2 = any(Ti.stim == s2);

        if has_s1 && has_s2 % if both are present, next iteration of p
            continue
        end

        if ~has_s1 && has_s2 % if first stim is missing
            corr = Ti.correct_response(Ti.stim == s2); % get correct response for second stim
            alt = 3 - corr; % calculate alternative response (which can only be 1 or 2)
            newRow = {id, s1, alt}; % create new row
            final_rows(end+1, :) = newRow; %add row
        end

        if ~has_s2 && has_s1 % if second stim is missing
            corr = Ti.correct_response(Ti.stim == s1);
            alt = 3 - corr;
            newRow = {id, s2, alt};
            final_rows(end+1, :) = newRow;
        end
    end
end


if ~isempty(final_rows) % check if table exist (if there were missing rows)
    
    T_missing = cell2table(final_rows, ... % if yes: convert to table
        'VariableNames', {'id','stim','correct_response'});
    % and append to original table
    correct_response_table = [correct_response_table; T_missing];
end

% to check if  that worked (i.e. if correct response is included for each 
% of the six stimuli for each participant), uncomment & execute the
% following:
%sum(ismissing(T_final.correct_response)) % any missing rows in the correct response column?
%sum(crosstab(T_final.id,T_final.stim)) % each stimulus exactly as often as the table as there are participants?


% 3. (and last) step: rename ids in learning context and correct response 
% mapping table to anonymous ids according to the order in
% the cell filenames created in dataprep

load('interim_datasets\filenames.mat'); % load list
fname = cellfun(@(x)x{1}, filenames, 'UniformOutput', false);
fname = cellstr(fname);  % ensure filenames is a cellstr

for i = 1:numel(fname) % loop through filenames
    
    old_id = fname{i};  % extract filename
    new_id = sprintf('subj_%02d', i);  % create new anonymous filename
    
    % use both variables to replace filename in
    % (a) correct response table
    correct_response_table.id(correct_response_table.id == old_id) = categorical({new_id});
    % (b) stim to context mapping
    stim_context_mapping.id(stim_context_mapping.id == old_id) = categorical({new_id});
end

% sort by id and stim
correct_response_table = sortrows(correct_response_table, {'id','stim'});
% save for later reuse
correct_response_table.id = removecats(correct_response_table.id); % drop old codes
save('interim_datasets\datastruc_correct_response_key.mat', 'correct_response_table');

% sort by id and stim
stim_context_mapping = sortrows(stim_context_mapping, {'id','stim'});
% save for later reuse
stim_context_mapping.id = removecats(stim_context_mapping.id); % drop old codes
save('interim_datasets\datastruc_stim_context_mapping.mat', 'stim_context_mapping');

% add "new id" also to filenames table
new_ids = arrayfun(@(i) sprintf('subj_%02d', i), 1:numel(fname), 'UniformOutput', false);
filenames = table(fname(:), new_ids(:), 'VariableNames', {'original_id','new_id'});
% save for later reuse
save('interim_datasets\filenames_old_new.mat', 'filenames');

% save also as csv to anonnymize behav data in R
writetable(filenames, 'interim_datasets\filenames_old_new.csv');

% anonymized behavioural data for the R analyses
behav = readtable('interim_datasets/FBOmiss_behaviour_immediate.csv', 'TextType', 'string');
% run-2 logfiles carry the participant code plus "2": keep the first 6 characters
behav.id = extractBefore(behav.id, 7);
% replace each code by its anonymous id; stop if any code is not in the key
[found, idx] = ismember(behav.id, string(fname));
assert(all(found), 'Participant code without entry in the key');
ids = string(new_ids);
behav.id = reshape(ids(idx), [], 1);
% check: exactly 48 ids, 480 trials each
% check: exactly 48 ids, 480 trials each
[~, ~, g] = unique(behav.id);      % group index per row
counts = accumarray(g, 1);         % number of rows per id
assert(numel(counts) == 48 && all(counts == 480), ...
       'Expected 48 participants with 480 trials each');
   
writetable(behav, 'interim_datasets/FBOmiss_behaviour_immediate_anonymized.csv');

% overwrite non-anonymized filenames and save in wd
filenames = filenames.new_id;
save('interim_datasets\filenames.mat', 'filenames');
