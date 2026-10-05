% Create a list of filenames in the given folder containing the string 
% 'omission_':
folder = 'logfiles';
files = getfilenames(folder, 'omission_');

% sort files by modification date
d = dir(fullfile(folder, '*'));
d = d(~[d.isdir]);                         % remove folders

% keep only the files that are in your filename list
d = d(ismember({d.name}', files));

[~, idx] = sort([d.datenum]);              % oldest to newest
files = {d(idx).name}';                    % sorted filename list

subjCodes = cellfun(@(x) regexp(x,'^[^-]+','match','once'), files, 'UniformOutput', false);
baseIDs   = cellfun(@(x) x(1:6), subjCodes, 'UniformOutput', false);
uniqueIDs = unique(baseIDs, 'stable');

for s = 1:numel(uniqueIDs)
    id  = uniqueIDs{s};      % e.g. AB12DE
    id2 = [id '2'];          % e.g. AB12DE2
    vpName = sprintf('VP%02d', s);

    % find filenames belonging to run 1 / run 2
    idx1 = find(strcmp(subjCodes, id));
    idx2 = find(strcmp(subjCodes, id2));

    if numel(idx1) ~= 1 || numel(idx2) ~= 1
        error('Expected exactly one file for %s and one for %s', id, id2);
    end

    file1 = files{idx1};
    file2 = files{idx2};

    % extract suffix after "-omission_" and before ".log"
    suf1 = regexp(file1, '(?<=-omission_)[ab][12](?=\.log$)', 'match', 'once');
    suf2 = regexp(file2, '(?<=-omission_)[ab][12](?=\.log$)', 'match', 'once');

    % first run of current subject
    if endsWith(id, '2')
        % if your run-1 code could end with 2, use the other function
        eval(sprintf([ ...
            '[%s_all_data_1,%s_learn_data_1,%s_resp_ti_data_1,%s_stat_outcome_1] = ' ...
            'read_logfile_om_aktiv_2(''%s'', ''%s-omission_%s'');'], ...
            vpName, vpName, vpName, vpName, id, id, suf1));
    else
        eval(sprintf([ ...
            '[%s_all_data_1,%s_learn_data_1,%s_resp_ti_data_1,%s_stat_outcome_1] = ' ...
            'read_logfile_om_aktiv_1(''%s'', ''%s-omission_%s'');'], ...
            vpName, vpName, vpName, vpName, id, id, suf1));
    end

    % second run of current subject
    if endsWith(id2, '2')
        eval(sprintf([ ...
            '[%s_all_data_2,%s_learn_data_2,%s_resp_ti_data_2,%s_stat_outcome_2] = ' ...
            'read_logfile_om_aktiv_2(''%s'', ''%s-omission_%s'');'], ...
            vpName, vpName, vpName, vpName, id2, id2, suf2));
    else
        eval(sprintf([ ...
            '[%s_all_data_2,%s_learn_data_2,%s_resp_ti_data_2,%s_stat_outcome_2] = ' ...
            'read_logfile_om_aktiv_1(''%s'', ''%s-omission_%s'');'], ...
            vpName, vpName, vpName, vpName, id2, id2, suf2));
    end

    % combine 
    eval(sprintf('%s_all_data=[%s_all_data_1;%s_all_data_2];', vpName, vpName, vpName));
    eval(sprintf('%s_learn_data=[%s_learn_data_1;%s_learn_data_2];', vpName, vpName, vpName));
    eval(sprintf('learn_data_all(:,:,%d)=%s_learn_data;', s, vpName));
    eval(sprintf('%s_resp_ti_data=[%s_resp_ti_data_1;%s_resp_ti_data_2];', vpName, vpName, vpName));
    eval(sprintf('learn_resp_ti_all(:,:,%d)=%s_resp_ti_data;', s, vpName));

    % clear
    eval(sprintf('clear %s_all_data_1', vpName));
    eval(sprintf('clear %s_all_data_2', vpName));
    eval(sprintf('clear %s_learn_data_1', vpName));
    eval(sprintf('clear %s_learn_data_2', vpName));
    eval(sprintf('clear %s_resp_ti_data_1', vpName));
    eval(sprintf('clear %s_resp_ti_data_2', vpName));
    
    
end

clearvars -except VP* learn_data_all learn_resp_ti_all
save('interim_datasets/all_vp_behav_paper.mat')

function [filenames]=getfilenames(folder, varargin)
% Read filenames from a folder and output a cell with all names matching
% the conditions set in varargin.
% by Alexander Seidel - 2017
%
% Note: I adjusted the function name to include as a local function within
% my script (not possible if function call and created variable (filenames)
% is identical. CW/2023
%
% INPUT
% folder [string]           Total or relative path to the folder.
% varargin [string]         Arbitrary amount of strings to specifiy
%                           conditions. Condition strings starting with '-'
%                           are used to exclude entries from the file name
%                           list, all others are form a requirement each
%                           entry must meet. Each entry must meet all
%                           requirements set by the conditions.
%
% OUTPUT
% files [cell]              vertical cell array with all filenames that
%                           meet the specified conditions.

    % Read filenames and remove all entries that are folders
    filenames = dir(folder);
    filenames = {filenames([filenames.isdir] == 0).name}';
    
    if isempty(filenames)
        error('The Folder is empty.');
    end
    
    % Remove entries not matching the conditions
    if ~isempty(varargin)
        for i=1:length(varargin)
            
            if varargin{i}(1) == '/'
                rows = ~strfindl(filenames,varargin{i}(2:end));
            else
                rows = strfindl(filenames,varargin{i});
            end
            
            filenames = filenames(rows);
        end
    end
    
    if isempty(filenames)
        error('No files matching the criteria were found.');
    end

end

function index=strfindl(str, pattern)
% A logical version of strfind

    strfound = strfind(str,pattern);
    index = cell2mat(cellfun(@(x) ~isempty(x),strfound,'uni', false));

end