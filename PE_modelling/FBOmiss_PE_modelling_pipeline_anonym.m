%% FB Omission - PE modelling pipeline

% CW; last modified in 01/2026

clear all; % clear workspace
tic % start stopwatch | duration for 50 iterations: ~313 minutes

%% 1. Read data

% Omitted for anonymization
% everything below is reproducable with anonymized matlab data files

%% 2. Fit rl models with varying complexity to the behavioural data (extract
% best fit and compare indices)

% Change path to folder with rl-model functions
addpath(genpath('./comp_model_scripts/'));

% Create empty table to collect fit indices
variable_names_types = [["model", "string"]; ...
			["no_free_parameter", "double"]; ...
			["features", "string"]; ...
			["ll", "double"]; ...
			["BIC", "double"]];
		
no_of_models = 36; % only for setup of empty modelfits table
r=1; % set counter to 1 for inserting model info to the table

% Create table using fieldnames & value types from above
modelfits = table('Size',[no_of_models,size(variable_names_types,1)],... 
	'VariableNames', variable_names_types(:,1),...
	'VariableTypes', variable_names_types(:,2));

% Specifiy number of fit iterations of fmincon for each function
niter = 50;

% Model 1: one global learning rate across stimuli/valence/feedback
% type/context/reward probability + update of chosen option only
    
probfb_funfit_model1 % calls function
%saves model info and fit indices
modelfits(r,:) = {'1',2,'global',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_1_ll','fit');
save('comp_model_fit_export\fit_model_1_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 1b: one global learning rate + update of chosen and unchosen option
    
probfb_funfit_model1b
modelfits(r,:) = {'1b',2,'global (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_1b_ll','fit');
save('comp_model_fit_export\fit_model_1b_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 1c: update of chosen and unchosen option + separate learning rates
% for chosen and unchosen action
    
probfb_funfit_model1c
modelfits(r,:) = {'1c',3,'choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_1c_ll','fit');
save('comp_model_fit_export\fit_model_1c_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 2: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability + update of chosen option only

probfb_funfit_model2
modelfits(r,:) = {'2',3,'valence',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2_ll','fit');
save('comp_model_fit_export\fit_model_2_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 2b: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability + update of chosen and unchosen
% option

probfb_funfit_model2b
modelfits(r,:) = {'2b',3,'valence (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2b_ll','fit');
save('comp_model_fit_export\fit_model_2b_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 2c: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) and chosen and 
% unchosen
% + update of chosen and unchosen option

probfb_funfit_model2c
modelfits(r,:) = {'2c',5,'valence + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2c_ll','fit');
save('comp_model_fit_export\fit_model_2c_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 3: separate learning rates for presented and omitted feedback
% trials across stimuli/valence/context/reward probability

probfb_funfit_model3 %calls function
modelfits(r,:) = {'3',3,'appearance',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_3_ll','fit');
save('comp_model_fit_export\fit_model_3_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 3b: separate learning rates for presented and omitted feedback
% trials across stimuli/valence/context/reward probability
% + update of chosen and unchosen option

probfb_funfit_model3b %calls function
modelfits(r,:) = {'3b',3,'appearance (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_3b_ll','fit');
save('comp_model_fit_export\fit_model_3b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 3c: separate learning rates for presented and omitted feedback
% trials, and chosen and unchosen action (but across stimuli/valence/
% context/reward probability)
% + update of chosen and unchosen option

probfb_funfit_model3c %calls function
modelfits(r,:) = {'3c',5,'appearance + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_3c_ll','fit');
save('comp_model_fit_export\fit_model_3c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params
    
% Model 4: separate learning rates for get reward and avoid loss 
% context  
    
probfb_funfit_model4 %calls function
modelfits(r,:) = {'4',3,'context',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_4_ll','fit');
save('comp_model_fit_export\fit_model_4_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 4b: separate learning rates for get reward and avoid loss 
% context + update of chosen and unchosen option
    
probfb_funfit_model4b %calls function
modelfits(r,:) = {'4b',3,'context (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_4b_ll','fit');
save('comp_model_fit_export\fit_model_4b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 4c: separate learning rates for get reward and avoid loss 
% context, and chosen and unchosen action    
% + update of chosen and unchosen option
    
probfb_funfit_model4c %calls function
modelfits(r,:) = {'4c',5,'context + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_4c_ll','fit');
save('comp_model_fit_export\fit_model_4c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 5: separate learning rates for valence + appearance // each
% possible feedback (i.e. positive presented, positive omitted, negative
% presented, negative omitted (inherently also separately for contexts))

probfb_funfit_model5 %calls function
modelfits(r,:) = {'5',5,'valence + appearance',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5_ll','fit');
save('comp_model_fit_export\fit_model_5_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 5b: separate learning rates for valence + appearance // each
% possible feedback (i.e. positive presented, positive omitted, negative
% presented, negative omitted (inherently also separately for contexts))
% + update of chosen and unchosen option

probfb_funfit_model5b %calls function
modelfits(r,:) = {'5b',5,'valence + appearance (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5b_ll','fit');
save('comp_model_fit_export\fit_model_5b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 5c: separate learning rates for valence + appearance + choice
% (i.e. positive presented chosen, positive presented unchosen, and so on)
% + update of chosen and unchosen option

probfb_funfit_model5c %calls function
modelfits(r,:) = {'5c',9,'valence + appearance + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5c_ll','fit');
save('comp_model_fit_export\fit_model_5c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 6: one learning rate per stimulus reward probabilities
% (i.e. a total of 3 learning rates)

probfb_funfit_model6 %calls function
modelfits(r,:) = {'6',4,'reward probability',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_6_ll','fit');
save('comp_model_fit_export\fit_model_6_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 6b: one learning rate per stimulus reward probabilities
% (i.e. a total of 3 learning rates)
% + update of chosen and unchosen option

probfb_funfit_model6b %calls function
modelfits(r,:) = {'6b',4,'reward probability (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_6b_ll','fit');
save('comp_model_fit_export\fit_model_6b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 6c: one learning rate per stimulus reward probabilities
% and separately for chosen and unchosen action
% + update of unchosen option

probfb_funfit_model6c %calls function
modelfits(r,:) = {'6c',7,'reward probability + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_6c_ll','fit');
save('comp_model_fit_export\fit_model_6c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 7: learning rates separately for stimulus reward probabilities
% and separately for valence but across feedback type/context (i.e. a 
% total of 6 learning rates)
    
probfb_funfit_model7 %calls function
modelfits(r,:) = {'7',7,'valence + reward probability',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_7_ll','fit');
save('comp_model_fit_export\fit_model_7_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 7b: learning rates separately for stimulus reward probabilities
% and separately for valence but across feedback type/context (i.e. a 
% total of 6 learning rates)
% + update of chosen and unchosen option

probfb_funfit_model7b %calls function
modelfits(r,:) = {'7b',7,'valence + reward probability (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_7b_ll','fit');
save('comp_model_fit_export\fit_model_7b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 7c: learning rates separately for stimulus reward probabilities,
% valence and choice
% + update of chosen and unchosen option

probfb_funfit_model7c %calls function
modelfits(r,:) = {'7c',13,'valence + reward probability + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_7c_ll','fit');
save('comp_model_fit_export\fit_model_7c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 8: learning rates separately for stimulus reward probabilities
% and feedback appearance
    
probfb_funfit_model8 %calls function
modelfits(r,:) = {'8',7,'appearance + reward probability',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_8_ll','fit');
save('comp_model_fit_export\fit_model_8_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 8b: learning rates separately for stimulus reward probabilities
% and feedback appearance
% + update of chosen and unchosen option

probfb_funfit_model8b %calls function
modelfits(r,:) = {'8b',7,'appearance + reward probability (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_8b_ll','fit');
save('comp_model_fit_export\fit_model_8b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 8c: learning rates separately for stimulus reward probabilities,
% feedback appearance, and choice
% + update of chosen and unchosen option

probfb_funfit_model8c %calls function
modelfits(r,:) = {'8c',13,'appearance + reward probability + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_8c_ll','fit');
save('comp_model_fit_export\fit_model_8c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 9: learning rates separately for stimulus reward probabilities
% and separately for context (i.e. a total of 6 learning rates; this 
% corresponds to a learning rate for each individual stimulus of the six
% used stimuli as the two stimuli assigned to one stimulus reward
% probability were shown either in the reward OR in the loss context
    
probfb_funfit_model9 %calls function
modelfits(r,:) = {'9',7,'learning context + reward probability',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_9_ll','fit');
save('comp_model_fit_export\fit_model_9_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 9b: learning rates separately for stimulus reward probabilities
% and separately for context (i.e. a total of 6 learning rates; this 
% corresponds to a learning rate for each individual stimulus of the six
% used stimuli as the two stimuli assigned to one stimulus reward
% probability were shown either in the reward OR in the loss context
% + update of chosen and unchosen option
    
probfb_funfit_model9b %calls function
modelfits(r,:) = {'9b',7,'learning context + reward probability (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_9b_ll','fit');
save('comp_model_fit_export\fit_model_9b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 9c: learning rates separately for stimulus reward probabilities,
% context and choice (i.e. a total of 12 learning rates; this 
% corresponds to a learning rate for each individual action to each
% stimulus of the six used stimuli as the two stimuli assigned to one
% stimulus reward probability were shown either in the reward OR in the
% loss context)
% + update of chosen and unchosen option
    
probfb_funfit_model9c %calls function
modelfits(r,:) = {'9c',13,'learning context + reward probability + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_9c_ll','fit');
save('comp_model_fit_export\fit_model_9c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10: learning rates separately for stimulus reward probabilities
% and separately for each possible feedback type (i.e. positive presented,
% positive omitted, negative presented, negative omitted; this is 
% inherently also separately for contexts; i.e. a total of 12 learning
% rates) - this corresponds to individual learning rates for each
% individual stimulus separately for positive and negative feedback as the
% two stimuli assigned to one reward probability were shown either in the
% reward OR in the loss context
  
probfb_funfit_model10 %calls function
modelfits(r,:) = {'10',13,'valence + reward probability + appearance',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10_ll','fit');
save('comp_model_fit_export\fit_model_10_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10a:  learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.)
% PLUS update of unchosen option
  
probfb_funfit_model10a %calls function
modelfits(r,:) = {'10a',13,'valence + reward probability + appearance (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10a_ll','fit');
save('comp_model_fit_export\fit_model_10a_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10aac:  learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.) and separately for chosen/unchosen
% PLUS update of unchosen option
  
probfb_funfit_model10aac %calls function
modelfits(r,:) = {'10aac',25,'valence + reward probability + appearance + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10aac_ll','fit');
save('comp_model_fit_export\fit_model_10aac_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10b: learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.)
% PLUS learning rate update in each trial with an exponential decay
% parameter (lambda)
  
probfb_funfit_model10b %calls function
modelfits(r,:) = {'10b',14,'valence + reward probability + appearance (+ exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10b_ll','fit');
save('comp_model_fit_export\fit_model_10b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10ab: learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.)
% PLUS learning rate update in each trial with an exponential decay
% parameter (lambda)
% PLUS update of unchosen option

probfb_funfit_model10ab %calls function
modelfits(r,:) = {'10ab',14,'valence + reward probability + appearance (+ exponential decay + unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10ab_ll','fit');
save('comp_model_fit_export\fit_model_10ab_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10abc: learning rates separately for reward probabilities,
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.) and choice (unchosen vs chosen)
% PLUS learning rate update in each trial with an exponential decay
% parameter (lambda)
% PLUS update of unchosen option

probfb_funfit_model10abc %calls function
modelfits(r,:) = {'10abc',26,'valence + reward probability + appearance + choice (+ exponential decay + unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10abc_ll','fit');
save('comp_model_fit_export\fit_model_10abc_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10c: learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.)
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE
  
probfb_funfit_model10c %calls function
modelfits(r,:) = {'10c',15,'valence + reward probability + appearance (+ RL-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10c_ll','fit');
save('comp_model_fit_export\fit_model_10c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10ac: learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.)
% PLUS learning rate update in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE
% PLUS update of unchosen option

probfb_funfit_model10ac %calls function
modelfits(r,:) = {'10ac',15,'valence + reward probability + appearance (+ RL-PH dynamic + unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10ac_ll','fit');
save('comp_model_fit_export\fit_model_10ac_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10acc: learning rates separately for reward probabilities,
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.) and separately for chosen and unchosen
% PLUS learning rate update in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE
% PLUS update of unchosen option
  
probfb_funfit_model10acc %calls function
modelfits(r,:) = {'10acc',27,'valence + reward probability + appearance + choice (+ RL-PH dynamic + unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10acc_ll','fit');
save('comp_model_fit_export\fit_model_10acc_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% save table with modelfits of each model       
writetable(modelfits,'comp_model_fit_export\modelfits_immediate.csv');  
 

%% 3. Simulate PE for each trial with fitted parameters of best fitting
% model and save simulated values as well as parameters (caution this is 
% hardcoded/manually prepared once the best model was identified)

probfb_simulation_best_model

%% 4. Validate model by recovering parameter based on simulated data

% 1. Simulate data
probfb_simulation_data_parameter_recovery

% 2. Fit simulated data to previously identified model
probfb_funfit_model10a_simulated_data
save('comp_model_fit_export\fit_model_recovery_ll','fit');
save('comp_model_fit_export\fit_model_recovery_bic','fit_BIC');

% 3. Rebuild fit object and export fits and recovered parameter as csv

% concatenate filenames, BICs and recovered alphas & beta
load filenames;
parameter_export = cat(2, filenames,num2cell(transpose(fit_BIC.bic)), ...
    num2cell(transpose(fit_BIC.alpha_90_presented_pos)), ...
    num2cell(transpose(fit_BIC.alpha_70_presented_pos)), ...
    num2cell(transpose(fit_BIC.alpha_50_presented_pos)), ...
    num2cell(transpose(fit_BIC.alpha_90_omitted_pos)), ...
    num2cell(transpose(fit_BIC.alpha_70_omitted_pos)), ...
    num2cell(transpose(fit_BIC.alpha_50_omitted_pos)), ...
    num2cell(transpose(fit_BIC.alpha_90_presented_neg)), ...
    num2cell(transpose(fit_BIC.alpha_70_presented_neg)), ...
    num2cell(transpose(fit_BIC.alpha_50_presented_neg)), ...
    num2cell(transpose(fit_BIC.alpha_90_omitted_neg)), ...
    num2cell(transpose(fit_BIC.alpha_70_omitted_neg)), ...
    num2cell(transpose(fit_BIC.alpha_50_omitted_neg)), ...    
    num2cell(transpose(fit.beta)));
% columns: filenames, BIC, alphas..., beta

% assign column labels
baseNames = {'BIC', ...
    'alpha_90_presented_pos', 'alpha_70_presented_pos', 'alpha_50_presented_pos', ...
    'alpha_90_omitted_pos',   'alpha_70_omitted_pos',   'alpha_50_omitted_pos', ...
    'alpha_90_presented_neg', 'alpha_70_presented_neg', 'alpha_50_presented_neg', ...
    'alpha_90_omitted_neg',   'alpha_70_omitted_neg',   'alpha_50_omitted_neg', ...
    'beta'};

expandedNames = cellfun(@(b) arrayfun( ...
    @(x) sprintf('%s_%02d', b, x), 1:25, 'UniformOutput', false), ...
    baseNames, 'UniformOutput', false);

expandedNames = [expandedNames{:}];   % flatten into 1×N cell
varNames = ['filename', expandedNames];
parameter_export = cell2table(parameter_export, 'VariableNames', varNames);

% export as csv
writetable(parameter_export,'comp_model_fit_export\FBOmiss_learning_parameter_recovered.csv');


% 4. Export simulated_data to calculate accuracy (posterior predictive check)

% simulated_data is a 1x1 struct with 48 fields (nsubs) with each 25 fields
% with each having a 480x9 table
firstlevel = fieldnames(simulated_data); 
nfirst = numel(firstlevel); % get number of fields at first level
nsecond = 25; % number of second level fields (= number of sim data per participant)
ntrials = 480; % rows per table in second level field
nrows = nfirst * nsecond * ntrials; % number of total rows in final table

out = cell(nrows, 11); % prepare output object
row = 1; % initialize counter

% loop through simulated data struct
for i = 1:nfirst
    fieldname_firstlevel = firstlevel{i};
    substruct = simulated_data.(fieldname_firstlevel);
    secondlevel = fieldnames(substruct);

    for j = 1:nsecond
        fieldname_secondlevel = secondlevel{j}; % get current sim label
        table_tmp = substruct.(fieldname_secondlevel); % get table
        
        % prepare row index
        row_index = row:row+ntrials-1;
        
        % fill rows
        out(row_index,1) = {fieldname_firstlevel};
        out(row_index,2) = {fieldname_secondlevel};
        out(row_index,3:11) = table2cell(table_tmp);
        
        row = row + ntrials; % add to counter
    end
end

% convert to table
table_colnames = table_tmp.Properties.VariableNames;

columnnames_all = [{'subj','sim_iteration'}, table_colnames];
behav_table = cell2table(out, 'VariableNames', columnnames_all);

% export as csv
writetable(behav_table,'comp_model_fit_export\FBOmiss_behav_recovered.csv');


toc % stopwatch end