%% FB Omission - Batch file

% CW; last modified in 10/2026

% This file calls several functions to (1) read raw logfiles and prepare
% the data to be fed in (2) several reinforcement learning models and save
% their fits. The two best fitting model are then used to (3) simulate action 
% values for each trial and export these for further analyses.

clear all; % clear workspace
rng(17, 'twister'); % fix random number generator for reproducibility
tic % start stopwatch

% total duration with 50 iterations: ~20 h (step 2 ~5.8 h, step 3 & 4 ~
% 14.8h)

%% 1. Reading logfiles and data preparation

% Change path to folder with data preparation scripts
addpath(genpath('./read_logfile_scripts_and_data_preparation/'));

% Call batch that runs the functions to read all relevant logfiles
% read_logfile_om_aktiv_batch_paper %% included hard coded filenames
% changed to loop
read_data;

% Create data structure expected by rl-model functions and save in folder
% 'interim_datasets' (specified in dataprep script)
[datastruc_stimuli, datastruc_choice,datastruc_not_choice,datastruc_feedback,datastruc_feedback_type,filenames,stim_rew_prob]=dataprep_PE_paper;

% Calculate means of feedback for each stimuli and each participant to
% identify which stimulus had which reward probability
stim_rew_prob.id = categorical(cellstr(stim_rew_prob.id));
stim_rew_prob.stim = cell2mat(stim_rew_prob.stim);
stim_rew_prob.feedback = cell2mat(stim_rew_prob.feedback);
stim_rew_prob.modality = cell2mat(stim_rew_prob.modality);
stim_rew_prob.choice = cell2mat(stim_rew_prob.choice);
stim_rew_prob.correct = cell2mat(stim_rew_prob.correct);

% Custom function to calculate mean while omitting NaN values in grpstats
omitnan_mean = @(x) mean(x, 'omitnan');


% Calculate the mean of feedback grouped by id and correct answers to
% verify stim numbering and positive outcome probability
stim_rew_prob_proportion_rewards = grpstats(stim_rew_prob, {'stim', 'correct'}, omitnan_mean, 'DataVars', 'feedback');
% stim 1&2 = 90% reward probability; stim 3&4 = 70%; stim 5&6 = 50% -
% confirmed

% Create tables in which for each participant mappings of correct responses
% and learning contexts to each stimuli are saved for later use in data
% simulation for parameter recovery:
create_correct_response_key;

%% 2. Fit rl models with varying complexity to the behavioural data

% 3 steps (as described in the manuscript): 

% 1. test learning rates with increasing complexity (separated by valence,
% appearance, context, rew prob)
% 2. test dynamic learning rates for (the two) best model (according to min
% BIC) identified above
% 3. test unchosen updates (+ separate learning rates for chosen and
% unchosen) for all of the above

% (fitting order below corresponds to order in Table S1; 
% to identify model for which dynamic learning rates were additionally 
% tested, all base models were fitted at first)

% model naming logic:

% first character: number indicating learning rate complexity (1 to 10; 
% same as in Table S1)

% second character:
% none or x = no update of unchosen action (in Table S1: X or X.0.X),
% b = update of unchosen action (in Table S1: X.1 or X.1.X),
% c = update of unchosen action and separate learning rates for chosen and 
% unchosen action (in Table S1: X.2 or X.2.X),

% third character:
% none = constant learning rates (in Table S1: X or X.X),
% a = exponential decay (in Table S1: X.X.1),
% b = RW-PH-dynamic (in Table S1: X.X.2).


% Change path to folder with rl-model functions
addpath(genpath('./comp_model_scripts/'));

% Create empty table to collect fit indices
variable_names_types = [["model", "string"]; ...
			["no_free_parameter", "double"]; ...
			["features", "string"]; ...
			["ll", "double"]; ...
			["BIC", "double"]];
		
no_of_models = 42; % only for setup of empty modelfits table
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
% feedback type/context/reward probability

probfb_funfit_model2 %calls function
modelfits(r,:) = {'2',3,'valence',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2_ll','fit');
save('comp_model_fit_export\fit_model_2_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 2xa: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability

% PLUS exponential decay of learning rates over presentations of the
% current stimulus (half-life lambda)

probfb_funfit_model2xa %calls function
modelfits(r,:) = {'2xa',4,'valence (+ exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2xa_ll','fit');
save('comp_model_fit_export\fit_model_2xa_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params


% Model 2xb: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability

% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

probfb_funfit_model2xb %calls function
modelfits(r,:) = {'2xb',5,'valence (+ RW-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2xb_ll','fit');
save('comp_model_fit_export\fit_model_2xb_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 2b: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability
% + update of chosen and unchosen option

probfb_funfit_model2b %calls function
modelfits(r,:) = {'2b',3,'valence (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2b_ll','fit');
save('comp_model_fit_export\fit_model_2b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 2ba: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability
% + update of chosen and unchosen option
% PLUS exponential decay of learning rates over presentations of the
% current stimulus (half-life lambda)

probfb_funfit_model2ba
modelfits(r,:) = {'2ba',4,'valence (+ unchosen update + exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2ba_ll','fit');
save('comp_model_fit_export\fit_model_2ba_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 2bb: separate learning rates for confirmatory and disconfirmatory
% trials (i.e. positive and negative feedback valence) across stimuli/
% feedback type/context/reward probability
% + update of chosen and unchosen option
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

probfb_funfit_model2bb
modelfits(r,:) = {'2bb',5,'valence (+ unchosen update + RW-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2bb_ll','fit');
save('comp_model_fit_export\fit_model_2bb_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 2c: separate learning rates for valence + choice
% (i.e. positive chosen, positive unchosen, and so on)
% + update of chosen and unchosen option

probfb_funfit_model2c %calls function
modelfits(r,:) = {'2c',5,'valence + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2c_ll','fit');
save('comp_model_fit_export\fit_model_2c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 2ca: separate learning rates for valence + choice
% (i.e. positive chosen, positive unchosen, and so on)
% + update of chosen and unchosen option
% PLUS exponential decay of learning rates over presentations of the
% current stimulus (half-life lambda)

probfb_funfit_model2ca %calls function
modelfits(r,:) = {'2ca',6,'valence + choice (+ unchosen update + exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2ca_ll','fit');
save('comp_model_fit_export\fit_model_2ca_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 2cb: separate learning rates for valence + choice
% (i.e. positive chosen, positive unchosen, and so on)
% + update of chosen and unchosen option
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

probfb_funfit_model2cb %calls function
modelfits(r,:) = {'2cb',7,'valence + choice (+ unchosen update + RW-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_2cb_ll','fit');
save('comp_model_fit_export\fit_model_2cb_bic','fit_BIC');
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

% Model 5xa: separate learning rates for valence + appearance // each
% possible feedback (i.e. positive presented, positive omitted, negative
% presented, negative omitted (inherently also separately for contexts))

% PLUS exponential decay of learning rates over presentations of the
% current stimulus (half-life lambda)

probfb_funfit_model5xa %calls function
modelfits(r,:) = {'5xa',6,'valence + appearance (+ exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5xa_ll','fit');
save('comp_model_fit_export\fit_model_5xa_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params


% Model 5xb: separate learning rates for valence + appearance // each
% possible feedback (i.e. positive presented, positive omitted, negative
% presented, negative omitted (inherently also separately for contexts))

% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

probfb_funfit_model5xb %calls function
modelfits(r,:) = {'5xb',7,'valence + appearance (+ RW-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5xb_ll','fit');
save('comp_model_fit_export\fit_model_5xb_bic','fit_BIC');
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

% Model 5ba: separate learning rates for valence + appearance // each
% possible feedback (i.e. positive presented, positive omitted, negative
% presented, negative omitted (inherently also separately for contexts))
% + update of chosen and unchosen option
% PLUS exponential decay of learning rates over presentations of the
% current stimulus (half-life lambda)

probfb_funfit_model5ba
modelfits(r,:) = {'5ba',6,'valence + appearance (+ unchosen update + exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5ba_ll','fit');
save('comp_model_fit_export\fit_model_5ba_bic','fit_BIC'); 
r=r+1;
clear fit fit_BIC params

% Model 5bb: separate learning rates for valence + appearance // each
% possible feedback (i.e. positive presented, positive omitted, negative
% presented, negative omitted (inherently also separately for contexts))
% + update of chosen and unchosen option
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

probfb_funfit_model5bb
modelfits(r,:) = {'5bb',7,'valence + appearance (+ unchosen update + RW-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5bb_ll','fit');
save('comp_model_fit_export\fit_model_5bb_bic','fit_BIC'); 
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

% Model 5ca: separate learning rates for valence + appearance + choice
% (i.e. positive presented chosen, positive presented unchosen, and so on)
% + update of chosen and unchosen option
% PLUS exponential decay of learning rates over presentations of the
% current stimulus (half-life lambda)

probfb_funfit_model5ca %calls function
modelfits(r,:) = {'5ca',10,'valence + appearance + choice (+ unchosen update + exponential decay)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5ca_ll','fit');
save('comp_model_fit_export\fit_model_5ca_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 5cb: separate learning rates for valence + appearance + choice
% (i.e. positive presented chosen, positive presented unchosen, and so on)
% + update of chosen and unchosen option
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

probfb_funfit_model5cb %calls function
modelfits(r,:) = {'5cb',11,'valence + appearance + choice (+ unchosen update + RW-PH dynamic)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_5cb_ll','fit');
save('comp_model_fit_export\fit_model_5cb_bic','fit_BIC');
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

% Model 10b:  learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.)
% PLUS update of unchosen option
  
probfb_funfit_model10b %calls function
modelfits(r,:) = {'10b',13,'valence + reward probability + appearance (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10b_ll','fit');
save('comp_model_fit_export\fit_model_10b_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params

% Model 10c:  learning rates separately for reward probabilities, and
% each possible feedback type (e.g. positive presented, positive omitted 
% etc.) and separately for chosen/unchosen
% PLUS update of unchosen option
  
probfb_funfit_model10c %calls function
modelfits(r,:) = {'10c',25,'valence + reward probability + appearance + choice (+ unchosen update)',mean(fit.ll), mean(fit_BIC.bic)};
save('comp_model_fit_export\fit_model_10c_ll','fit');
save('comp_model_fit_export\fit_model_10c_bic','fit_BIC');
r=r+1;
clear fit fit_BIC params


% add AIC to model fits
modelfits.AIC = 2*modelfits.ll + 2*modelfits.no_free_parameter;

% save table with modelfits of each model       
writetable(modelfits,'comp_model_fit_export\modelfits.csv');


% model with the lowest BIC (fit penalized for number of parameters/trials)
modelfits.BIC
[~, idx_bic] = min(modelfits.BIC);
disp('Model with lowest BIC:');
disp(modelfits(idx_bic,:)); % best model: 2b, very near 5b
 


 %% 3. Simulate PE for each trial with fitted parameters of the best models
% (here: model 2b (lowest BIC), model 5b (second lowest BIC) and save
% simulated values and parameters)

% to only run from hereon, uncomment the following
% tic
% Change path to folder with rl-model functions
% addpath(genpath('./comp_model_scripts/'));
% Specifiy number of fit iterations of fmincon for each function
% niter = 50;

probfb_simulation_best_model_BIC
probfb_simulation_best_model_5b

%% 4. Validate models by recovering parameter based on simulated data

% ---- Best BIC model (2b) ----

% 1. Simulate data
rng(17, 'twister');
probfb_simulation_data_parameter_recovery_BIC_model % for BIC model

% 2. Fit simulated data to the best model
probfb_funfit_model2b_simulated_data

% 3. Export recovered parameters as csv
% (one column per parameter and simulated data set, e.g. beta_01 ... beta_25)

load filenames;
export_names = param_names;
parameter_export = filenames;
varNames = {'filename'};
for p = 1:numel(export_names)
    vals = transpose(fit.(export_names{p})); % participants x simulations
    parameter_export = [parameter_export, num2cell(vals)];
    base = export_names{p};
    varNames = [varNames, arrayfun(@(x) sprintf('%s_%02d', base, x), 1:size(vals,2), 'UniformOutput', false)];
end
parameter_export = cell2table(parameter_export, 'VariableNames', varNames);

% export as csv
writetable(parameter_export,'comp_model_fit_export\FBOmiss_learning_parameter_recovered_BIC.csv');


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
writetable(behav_table,'comp_model_fit_export\FBOmiss_behav_recovered_BIC.csv');




% ---- Second Best BIC model (5b) ----

% duration ~ 2.2h

% 1. Simulate data
rng(17, 'twister');
probfb_simulation_data_parameter_recovery_5b_model % for model 5b (second-lowest BIC)

% 2. Fit simulated data to model 5b
probfb_funfit_model5b_simulated_data

% 3. Export recovered parameters as csv
% (one column per parameter and simulated data set, e.g. beta_01 ... beta_25)

load filenames;
export_names = param_names;
parameter_export = filenames;
varNames = {'filename'};
for p = 1:numel(export_names)
    vals = transpose(fit.(export_names{p})); % participants x simulations
    parameter_export = [parameter_export, num2cell(vals)];
    base = export_names{p};
    varNames = [varNames, arrayfun(@(x) sprintf('%s_%02d', base, x), 1:size(vals,2), 'UniformOutput', false)];
end
parameter_export = cell2table(parameter_export, 'VariableNames', varNames);

% export as csv
writetable(parameter_export,'comp_model_fit_export\FBOmiss_learning_parameter_recovered_5b.csv');


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
writetable(behav_table,'comp_model_fit_export\FBOmiss_behav_recovered_5b.csv');

toc % stopwatch end