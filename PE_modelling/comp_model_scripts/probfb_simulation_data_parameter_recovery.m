% Simulate data for parameter recovery;

% last modified, CW, 10/2026

%clear;

% load data and set hardcoded specifications

    % load fit and parameter information of best model
    addpath(genpath('./comp_model_fit_export/'));
    load('fit_model_5b_bic');

    % load overview of
    % (a) correct response mapping for stimuli and participants
    % (b) stimulus to learning context mapping
    addpath(genpath('./interim_datasets/'));
    load('datastruc_correct_response_key.mat');
    load('datastruc_stim_context_mapping.mat');

    nsubs = length(fit_BIC.bic); % get number of subjects
    nsim = 25; % specify number of simulated data sets (per participant)
    ntrials = 480; % specify maximum number of trials per participant
    
% simulate data and save in wd

for sub_i =1:nsubs % simulate data for each participant
          
    % 1. Extract individual parameter estimates of winning model of current
    % participant and concatenate:
    
    alpha_pospresented = fit_BIC.alpha_pospresented(sub_i);
    alpha_posomitted = fit_BIC.alpha_posomitted(sub_i);
    alpha_negpresented = fit_BIC.alpha_negpresented(sub_i);
    alpha_negomitted = fit_BIC.alpha_negomitted(sub_i);

    beta = fit_BIC.beta(sub_i);
    
    params_est = [alpha_pospresented, alpha_posomitted, ...
       alpha_negpresented, alpha_negomitted, beta];
    
    % 2. Simulate data (stimuli, feedback, choices, action values based on
    % stimulus-context-outcome mapping and estimated learning parameter):
   
    sim_n = sprintf('sim_data_%02d', sub_i) % print current participant
    
    for s=1:nsim
        
        sim_i = sprintf('sim_%02d', s); % current iteration
        
        % call simulation function (see below)
        simulated_data.(sim_n).(sim_i) = recovery_simulation(sub_i,...
            params_est, correct_response_table, stim_context_mapping);
     
    end
    
end

save('interim_datasets\simulated_data','simulated_data'); % save in wd



function sim_data = recovery_simulation(sub_i, params, correct_response_table, stim_context_mapping)


    % 1. Simulate stimulus presentation order (datastruc_stimuli)

    for b = 1:8 % total number of blocks

        block{b} = repelem(b,60); % number of block
        patterns = repmat([1:6],1,10); % create number of trials per stimulus id
        stim{b} = patterns(randperm(length(patterns))); % shuffle

    end

    % concatenate blocks and add trial number
    sim_data.trial = 1:480;
    sim_data.block = cat(2,block{:});
    sim_data.stim = cat(2,stim{:});

    % 2. Simulate pattern of chosen (datastruc_choice), and unchosen actions
    % (datastruc_not_choice), and probabilistic feedback (datastruc_feedback)
    % based on modelled learning parameter

    alpha_pospresented = params(1);
    alpha_posomitted = params(2);
    alpha_negpresented = params(3);
    alpha_negomitted = params(4);

    beta = params(5);

    % initialise action values for each action (right or left) in response to
    % each of the six stimuli (i.e. a total of 12)
    Q = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    p = NaN(1,length(sim_data.trial)); % NaN(M,N) is an M-by-N matrix of NaNs.

    for t = sim_data.trial % iterate through trials (here: 480)

       current_stim = sim_data.stim(t); % stim of current trial

       % first: check which action (left or right; coded as 1 or 2) was correct
       % in the experimental version (i.e. mapping of action to outcomes for 
       % the respective participant and stimulus):
       correct = correct_response_table.correct_response( ...
            correct_response_table.id ==  categorical(cellstr(sprintf('subj_%02d', sub_i))) ... % current participant
            & correct_response_table.stim == current_stim); % current stim

       % save in temporal variable q_i the current action values for the two
       % possible actions in response to the stimulus of the current trial and
       % store them such that the correct action is the first element in q_i:
       if correct == 1
           correct_index = current_stim+current_stim-1; % index for q-value of correct action
           incorrect_index = current_stim+current_stim; % index for q-value of incorrect action
           q_i = [Q(correct_index) Q(incorrect_index)]; % select qs
        elseif correct == 2
           correct_index = current_stim+current_stim; % index for q-value of correct action
           incorrect_index = current_stim+current_stim-1; % index for q-value of incorrect action
           q_i = [Q(correct_index) Q(incorrect_index)];
       end

       % calculate probability of making the correct decision based on fitted
       % exploration parameter and current action values:
       p = exp(beta*q_i(1))/(exp(beta*q_i(1)) + exp(beta*q_i(2))); 

       % generate random number between 0 and 1 and test if probability p (i.e.
       % probability to select correct action) is larger:  
       if rand(1) < p 
            response = correct; % if yes, select correct action
        else
            response = 3-correct; % if no, select incorrect
       end
       % logic: for high probabilities, random numbers will be more often
       % smaller and thus the correct action is more often selected and so on
       % (cf. Raab et al., 2020; https://doi.org/10.1038/s41598-020-72628-w)

       % save chosen, unchosen, and correct action in sim_data object
       sim_data.chosen(t) = response;
       sim_data.unchosen(t) = 3-response; % alternative response
       sim_data.best_choice(t) = correct;

       % add whether selected action is correct
       sim_data.accuracy(t) = NaN;
       if response == correct
          sim_data.accuracy(t) = 1;
        else
          sim_data.accuracy(t) = 0;
       end

       % simulate probabilistic feedback based on accuracy and outcome
       % probabilities (stim 1 & 2: 90% for correct action; stim 3 & 4: 70%;
       % stim 5 & 6: 50%; and for each inverse probability for incorrect 
       % action):
       if current_stim == 1 || current_stim == 2

           if sim_data.accuracy(t) == 1
                sim_data.feedback(t) = double(rand<0.9);
            else 
                sim_data.feedback(t) = double(rand<0.1);
           end

       elseif current_stim == 3 || current_stim == 4

           if sim_data.accuracy(t) == 1
                sim_data.feedback(t) = double(rand<0.7);
            else 
                sim_data.feedback(t) = double(rand<0.3);
           end

       elseif current_stim == 5 || current_stim == 6

           if sim_data.accuracy(t) == 1
                 sim_data.feedback(t) = double(rand<0.5);
            else % not sure if necessary
                 sim_data.feedback(t) = double(rand<0.5);
           end

       end


    % 3. Update action values

    % create an index for chosen and unchosen action
    chosen_index = NaN; % initialize
    unchosen_index = NaN;

    if sim_data.accuracy(t) == 1
        chosen_index = correct_index;
        unchosen_index = incorrect_index;
    elseif sim_data.accuracy(t) == 0
        chosen_index = incorrect_index;
        unchosen_index = correct_index;
    end

    % calculate prediction error of current trial for chosen and unchosen
    % action:
    PE_c = sim_data.feedback(t) - Q(chosen_index); % chosen PE
    PE_u = 1 - sim_data.feedback(t) - Q(unchosen_index); % unchosen PE

    % update action values

    % check learning context to select correct learning rate
    context = stim_context_mapping.context( ...
            stim_context_mapping.id ==  sprintf('subj_%02d', sub_i) ... % current participant
            & stim_context_mapping.stim == current_stim); % current stim
    context = string(context); % convert to string for if clause

    % initialize feedback type variable to add that in the loop
    sim_data.feedback_type(t) = NaN;

    % feedback type follows from learning context and outcome; learning
    % rates differ only by valence x appearance (not by stimulus)
    if context == "getreward" & sim_data.feedback(t) == 1 % presented pos FB
        Q(chosen_index) = Q(chosen_index) + alpha_pospresented*PE_c;
        Q(unchosen_index) = Q(unchosen_index) + alpha_pospresented*PE_u;
        sim_data.feedback_type(t) = 1;
    elseif context == "getreward" & sim_data.feedback(t) == 0 % omitted neg FB
        Q(chosen_index) = Q(chosen_index) + alpha_negomitted*PE_c;
        Q(unchosen_index) = Q(unchosen_index) + alpha_negomitted*PE_u;
        sim_data.feedback_type(t) = 0;
    elseif context == "avoidloss" & sim_data.feedback(t) == 1 % omitted pos FB
        Q(chosen_index) = Q(chosen_index) + alpha_posomitted*PE_c;
        Q(unchosen_index) = Q(unchosen_index) + alpha_posomitted*PE_u;
        sim_data.feedback_type(t) = 0;
    elseif context == "avoidloss" & sim_data.feedback(t) == 0 % presented neg FB
        Q(chosen_index) = Q(chosen_index) + alpha_negpresented*PE_c;
        Q(unchosen_index) = Q(unchosen_index) + alpha_negpresented*PE_u;
        sim_data.feedback_type(t) = 1;
    end

    end

    % convert sim_data (array of structs) to table
    fns = fieldnames(sim_data);

    for f = 1:length(fns)
        sim_data.(fns{f}) = sim_data.(fns{f})';
    end

    sim_data = struct2table(sim_data);

end
