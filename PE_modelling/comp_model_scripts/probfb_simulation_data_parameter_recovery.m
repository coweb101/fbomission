% 4. Simulate data for parameter recovery and posterior predictive check
% Best model (lowest BIC): model 2b: separate learning rates for positive and negative feedback (+ update of chosen and unchosen option)
% (generated from probfb_fun_model2b.m; model update copied verbatim)

% load fit and parameter information of best model
addpath(genpath('./comp_model_fit_export/'));
load('fit_model_2b_bic');

% names of the free parameters in the order of params(1), params(2), ...
param_names = {'alpha_pos', 'alpha_neg', 'beta'};

% load overview of
% (a) correct response mapping for stimuli and participants
% (b) stimulus to learning context mapping
addpath(genpath('./interim_datasets/'));
load('datastruc_correct_response_key.mat');
load('datastruc_stim_context_mapping.mat');

nsubs = length(fit_BIC.bic); % number of participants
nsim = 25; % number of simulated data sets per participant
ntrials = 480; % number of trials per data set

for sub_i = 1:nsubs % simulate data for each participant

    % parameter estimates of the current participant
    params_est = zeros(1, numel(param_names));
    for k = 1:numel(param_names)
        params_est(k) = fit_BIC.(param_names{k})(sub_i);
    end

    sim_n = sprintf('sim_data_%02d', sub_i) % print current participant

    for s = 1:nsim
        sim_i = sprintf('sim_%02d', s); % current iteration
        simulated_data.(sim_n).(sim_i) = recovery_simulation(sub_i, ...
            params_est, correct_response_table, stim_context_mapping, ntrials);
    end

end

save('interim_datasets\simulated_data','simulated_data');


function sim_data = recovery_simulation(sub_i, params, correct_response_table, stim_context_mapping, ntrials)
% Simulates one data set (stimulus order, choices, probabilistic feedback)
% with model 2b. Parameter assignment and update of action values are
% copied verbatim from probfb_fun_model2b.m.

    % 1. Stimulus order: 8 blocks of 60 trials, each stimulus 10 times per block
    for b = 1:8
        block{b} = repelem(b,60);
        patterns = repmat(1:6,1,10);
        stim{b} = patterns(randperm(length(patterns))); % shuffle
    end
    sim_data.trial = 1:ntrials;
    sim_data.block = cat(2,block{:});
    sim_data.stim = cat(2,stim{:});

    % 2. Parameters
    % Separate learning rates for positive and negative FB
    % + chosen and unchosen action updated

    % set parameters
    alpha_pos = params(1);
    alpha_neg = params(2);
    beta = params(3);

    % action values for each action (left/right) for each of the six stimuli
    Q = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];

    % trial-wise inputs in the format used by the model function
    sub_stimuli = sim_data.stim;
    sub_fb = NaN(1,ntrials);
    sub_feedback_type = NaN(1,ntrials);

    % participant-specific mappings
    sub_key = correct_response_table(correct_response_table.id == categorical(cellstr(sprintf('subj_%02d', sub_i))), :);
    sub_context = stim_context_mapping(stim_context_mapping.id == sprintf('subj_%02d', sub_i), :);

    for i = 1:ntrials

        current_stim = sub_stimuli(i);


        % correct action (1 = left, 2 = right) for this participant and stimulus
        correct = sub_key.correct_response(sub_key.stim == current_stim);
        if correct == 1
            correct_index = current_stim+current_stim-1;
            incorrect_index = current_stim+current_stim;
        else
            correct_index = current_stim+current_stim;
            incorrect_index = current_stim+current_stim-1;
        end

        % probability of the correct action (softmax) and simulated response
        p_correct = exp(beta*Q(correct_index))/(exp(beta*Q(correct_index)) + exp(beta*Q(incorrect_index)));
        if rand(1) < p_correct
            response = correct;
        else
            response = 3-correct;
        end
        sim_data.chosen(i) = response;
        sim_data.unchosen(i) = 3-response;
        sim_data.best_choice(i) = correct;
        sim_data.accuracy(i) = double(response == correct);

        % simulate probabilistic feedback based on accuracy and outcome
        % probabilities (stim 1 & 2: 90% for correct action; stim 3 & 4: 70%;
        % stim 5 & 6: 50%; and for each inverse probability for incorrect 
        % action):
        if current_stim == 1 || current_stim == 2

            if sim_data.accuracy(i) == 1
                 sim_data.feedback(i) = double(rand<0.9);
             else 
                 sim_data.feedback(i) = double(rand<0.1);
            end

        elseif current_stim == 3 || current_stim == 4

            if sim_data.accuracy(i) == 1
                 sim_data.feedback(i) = double(rand<0.7);
             else 
                 sim_data.feedback(i) = double(rand<0.3);
            end

        elseif current_stim == 5 || current_stim == 6

            if sim_data.accuracy(i) == 1
                  sim_data.feedback(i) = double(rand<0.5);
             else % not sure if necessary
                  sim_data.feedback(i) = double(rand<0.5);
            end

        end

        % feedback type follows from learning context and outcome:
        % get reward: positive = presented, negative = omitted
        % avoid loss: positive = omitted, negative = presented
        context = string(sub_context.context(sub_context.stim == current_stim));
        sub_fb(i) = sim_data.feedback(i);
        if context == "getreward"
            sub_feedback_type(i) = double(sub_fb(i) == 1);
        else
            sub_feedback_type(i) = double(sub_fb(i) == 0);
        end
        sim_data.feedback_type(i) = sub_feedback_type(i);

        % chosen and unchosen action (index of their action values)
        if response == 1
            chosen_index = current_stim+current_stim-1;
            unchosen_index = current_stim+current_stim;
        else
            chosen_index = current_stim+current_stim;
            unchosen_index = current_stim+current_stim-1;
        end

        % prediction errors of chosen and unchosen action
        PE_c = sub_fb(i) - Q(chosen_index);
        PE_u = 1 - sub_fb(i) - Q(unchosen_index);

        % ---- update (copied from probfb_fun_model2b.m) ----
        % separate learning rates for confirmatory and disconfirmatory
        % trials (i.e., positive and negative fb which is associated with a
        % positive and negative PE, repsectively)
        if sub_fb(i) == 1 % if fb is positive

            % update action value of chosen action with alpha_pos
            Q(chosen_index) = Q(chosen_index) + alpha_pos*PE_c;
            Q(unchosen_index) = Q(unchosen_index) + alpha_pos*PE_u;

        elseif sub_fb(i) == 0 % if fb is negative

            % use alpha_neg to update action value of chosen 
            Q(chosen_index) = Q(chosen_index)+ alpha_neg*PE_c;
            Q(unchosen_index) = Q(unchosen_index) + alpha_neg*PE_u;

        end
    end

    % convert sim_data (struct of row vectors) to table
    fns = fieldnames(sim_data);
    for f = 1:length(fns)
        sim_data.(fns{f}) = sim_data.(fns{f})';
    end
    sim_data = struct2table(sim_data);
end
