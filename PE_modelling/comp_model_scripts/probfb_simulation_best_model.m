% 3. Simulate action values and prediction errors for each trial with the
% fitted parameters of the best model (lowest BIC: model 2b)
% Model 2b: separate learning rates for positive and negative feedback (+ update of chosen and unchosen option)
% (generated from probfb_fun_model2b.m; model update copied verbatim)

% set directory for model fit
addpath(genpath('./comp_model_fit_export/'));

% load parameters of best fit of model 2b (fit_BIC: within a model the
% iteration with the lowest BIC is also the one with the lowest -LL and AIC)
load fit_model_2b_bic;

% names of the free parameters in the order of params(1), params(2), ...
param_names = {'alpha_pos', 'alpha_neg', 'beta'};

% set directory for behavioural data (prepared with dataprep)
addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep
load datastruc_feedback_type; % created with dataprep
load filenames; % anonymous ids (created with create_correct_response_key)

ntrials = 480;

for i = 1:length(fit_BIC.bic) % loop through participants

    i

    % retrieve parameters of best fit
    params = zeros(1, numel(param_names));
    for k = 1:numel(param_names)
        params(k) = fit_BIC.(param_names{k})(i);
    end

    % prepare choice and outcome for simulation function
    sub_choice = cell2mat(datastruc_choice(i,:));
    sub_not_choice = cell2mat(datastruc_not_choice(i,:));
    sub_outcome = cell2mat(datastruc_feedback(i,:));
    sub_stimuli = cell2mat(datastruc_stimuli(i,:));
    sub_feedback_type = cell2mat(datastruc_feedback_type(i,:));

    % simulate values based on parameters of best fit
    currsim = sim_values_best_model(params, sub_choice, sub_outcome, sub_stimuli, sub_feedback_type);

    % save simulated values of current participant
    datastruc.p1l(i,:) = currsim.p(1,:);
    datastruc.p1r(i,:) = currsim.p(2,:);
    datastruc.p2l(i,:) = currsim.p(3,:);    
    datastruc.p2r(i,:) = currsim.p(4,:);    
    datastruc.p3l(i,:) = currsim.p(5,:);
    datastruc.p3r(i,:) = currsim.p(6,:);
    datastruc.p4l(i,:) = currsim.p(7,:);
    datastruc.p4r(i,:) = currsim.p(8,:);
    datastruc.p5l(i,:) = currsim.p(9,:);
    datastruc.p5r(i,:) = currsim.p(10,:);
    datastruc.p6l(i,:) = currsim.p(11,:);
    datastruc.p6r(i,:) = currsim.p(12,:);
    
    datastruc.Q1l(i,:) = currsim.Q(1,:);
    datastruc.Q1r(i,:) = currsim.Q(2,:);
    datastruc.Q2l(i,:) = currsim.Q(3,:);
    datastruc.Q2r(i,:) = currsim.Q(4,:);
    datastruc.Q3l(i,:) = currsim.Q(5,:);
    datastruc.Q3r(i,:) = currsim.Q(6,:);
    datastruc.Q4l(i,:) = currsim.Q(7,:);
    datastruc.Q4r(i,:) = currsim.Q(8,:);
    datastruc.Q5l(i,:) = currsim.Q(9,:);
    datastruc.Q5r(i,:) = currsim.Q(10,:);
    datastruc.Q6l(i,:) = currsim.Q(11,:);
    datastruc.Q6r(i,:) = currsim.Q(12,:);

    datastruc.PE(i,:) = currsim.PE(1,:);

    
   
    
end

% to save object inbetween as a matlab file
save('comp_model_fit_export\FBOmiss_immediate_Q_values_and_PEs.mat','datastruc');

% Restructure for writing csv
datastruc.p1l_ = transpose(datastruc.p1l);
datastruc.p2l_ = transpose(datastruc.p2l);
datastruc.p3l_ = transpose(datastruc.p3l);
datastruc.p4l_ = transpose(datastruc.p4l);
datastruc.p5l_ = transpose(datastruc.p5l);
datastruc.p6l_ = transpose(datastruc.p6l);

datastruc.p1r_ = transpose(datastruc.p1r);
datastruc.p2r_ = transpose(datastruc.p2r);
datastruc.p3r_ = transpose(datastruc.p3r);
datastruc.p4r_ = transpose(datastruc.p4r);
datastruc.p5r_ = transpose(datastruc.p5r);
datastruc.p6r_ = transpose(datastruc.p6r);

datastruc.Q1l_ = transpose(datastruc.Q1l);
datastruc.Q2l_ = transpose(datastruc.Q2l);
datastruc.Q3l_ = transpose(datastruc.Q3l);
datastruc.Q4l_ = transpose(datastruc.Q4l);
datastruc.Q5l_ = transpose(datastruc.Q5l);
datastruc.Q6l_ = transpose(datastruc.Q6l);

datastruc.Q1r_ = transpose(datastruc.Q1r);
datastruc.Q2r_ = transpose(datastruc.Q2r);
datastruc.Q3r_ = transpose(datastruc.Q3r);
datastruc.Q4r_ = transpose(datastruc.Q4r);
datastruc.Q5r_ = transpose(datastruc.Q5r);
datastruc.Q6r_ = transpose(datastruc.Q6r);

datastruc.PE_ = transpose(datastruc.PE);

p1l = reshape(datastruc.p1l_,[],1);
p2l = reshape(datastruc.p2l_,[],1);
p3l = reshape(datastruc.p3l_,[],1);
p4l = reshape(datastruc.p4l_,[],1);
p5l = reshape(datastruc.p5l_,[],1);
p6l = reshape(datastruc.p6l_,[],1);

Q1l = reshape(datastruc.Q1l_,[],1);
Q2l = reshape(datastruc.Q2l_,[],1);
Q3l = reshape(datastruc.Q3l_,[],1);
Q4l = reshape(datastruc.Q4l_,[],1);
Q5l = reshape(datastruc.Q5l_,[],1);
Q6l = reshape(datastruc.Q6l_,[],1);

p1r = reshape(datastruc.p1r_,[],1);
p2r = reshape(datastruc.p2r_,[],1);
p3r = reshape(datastruc.p3r_,[],1);
p4r = reshape(datastruc.p4r_,[],1);
p5r = reshape(datastruc.p5r_,[],1);
p6r = reshape(datastruc.p6r_,[],1);

Q1r = reshape(datastruc.Q1r_,[],1);
Q2r = reshape(datastruc.Q2r_,[],1);
Q3r = reshape(datastruc.Q3r_,[],1);
Q4r = reshape(datastruc.Q4r_,[],1);
Q5r = reshape(datastruc.Q5r_,[],1);
Q6r = reshape(datastruc.Q6r_,[],1);

PE = reshape(datastruc.PE_,[],1);

% double to cell (necessary for later commands)
p1l = num2cell(p1l);
p2l = num2cell(p2l);
p3l = num2cell(p3l);
p4l = num2cell(p4l);
p5l = num2cell(p5l);
p6l = num2cell(p6l);

Q1l = num2cell(Q1l);
Q2l = num2cell(Q2l);
Q3l = num2cell(Q3l);
Q4l = num2cell(Q4l);
Q5l = num2cell(Q5l);
Q6l = num2cell(Q6l);

p1r = num2cell(p1r);
p2r = num2cell(p2r);
p3r = num2cell(p3r);
p4r = num2cell(p4r);
p5r = num2cell(p5r);
p6r = num2cell(p6r);

Q1r = num2cell(Q1r);
Q2r = num2cell(Q2r);
Q3r = num2cell(Q3r);
Q4r = num2cell(Q4r);
Q5r = num2cell(Q5r);
Q6r = num2cell(Q6r);

PE = num2cell(PE);


% Create column filenames_
%filenames_ = transpose(repelem(filenames,ntrials));
filenames_ = repelem(filenames,ntrials);

% Create table and add variable names
results = {filenames_,p1l,p2l,p3l,p4l,p5l,p6l,p1r,p2r,p3r,p4r,p5r,p6r,Q1l,Q2l,Q3l,Q4l,Q5l,Q6l,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,PE};
results = horzcat(results{:});
results = cell2table(results,...
    'VariableNames',{'filename' 'p1l' 'p2l' 'p3l' 'p4l' 'p5l' 'p6l' 'p1r' 'p2r' 'p3r' 'p4r' 'p5r' 'p6r' 'Q1l' 'Q2l' 'Q3l' 'Q4l' 'Q5l' 'Q6l' 'Q1r' 'Q2r' 'Q3r' 'Q4r' 'Q5r' 'Q6r' 'PE'});

% Export simulated values and PEs as csv
writetable(results,'comp_model_fit_export\FBOmiss_immediate_Q_values_and_PEs.csv');

% Export fits and parameters (BIC + all free parameters of the model)
parameter_export = [filenames, num2cell(fit_BIC.bic(:))];
for k = 1:numel(param_names)
    parameter_export = [parameter_export, num2cell(fit_BIC.(param_names{k})(:))];
end
parameter_export = cell2table(parameter_export, 'VariableNames', [{'filename', 'BIC'}, param_names]);

% export as csv
writetable(parameter_export,'comp_model_fit_export\FBOmiss_learning_parameter.csv');


function sim_data = sim_values_best_model(params, sub_choice, sub_fb, sub_stimuli, sub_feedback_type)
% Trial-wise action values, choice probabilities and prediction errors of
% model 2b for the observed choices and outcomes of one participant.
% The parameter assignment and the update of action values are copied
% verbatim from probfb_fun_model2b.m.

    % Separate learning rates for positive and negative FB
    % + chosen and unchosen action updated

    % set parameters
    alpha_pos = params(1);
    alpha_neg = params(2);
    beta = params(3);

    ntrials = length(sub_choice);

    % action values for each action (left/right) for each of the six stimuli
    Q = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];

    % create empty objects
    sim_data.PE = NaN(1,ntrials);
    sim_data.p = NaN(12,ntrials); % six stimuli with two possible actions each
    sim_data.Q = NaN(12,ntrials);

    for i = 1:ntrials
        if isnan(sub_choice(i)) || isnan(sub_fb(i))
            continue % invalid trial (no response or no feedback)
        end


        % action values of the two options (chosen action first)
        if sub_choice(i) == 1
            chosen_index = sub_stimuli(i)+sub_stimuli(i)-1;
            unchosen_index = sub_stimuli(i)+sub_stimuli(i);
        else
            chosen_index = sub_stimuli(i)+sub_stimuli(i);
            unchosen_index = sub_stimuli(i)+sub_stimuli(i)-1;
        end
        q_i = [Q(chosen_index) Q(unchosen_index)];

        % choice probabilities (softmax)
        p_i = probfb_softmax(q_i,beta);

        % prediction errors of chosen and unchosen action
        PE_c = sub_fb(i) - Q(chosen_index);
        PE_u = 1 - sub_fb(i) - Q(unchosen_index);

        % save probabilities, action values (before the update) and PE
        sim_data.p(chosen_index,i) = p_i(1);
        sim_data.p(unchosen_index,i) = p_i(2);
        sim_data.Q(:,i) = Q;
        sim_data.PE(1,i) = PE_c;

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
end
