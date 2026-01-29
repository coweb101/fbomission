% 2. Simulate prediction error for each trial based on best model (model
% 10)

% set directory for model fit
addpath(genpath('./comp_model_fit_export/'));

% load parameters of best BIC fit of model 10a which performed best
load fit_model_10a_bic; 


% set directory for behavioural data (prepared with dataprep
addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep
load datastruc_feedback_type; % created with dataprep
load filenames; %(anonymous; created with with create_correct_response_key;)

ntrials=480;

for i = 1:length(fit_BIC.bic) % loop through participants
    
    i
    
    % retrieve parameters of best fit
    params(1) = fit_BIC.alpha_90_presented_pos(i);
    params(2) = fit_BIC.alpha_70_presented_pos(i);
    params(3) = fit_BIC.alpha_50_presented_pos(i);
    params(4) = fit_BIC.alpha_90_omitted_pos(i);
    params(5) = fit_BIC.alpha_70_omitted_pos(i);
    params(6) = fit_BIC.alpha_50_omitted_pos(i);
    params(7) = fit_BIC.alpha_90_presented_neg(i);
    params(8) = fit_BIC.alpha_70_presented_neg(i);
    params(9) = fit_BIC.alpha_50_presented_neg(i);
    params(10) = fit_BIC.alpha_90_omitted_neg(i);
    params(11) = fit_BIC.alpha_70_omitted_neg(i);
    params(12) = fit_BIC.alpha_50_omitted_neg(i);
    
    params(13) = fit_BIC.beta(i);    
      
    % prepare choice and outcome for simulation function
    sub_choice = datastruc_choice(i,:);
    sub_not_choice = datastruc_not_choice(i,:);
    sub_outcome = datastruc_feedback(i,:);
    sub_stimuli = datastruc_stimuli(i,:);
    sub_feedback_type = datastruc_feedback_type(i,:);
    
    sub_choice = cell2mat(sub_choice);
    sub_not_choice = cell2mat(sub_not_choice);
    sub_outcome = cell2mat(sub_outcome);
    sub_stimuli = cell2mat(sub_stimuli);
    sub_feedback_type = cell2mat(sub_feedback_type);
   
    % simulate values based on parameters of best fit
    currsim = sim_values_model10a(params,sub_choice,sub_not_choice, sub_outcome, sub_stimuli, sub_feedback_type);
    
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

% Export fits and learning rates

% concatenate filenames and ll/alphas/beta
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
    num2cell(transpose(fit_BIC.beta)));
% columns: filenames, BIC, alphas..., beta

parameter_export = cell2table(parameter_export,...
    'VariableNames',{'filename' ...
    'BIC' ...
    'alpha_90_presented_pos'  ...
    'alpha_70_presented_pos'  ...
    'alpha_50_presented_pos'  ...
    'alpha_90_omitted_pos'  ...
    'alpha_70_omitted_pos'  ...
    'alpha_50_omitted_pos'  ...
    'alpha_90_presented_neg'  ...
    'alpha_70_presented_neg'  ...
    'alpha_50_presented_neg'  ...
    'alpha_90_omitted_neg'  ...
    'alpha_70_omitted_neg'  ...
    'alpha_50_omitted_neg'  ...
    'beta'});

% export as csv
writetable(parameter_export,'comp_model_fit_export\FBOmiss_learning_parameter.csv');


function sim_data = sim_values_model10a(params,sub_choice,sub_not_choice, sub_fb, sub_stimuli, sub_feedback_type)

% set parameters

   alpha_90_presented_pos = params(1);
   alpha_70_presented_pos =  params(2);
   alpha_50_presented_pos = params(3);
   alpha_90_omitted_pos = params(4);
   alpha_70_omitted_pos = params(5);
   alpha_50_omitted_pos = params(6);
   alpha_90_presented_neg = params(7);
   alpha_70_presented_neg = params(8);
   alpha_50_presented_neg = params(9);
   alpha_90_omitted_neg = params(10);
   alpha_70_omitted_neg = params(11);
   alpha_50_omitted_neg = params(12);
   
   beta = params(13);

   % save number of trials
   ntrials = length(sub_choice);
    
   % set action values for each action (right or left) for each of the six
   % stimuli (i.e. a total of 12) initially to .5
   Q = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
   p = NaN(1,ntrials);
   % NaN(M,N) is an M-by-N matrix of NaNs.

   % create empty objects
   sim_data.PE = NaN(1,ntrials);
   sim_data.p = NaN(12,ntrials); % six stimuli with two possible actions each
   sim_data.Q = NaN(12,ntrials);
  
   % fit loop
    for i = 1:ntrials
        % trial valid?
        if isnan(sub_choice(i)) || isnan(sub_fb(i))
            p(1,i) = NaN;
        else
            % save in temporal variable q_i the current expectation values of
            % the two choice options in current trial (order depending on
            % choice:        
            if sub_choice(i) == 1
                q_i = [Q(sub_stimuli(i)+sub_stimuli(i)-1) Q(sub_stimuli(i)+sub_stimuli(i))];
                chosen_index = sub_stimuli(i)+sub_stimuli(i)-1; % save index for chosen action Q-value
                unchosen_index = sub_stimuli(i)+sub_stimuli(i);
            else
                q_i = [Q(sub_stimuli(i)+sub_stimuli(i)) Q(sub_stimuli(i)+sub_stimuli(i)-1)];
                chosen_index = sub_stimuli(i)+sub_stimuli(i); % save index for chosen action Q-value
                unchosen_index = sub_stimuli(i)+sub_stimuli(i)-1;
            end

            % estimate probablity of each choice based on softmax function
            p_i = probfb_softmax(q_i,beta);

            % calculate prediction error of current trial with expected value
            % of chosen action
            PE_c = sub_fb(i) - Q(chosen_index);
            % for unchosen option
			PE_u = 1 - sub_fb(i) - Q(unchosen_index);
			
            % save probability of each choice option and action value for this trial

            sim_data.p(chosen_index,i) = p_i(1); % chosen action
            sim_data.p(unchosen_index,i) = p_i(2); % unchosen action
            
            sim_data.Q(:,i) = Q; % save all action values
            
            % save prediction error
            sim_data.PE(1,i) = PE_c;
       

            % update action value of chosen action for next trial

            if sub_stimuli(i) == 1 | sub_stimuli(i) == 2 % if current stimulus is 90%

                % if fb is presented and the feedback is positive 
                if sub_feedback_type(i) == 1 & sub_fb(i) == 1 

                    Q(chosen_index) = Q(chosen_index) + alpha_90_presented_pos*PE_c;
					Q(unchosen_index) = Q(unchosen_index) + alpha_90_presented_pos*PE_u;

                elseif sub_feedback_type(i) == 0 & sub_fb(i) == 0
                % else if feedback is omitted and feedback is negative

                    Q(chosen_index) = Q(chosen_index) + alpha_90_omitted_neg*PE_c;
					Q(unchosen_index) = Q(unchosen_index) + alpha_90_omitted_neg*PE_u;
					
                elseif sub_feedback_type(i) == 0 & sub_fb(i) == 1
                % else if feedback is omitted and feedback is positive

                    % use alpha_posomitted to update action value of chosen
                    Q(chosen_index) = Q(chosen_index)+ alpha_90_omitted_pos*PE_c;
					Q(unchosen_index) = Q(unchosen_index)+ alpha_90_omitted_pos*PE_u;
					
                elseif sub_feedback_type(i) == 1 & sub_fb(i) == 0
                % else if feedback is presented and feedback is negative

                    % use alpha_negpresented to update action value of chosen
                    Q(chosen_index) = Q(chosen_index)+ alpha_90_presented_neg*PE_c;
					Q(unchosen_index) = Q(unchosen_index)+ alpha_90_presented_neg*PE_u;

                end


            elseif sub_stimuli(i) == 3 | sub_stimuli(i) == 4 % 70%

                % if fb is presented and the feedback is positive 
                if sub_feedback_type(i) == 1 & sub_fb(i) == 1 

                    Q(chosen_index) = Q(chosen_index) + alpha_70_presented_pos*PE_c;
					Q(unchosen_index) = Q(unchosen_index) + alpha_70_presented_pos*PE_u;

                elseif sub_feedback_type(i) == 0 & sub_fb(i) == 0
                % else if feedback is omitted and feedback is negative

                    Q(chosen_index) = Q(chosen_index) + alpha_70_omitted_neg*PE_c;
					Q(unchosen_index) = Q(unchosen_index) + alpha_70_omitted_neg*PE_u;

                elseif sub_feedback_type(i) == 0 & sub_fb(i) == 1
                % else if feedback is omitted and feedback is positive

                    % use alpha_posomitted to update action value of chosen
                    Q(chosen_index) = Q(chosen_index)+ alpha_70_omitted_pos*PE_c;
					Q(unchosen_index) = Q(unchosen_index)+ alpha_70_omitted_pos*PE_u;

                elseif sub_feedback_type(i) == 1 & sub_fb(i) == 0
                % else if feedback is presented and feedback is negative

                    % use alpha_negpresented to update action value of chosen
                    Q(chosen_index) = Q(chosen_index)+ alpha_70_presented_neg*PE_c;
					Q(unchosen_index) = Q(unchosen_index)+ alpha_70_presented_neg*PE_u;

                end


            elseif sub_stimuli(i) == 5 | sub_stimuli(i) == 6 % 50%

                        % if fb is presented and the feedback is positive 
                if sub_feedback_type(i) == 1 & sub_fb(i) == 1 

                    Q(chosen_index) = Q(chosen_index) + alpha_50_presented_pos*PE_c;
					Q(unchosen_index) = Q(unchosen_index) + alpha_50_presented_pos*PE_u;

                elseif sub_feedback_type(i) == 0 & sub_fb(i) == 0
                % else if feedback is omitted and feedback is negative

                    Q(chosen_index) = Q(chosen_index) + alpha_50_omitted_neg*PE_c;
					Q(unchosen_index) = Q(unchosen_index) + alpha_50_omitted_neg*PE_u;

                elseif sub_feedback_type(i) == 0 & sub_fb(i) == 1
                % else if feedback is omitted and feedback is positive

                    % use alpha_posomitted to update action value of chosen
                    Q(chosen_index) = Q(chosen_index)+ alpha_50_omitted_pos*PE_c;
					Q(unchosen_index) = Q(unchosen_index)+ alpha_50_omitted_pos*PE_u;

                elseif sub_feedback_type(i) == 1 & sub_fb(i) == 0
                % else if feedback is presented and feedback is negative

                    % use alpha_negpresented to update action value of chosen
                    Q(chosen_index) = Q(chosen_index)+ alpha_50_presented_neg*PE_c;
					Q(unchosen_index) = Q(unchosen_index)+ alpha_50_presented_neg*PE_u;

                end

            end
        end
    end
end

function [probs] = probfb_softmax(Q,beta)

% INPUT:
% Q is the two estimated action values, e.g. [0.5 0.5]
% beta is the inverse temperature
% OUTPUT:
% the probability of the choices

probs = (exp(beta*Q))/(exp(beta*Q(1))+exp(beta*Q(2)));
% Softmax/Boltzmann Function zur Berechnung der W'keit der Auswahl der Akt.
end
