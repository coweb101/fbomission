function error = probfb_fun_model4(params,sub_choice,sub_fb,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type)

% Separate learning rates for reward and loss context
% Only chosen stimulus is updated (but not unchosen)

% set parameters
alpha_rewardcontext = params(1);
alpha_losscontext = params(2);
beta = params(3);

% initialise

% set action values for each action (right or left) for each of the six
% stimuli (i.e. a total of 12)
Q = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
p = NaN(1,ntrials);
% NaN(M,N) is an M-by-N matrix of NaNs.


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
            chosen_index = sub_stimuli(i)+sub_stimuli(i)-1; % save index for Q-values
        else
            q_i = [Q(sub_stimuli(i)+sub_stimuli(i)) Q(sub_stimuli(i)+sub_stimuli(i)-1)];
            chosen_index = sub_stimuli(i)+sub_stimuli(i); % save index for Q-values
        end

        % probablity of each choice
        p_i = probfb_softmax(q_i,beta); %% hier statt Q, die beiden relevanten Qs
        
        % calculate prediction error of current trial with expected value
        % of chosen action
        PE_c = sub_fb(i) - Q(chosen_index);
        
        % separate learning rates for confirmatory and disconfirmatory
        % trials (i.e., positive and negative fb which is associated with a
        % positive and negative PE, repsectively)
        
        % if fb is presented and the feedback is positive OR feedback is
        % omitted and feedback is negative
        if sub_feedback_type(i) == 1 & sub_fb(i) == 1 | sub_feedback_type(i) == 0 & sub_fb(i) == 0
            
            % update action value of chosen action with alpha_rewardcontext
            Q(chosen_index) = Q(chosen_index) + alpha_rewardcontext*PE_c;
        
        else
            
            % use alpha_losscontext to update action value of chosen
            Q(chosen_index) = Q(chosen_index)+ alpha_losscontext*PE_c;

        end
            
            
        p(1,i) = p_i(1); % save probabilities
        % save first delivered probability from softmax function (which 
        % refers to the stimulus that was chosen because this was delivered
        % as the first argument (Q))
    end
end
p(p<=1e-10) = 1e-5; % avoid 0s or negatives, "underflow"

error = -sum(log(p),'omitnan');