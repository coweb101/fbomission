function error = probfb_fun_model8c(params,sub_choice,sub_fb,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type)

% Separate learning rates for appearance, reward probability, choice
% + update of unchosen action

% set parameters
alpha_90_presented = params(1);
alpha_70_presented = params(2);
alpha_50_presented = params(3);
alpha_90_omitted = params(4);
alpha_70_omitted = params(5);
alpha_50_omitted = params(6);
alpha_90_presented_unchosen = params(7);
alpha_70_presented_unchosen = params(8);
alpha_50_presented_unchosen = params(9);
alpha_90_omitted_unchosen = params(10);
alpha_70_omitted_unchosen = params(11);
alpha_50_omitted_unchosen = params(12);
beta = params(13);

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
            unchosen_index = sub_stimuli(i)+sub_stimuli(i); % save index for Q-value of unchosen stimulus
        else
            q_i = [Q(sub_stimuli(i)+sub_stimuli(i)) Q(sub_stimuli(i)+sub_stimuli(i)-1)];
            chosen_index = sub_stimuli(i)+sub_stimuli(i); % save index for Q-values
            unchosen_index = sub_stimuli(i)+sub_stimuli(i)-1; % save index for Q-value of unchosen stimulus
         end

        % probablity of each choice
        p_i = probfb_softmax(q_i,beta); %% hier statt Q, die beiden relevanten Qs
        
        % calculate prediction error of current trial with expected value
        % of chosen action
        PE_c = sub_fb(i) - Q(chosen_index);
        PE_u = 1 - sub_fb(i) - Q(unchosen_index);
        
        % update action value of chosen action
        
        if sub_stimuli(i) == 1 | sub_stimuli(i) == 2 % if current stimulus is 90%
            
            if sub_feedback_type(i) == 1 % if fb is presented
               
                Q(chosen_index) = Q(chosen_index) + alpha_90_presented*PE_c;
                Q(unchosen_index) = Q(unchosen_index) + alpha_90_presented_unchosen*PE_u;
            elseif sub_feedback_type(i) == 0 % if fb is omitted
                
                Q(chosen_index) = Q(chosen_index) + alpha_90_omitted*PE_c;
                Q(unchosen_index) = Q(unchosen_index) + alpha_90_omitted_unchosen*PE_u;
            end
            
        elseif sub_stimuli(i) == 3 | sub_stimuli(i) == 4 % 70%
            
            if sub_feedback_type(i) == 1 % if fb is presented
                
                Q(chosen_index) = Q(chosen_index) + alpha_70_presented*PE_c;
                Q(unchosen_index) = Q(unchosen_index) + alpha_70_presented_unchosen*PE_u;
            elseif sub_feedback_type(i) == 0 % if fb is omitted
                
                Q(chosen_index) = Q(chosen_index) + alpha_70_omitted*PE_c;
                Q(unchosen_index) = Q(unchosen_index) + alpha_70_omitted_unchosen*PE_u;
            end
            
        elseif sub_stimuli(i) == 5 | sub_stimuli(i) == 6 % 50%
            
            if sub_feedback_type(i) == 1 % if fb is presented
                
                Q(chosen_index) = Q(chosen_index) + alpha_50_presented*PE_c;
                Q(unchosen_index) = Q(unchosen_index) + alpha_50_presented_unchosen*PE_u;
            
            elseif sub_feedback_type(i) == 0 % if fb is omitted
                
                Q(chosen_index) = Q(chosen_index) + alpha_50_omitted*PE_c;
                Q(unchosen_index) = Q(unchosen_index) + alpha_50_omitted_unchosen*PE_u;
                
            end
        end
        
        p(1,i) = p_i(1); % save probabilities
        % save first delivered probability from softmax function (which 
        % refers to the stimulus that was chosen because this was delivered
        % as the first argument (Q))
    end
end
p(p<=1e-10) = 1e-5; % avoid 0s or negatives, "underflow"

error = -sum(log(p),'omitnan');