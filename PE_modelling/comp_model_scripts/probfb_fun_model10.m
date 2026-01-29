function error = probfb_fun_model10(params,sub_choice,sub_fb,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type)

% Separate learning rates for valence + reward probability + appearance

% set parameters
alpha_90_presented_pos = params(1);
alpha_70_presented_pos = params(2);
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
        p_i = probfb_softmax(q_i,beta);
        
        % calculate prediction error of current trial with expected value
        % of chosen action
        PE_c = sub_fb(i) - Q(chosen_index);
        
        % update action value of chosen action
        
        if sub_stimuli(i) == 1 | sub_stimuli(i) == 2 % if current stimulus is 90%
            
            % if fb is presented and the feedback is positive 
            if sub_feedback_type(i) == 1 & sub_fb(i) == 1 
               
                Q(chosen_index) = Q(chosen_index) + alpha_90_presented_pos*PE_c;
            
            elseif sub_feedback_type(i) == 0 & sub_fb(i) == 0
            % else if feedback is omitted and feedback is negative
                
                Q(chosen_index) = Q(chosen_index) + alpha_90_omitted_neg*PE_c;
            
            elseif sub_feedback_type(i) == 0 & sub_fb(i) == 1
            % else if feedback is omitted and feedback is positive
        
                % use alpha_posomitted to update action value of chosen
                Q(chosen_index) = Q(chosen_index)+ alpha_90_omitted_pos*PE_c;
        
            elseif sub_feedback_type(i) == 1 & sub_fb(i) == 0
            % else if feedback is presented and feedback is negative
        
                % use alpha_negpresented to update action value of chosen
                Q(chosen_index) = Q(chosen_index)+ alpha_90_presented_neg*PE_c;

            end
            
                   
        elseif sub_stimuli(i) == 3 | sub_stimuli(i) == 4 % 70%
            
            % if fb is presented and the feedback is positive 
            if sub_feedback_type(i) == 1 & sub_fb(i) == 1 
               
                Q(chosen_index) = Q(chosen_index) + alpha_70_presented_pos*PE_c;
            
            elseif sub_feedback_type(i) == 0 & sub_fb(i) == 0
            % else if feedback is omitted and feedback is negative
                
                Q(chosen_index) = Q(chosen_index) + alpha_70_omitted_neg*PE_c;
            
            elseif sub_feedback_type(i) == 0 & sub_fb(i) == 1
            % else if feedback is omitted and feedback is positive
        
                % use alpha_posomitted to update action value of chosen
                Q(chosen_index) = Q(chosen_index)+ alpha_70_omitted_pos*PE_c;
        
            elseif sub_feedback_type(i) == 1 & sub_fb(i) == 0
            % else if feedback is presented and feedback is negative
        
                % use alpha_negpresented to update action value of chosen
                Q(chosen_index) = Q(chosen_index)+ alpha_70_presented_neg*PE_c;

            end
            
            
        elseif sub_stimuli(i) == 5 | sub_stimuli(i) == 6 % 50%
            
                    % if fb is presented and the feedback is positive 
            if sub_feedback_type(i) == 1 & sub_fb(i) == 1 
               
                Q(chosen_index) = Q(chosen_index) + alpha_50_presented_pos*PE_c;
            
            elseif sub_feedback_type(i) == 0 & sub_fb(i) == 0
            % else if feedback is omitted and feedback is negative
                
                Q(chosen_index) = Q(chosen_index) + alpha_50_omitted_neg*PE_c;
            
            elseif sub_feedback_type(i) == 0 & sub_fb(i) == 1
            % else if feedback is omitted and feedback is positive
        
                % use alpha_posomitted to update action value of chosen
                Q(chosen_index) = Q(chosen_index)+ alpha_50_omitted_pos*PE_c;
        
            elseif sub_feedback_type(i) == 1 & sub_fb(i) == 0
            % else if feedback is presented and feedback is negative
        
                % use alpha_negpresented to update action value of chosen
                Q(chosen_index) = Q(chosen_index)+ alpha_50_presented_neg*PE_c;

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