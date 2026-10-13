function error = probfb_fun_model2cb(params,sub_choice,sub_fb,sub_not_choice,ntrials, sub_stimuli)

% Model 2cb: Separate learning rates for positive and negative FB (valence)
% separate learning rates for chosen and unchosen action
% + update of chosen and unchosen action value
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% hybrid model) and scaling of impact of PE (kappa)

% set parameters
alpha_pos = params(1);
alpha_neg = params(2);
alpha_pos_unchosen = params(3);
alpha_neg_unchosen = params(4);
beta = params(5);
kappa = params(6);
eta = params(7);

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
        % choice):
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
        % of chosen and unchosen action
        PE_c = sub_fb(i) - Q(chosen_index);
        PE_u = 1 - sub_fb(i) - Q(unchosen_index);

        % separate learning rates for confirmatory and disconfirmatory
        % trials (i.e., positive and negative fb which is associated with a
        % positive and negative PE, repsectively)

        if sub_fb(i) == 1 % if fb is positive
            
            % update action value of chosen action with scaled alpha_pos
            Q(chosen_index) = Q(chosen_index) + kappa*alpha_pos*PE_c;
            Q(unchosen_index) = Q(unchosen_index) + kappa*alpha_pos_unchosen*PE_u;
            % update learning rates for the next trial of this type (Pearce-Hall)
            alpha_pos = eta*abs(PE_c) + (1-eta)*alpha_pos;
            alpha_pos_unchosen = eta*abs(PE_u) + (1-eta)*alpha_pos_unchosen;
        
        elseif sub_fb(i) == 0 % if fb is negative
            
            % update action value of chosen action with scaled alpha_neg
            Q(chosen_index) = Q(chosen_index) + kappa*alpha_neg*PE_c;
            Q(unchosen_index) = Q(unchosen_index) + kappa*alpha_neg_unchosen*PE_u;
            % update learning rates for the next trial of this type (Pearce-Hall)
            alpha_neg = eta*abs(PE_c) + (1-eta)*alpha_neg;
            alpha_neg_unchosen = eta*abs(PE_u) + (1-eta)*alpha_neg_unchosen;
        
        end


        p(1,i) = p_i(1); % save probabilities
        % save first delivered probability from softmax function (which
        % refers to the stimulus that was chosen because this was delivered
        % as the first argument (Q))
    end
end
p(p<1e-10) = 1e-10; % avoid 0s or negatives, "underflow"

error = -sum(log(p),'omitnan');