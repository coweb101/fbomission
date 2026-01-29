% Model 10a: Fit model with learning rates separately for outcome
% probabilities, valence, and appearance (i.e. a total of 12 learning
% rates); Chosen and unchosen action is updated

addpath(genpath('./interim_datasets/'));

% load simulated data
load simulated_data;

nsubs = length(fieldnames(simulated_data)); % number of subjects
ntrials = height(simulated_data.sim_data_01.sim_01); % number of trials

% lower and upper bound for fit (alphas,beta)
LB = [0 0 0 0 0 0 0 0 0 0 0 0 0]; % lower bound
UB = [1 1 1 1 1 1 1 1 1 1 1 1 100]; % upper bound

% loop throught subjects
for i = 1:nsubs
    
    i
    % get simulated data of current participant
    data = eval(sprintf('simulated_data.sim_data_%02d', i));
    
    for k = 1:nsim % loop through simulated datasets per participant
        
        data_nsim = eval(sprintf('data.sim_%02d', k));
    
        sub_choice = data_nsim.chosen;
        sub_not_choice = data_nsim.unchosen;
        sub_outcome = data_nsim.feedback;
        sub_stimuli = data_nsim.stim;
        sub_feedback_type = data_nsim.feedback_type;

  
            for j = 1:niter % iterations of fmincon

                % selects random start value for alpha and beta    
                alpha_90_presented_pos = rand; 
                alpha_70_presented_pos = rand; 
                alpha_50_presented_pos = rand;
                alpha_90_omitted_pos = rand; 
                alpha_70_omitted_pos = rand;  
                alpha_50_omitted_pos = rand;
                alpha_90_presented_neg = rand; 
                alpha_70_presented_neg  = rand; 
                alpha_50_presented_neg  = rand;
                alpha_90_omitted_neg  = rand;
                alpha_70_omitted_neg  = rand;  
                alpha_50_omitted_neg  = rand;

                beta = rand*100;

                params = [alpha_90_presented_pos,alpha_70_presented_pos,...
                    alpha_50_presented_pos,alpha_90_omitted_pos,alpha_70_omitted_pos,...
                    alpha_50_omitted_pos,alpha_90_presented_neg,alpha_70_presented_neg,...
                    alpha_50_presented_neg,alpha_90_omitted_neg,alpha_70_omitted_neg,...
                    alpha_50_omitted_neg,beta];

                % set options for fmincon
                options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 

                % call fmincon to minimize output of function delivered function (first argument):
                [params, LL] = fmincon(@probfb_fun_model10a,params,[],[],[],[],LB,UB,[],options,......
                    sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type);

                % save optimized parameter of this iteration of fmincon fit
                all_alpha_90_presented_pos(j) = params(1);
                all_alpha_70_presented_pos(j) = params(2);
                all_alpha_50_presented_pos(j) = params(3);
                all_alpha_90_omitted_pos(j) = params(4);
                all_alpha_70_omitted_pos(j) = params(5);
                all_alpha_50_omitted_pos(j) = params(6);
                all_alpha_90_presented_neg(j) = params(7);
                all_alpha_70_presented_neg(j) = params(8);
                all_alpha_50_presented_neg(j) = params(9);
                all_alpha_90_omitted_neg(j) = params(10);
                all_alpha_70_omitted_neg(j) = params(11);
                all_alpha_50_omitted_neg(j) = params(12);  
                all_beta(j) = params(13);

                % save fit indices
                all_ll(j) = LL;
                all_bic(j) = aicbic(-LL,  length(params), ntrials);
            end
        
    % save best fit according to -LL
    [fit.ll(k,i),best] = min(all_ll); %% kann ich hier auch noch zeile indizieren? dann unten auch umsetzen
    
    % save parameter values of this fit
    fit.alpha_90_presented_pos(k,i) = all_alpha_90_presented_pos(best);
    fit.alpha_70_presented_pos(k,i) = all_alpha_70_presented_pos(best);
    fit.alpha_50_presented_pos(k,i) = all_alpha_50_presented_pos(best);
    fit.alpha_90_omitted_pos(k,i) = all_alpha_90_omitted_pos(best);
    fit.alpha_70_omitted_pos(k,i) = all_alpha_70_omitted_pos(best);
    fit.alpha_50_omitted_pos(k,i) = all_alpha_50_omitted_pos(best);
    fit.alpha_90_presented_neg(k,i) = all_alpha_90_presented_neg(best);
    fit.alpha_70_presented_neg(k,i) = all_alpha_70_presented_neg(best);
    fit.alpha_50_presented_neg(k,i) = all_alpha_50_presented_neg(best);
    fit.alpha_90_omitted_neg(k,i) = all_alpha_90_omitted_neg(best);
    fit.alpha_70_omitted_neg(k,i) = all_alpha_70_omitted_neg(best);
    fit.alpha_50_omitted_neg(k,i) = all_alpha_50_omitted_neg(best);
    fit.beta(k,i) = all_beta(best);
    
    % save best fit according to BIC
    [fit_BIC.bic(k,i),bestBIC] = min(all_bic);
    
    % save parameter values of this fit
    fit_BIC.alpha_90_presented_pos(k,i) = all_alpha_90_presented_pos(bestBIC);
    fit_BIC.alpha_70_presented_pos(k,i) = all_alpha_70_presented_pos(bestBIC);
    fit_BIC.alpha_50_presented_pos(k,i) = all_alpha_50_presented_pos(bestBIC);
    fit_BIC.alpha_90_omitted_pos(k,i) = all_alpha_90_omitted_pos(bestBIC);
    fit_BIC.alpha_70_omitted_pos(k,i) = all_alpha_70_omitted_pos(bestBIC);
    fit_BIC.alpha_50_omitted_pos(k,i) = all_alpha_50_omitted_pos(bestBIC);
    fit_BIC.alpha_90_presented_neg(k,i) = all_alpha_90_presented_neg(bestBIC);
    fit_BIC.alpha_70_presented_neg(k,i) = all_alpha_70_presented_neg(bestBIC);
    fit_BIC.alpha_50_presented_neg(k,i) = all_alpha_50_presented_neg(bestBIC);
    fit_BIC.alpha_90_omitted_neg(k,i) = all_alpha_90_omitted_neg(bestBIC);
    fit_BIC.alpha_70_omitted_neg(k,i) = all_alpha_70_omitted_neg(bestBIC);
    fit_BIC.alpha_50_omitted_neg(k,i) = all_alpha_50_omitted_neg(bestBIC);
    
    fit_BIC.beta(k,i) = all_beta(bestBIC);
    
    end        
end