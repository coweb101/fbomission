% Model 5b: Fit model with learning rates separately for valence and
% appearance (i.e. a total of 4 learning rates) to simulated data;
% chosen and unchosen action is updated

addpath(genpath('./interim_datasets/'));

% load simulated data
load simulated_data;

nsubs = length(fieldnames(simulated_data)); % number of subjects
ntrials = height(simulated_data.sim_data_01.sim_01); % number of trials

% lower and upper bound for fit (alphas,beta)
LB = [0 0 0 0 0]; % lower bound
UB = [1 1 1 1 100]; % upper bound

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

        n_valid = sum(~isnan(sub_choice) & ~isnan(sub_outcome));
  
            for j = 1:niter % iterations of fmincon

                % selects random start value for alpha and beta    
                alpha_pospresented = rand; 
                alpha_posomitted = rand; 
                alpha_negpresented = rand;
                alpha_negomitted = rand; 

                beta = rand*100;

                params = [alpha_pospresented, alpha_posomitted, ...
                    alpha_negpresented, alpha_negomitted, beta];

                % set options for fmincon
                options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 

                % call fmincon to minimize output of function delivered function (first argument):
                [params, LL] = fmincon(@probfb_fun_model5b,params,[],[],[],[],LB,UB,[],options,......
                    sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type);

                % save optimized parameter of this iteration of fmincon fit
                all_alpha_pospresented(j) = params(1);
                all_alpha_posomitted(j) = params(2);
                all_alpha_negpresented(j) = params(3);
                all_alpha_negomitted(j) = params(4);
                all_beta(j) = params(5);

                % save fit indices
                all_ll(j) = LL;
                all_bic(j) = 2*LL + length(params)*log(n_valid);
            end
        
    % save best fit according to -LL
    [fit.ll(k,i),best] = min(all_ll); %% kann ich hier auch noch zeile indizieren? dann unten auch umsetzen
    
    % save parameter values of this fit
    fit.alpha_pospresented(k,i) = all_alpha_pospresented(best);
    fit.alpha_posomitted(k,i) = all_alpha_posomitted(best);
    fit.alpha_negpresented(k,i) = all_alpha_negpresented(best);
    fit.alpha_negomitted(k,i) = all_alpha_negomitted(best);
    fit.beta(k,i) = all_beta(best);
    
    % save best fit according to BIC
    [fit_BIC.bic(k,i),bestBIC] = min(all_bic);
    
    % save parameter values of this fit
    fit_BIC.alpha_pospresented(k,i) = all_alpha_pospresented(bestBIC);
    fit_BIC.alpha_posomitted(k,i) = all_alpha_posomitted(bestBIC);
    fit_BIC.alpha_negpresented(k,i) = all_alpha_negpresented(bestBIC);
    fit_BIC.alpha_negomitted(k,i) = all_alpha_negomitted(bestBIC);
    
    fit_BIC.beta(k,i) = all_beta(bestBIC);
    
    end        
end
