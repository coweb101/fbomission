% Model 5c: Fit model with eight learning rates for each of the possible 
% four feedback types (i.e. positive presented, positive omitted, negative 
% presented, negative omitted, which is inherently also separately for 
% contexts) and chosen and unchosen action
% + update of unchosen action

addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep
load datastruc_stimuli; % created with dataprep
load datastruc_feedback_type; % created with dataprep

nsubs = length(datastruc_choice(:,1)); % number of subjects
ntrials = length(datastruc_choice(1,:)); % number of trials

% lower and upper bound for fit (alphas,beta)
LB = [0 0 0 0 0 0 0 0 0]; % lower bound
UB = [1 1 1 1 1 1 1 1 100]; % upper bound

% loop throught subjects
for i = 1:nsubs
    i
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
    
    % loop through iterations
    for j = 1:niter
        
        % selects random start value for alpha and beta    
        alpha_pospresented = rand; % separate learning rates for
        alpha_posomitted = rand; % each possible feedback type
        alpha_negpresented = rand; 
        alpha_negomitted = rand;
        alpha_pospresented_unchosen = rand; % separate learning rates for
        alpha_posomitted_unchosen = rand; % each possible feedback type
        alpha_negpresented_unchosen = rand; 
        alpha_negomitted_unchosen = rand;
        
        beta = rand*100; % one exploration parameter
        params = [alpha_pospresented, alpha_posomitted, alpha_negpresented, alpha_negomitted, alpha_pospresented_unchosen, alpha_posomitted_unchosen, alpha_negpresented_unchosen, alpha_negomitted_unchosen, beta];
        
        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 
        
        % call fmincon to minimize output of delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model5c,params,[],[],[],[],LB,UB,[],options,......
            sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type);
        
        % save optimized parameter of this iteration of fmincon fit
        all_alpha_pospresented(j) = params(1);
        all_alpha_posomitted(j) = params(2);
        all_alpha_negpresented(j) = params(3);
        all_alpha_negomitted(j) = params(4);
        all_alpha_pospresented_unchosen(j) = params(5);
        all_alpha_posomitted_unchosen(j) = params(6);
        all_alpha_negpresented_unchosen(j) = params(7);
        all_alpha_negomitted_unchosen(j) = params(8);
        all_beta(j) = params(9);
        
        % save fit indices
        all_ll(j) = LL;
        all_bic(j) = aicbic(-LL,  length(params), ntrials);
    end
    
    % save best fit according to -LL
    [fit.ll(i),best] = min(all_ll);
    
    % save parameter values of this fit
    fit.alpha_pospresented(i) = all_alpha_pospresented(best);
    fit.alpha_posomitted(i) = all_alpha_posomitted(best);
    fit.alpha_negpresented(i) = all_alpha_negpresented(best);
    fit.alpha_negomitted(i) = all_alpha_negomitted(best);
    fit.alpha_pospresented_unchosen(i) = all_alpha_pospresented_unchosen(best);
    fit.alpha_posomitted_unchosen(i) = all_alpha_posomitted_unchosen(best);
    fit.alpha_negpresented_unchosen(i) = all_alpha_negpresented_unchosen(best);
    fit.alpha_negomitted_unchosen(i) = all_alpha_negomitted_unchosen(best);
    fit.beta(i) = all_beta(best);
    
    % save best fit according to BIC
    [fit_BIC.bic(i),bestBIC] = min(all_bic);
    
    % save parameter values of this fit
    
    fit_BIC.alpha_pospresented(i) = all_alpha_pospresented(bestBIC);
    fit_BIC.alpha_posomitted(i) = all_alpha_posomitted(bestBIC);
    fit_BIC.alpha_negpresented(i) = all_alpha_negpresented(bestBIC);
    fit_BIC.alpha_negomitted(i) = all_alpha_negomitted(bestBIC);
    fit_BIC.alpha_pospresented_unchosen(i) = all_alpha_pospresented_unchosen(bestBIC);
    fit_BIC.alpha_posomitted_unchosen(i) = all_alpha_posomitted_unchosen(bestBIC);
    fit_BIC.alpha_negpresented_unchosen(i) = all_alpha_negpresented_unchosen(bestBIC);
    fit_BIC.alpha_negomitted_unchosen(i) = all_alpha_negomitted_unchosen(bestBIC);
    fit_BIC.beta(i) = all_beta(bestBIC);

end