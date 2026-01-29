% Model 8c: Fit rl model with separate learning rates for stimulus reward 
% probabilities, feedback appearance and choice + update of chosen and
% unchosen

addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep
load datastruc_feedback_type; % created with dataprep

nsubs = length(datastruc_choice(:,1)); % number of subjects
ntrials = length(datastruc_choice(1,:)); % number of trials

% lower and upper bound for fit (alpha,beta)
LB = [0 0 0 0 0 0 0 0 0 0 0 0 0]; % lower bound
UB = [1 1 1 1 1 1 1 1 1 1 1 1 100]; % upper bound

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
        alpha_90_presented = rand; % two learning rates
        alpha_70_presented = rand; % separately for the stimuli 
        alpha_50_presented = rand; % reward probabilities
        alpha_90_omitted = rand; % and feedback appearance
        alpha_70_omitted = rand;  
        alpha_50_omitted = rand;
        alpha_90_presented_unchosen = rand; % and for chosen/unchosen
        alpha_70_presented_unchosen = rand; 
        alpha_50_presented_unchosen = rand; 
        alpha_90_omitted_unchosen = rand; 
        alpha_70_omitted_unchosen = rand;  
        alpha_50_omitted_unchosen = rand;
        
        beta = rand*100; % one exploration parameter
        params = [alpha_90_presented,alpha_70_presented,alpha_50_presented,...
            alpha_90_omitted,alpha_70_omitted,alpha_50_omitted,...
            alpha_90_presented_unchosen,alpha_70_presented_unchosen,alpha_50_presented_unchosen,...
            alpha_90_omitted_unchosen,alpha_70_omitted_unchosen,alpha_50_omitted_unchosen,...
            beta];
        
        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 
        
        % call fmincon to minimize output of delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model8c,params,[],[],[],[],LB,UB,[],options,......
            sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type);
        
        % save optimized parameter of this iteration of fmincon fit
        all_alpha_90_presented(j) = params(1);
        all_alpha_70_presented(j) = params(2);
        all_alpha_50_presented(j) = params(3);
        all_alpha_90_omitted(j) = params(4);
        all_alpha_70_omitted(j) = params(5);
        all_alpha_50_omitted(j) = params(6);    
        all_alpha_90_presented_unchosen(j) = params(7);
        all_alpha_70_presented_unchosen(j) = params(8);
        all_alpha_50_presented_unchosen(j) = params(9);
        all_alpha_90_omitted_unchosen(j) = params(10);
        all_alpha_70_omitted_unchosen(j) = params(11);
        all_alpha_50_omitted_unchosen(j) = params(12);    
        all_beta(j) = params(13);
        
        % save fit indices
        all_ll(j) = LL;
        all_bic(j) = aicbic(-LL,  length(params), ntrials);
    end
    
    % save best fit according to -LL
    [fit.ll(i),best] = min(all_ll);
    
    % save parameter values of this fit
    fit.alpha_90_presented(i) = all_alpha_90_presented(best);
    fit.alpha_70_presented(i) = all_alpha_70_presented(best);
    fit.alpha_50_presented(i) = all_alpha_50_presented(best);
    fit.alpha_90_omitted(i) = all_alpha_90_omitted(best);
    fit.alpha_70_omitted(i) = all_alpha_70_omitted(best);
    fit.alpha_50_omitted(i) = all_alpha_50_omitted(best);
    fit.alpha_90_presented_unchosen(i) = all_alpha_90_presented_unchosen(best);
    fit.alpha_70_presented_unchosen(i) = all_alpha_70_presented_unchosen(best);
    fit.alpha_50_presented_unchosen(i) = all_alpha_50_presented_unchosen(best);
    fit.alpha_90_omitted_unchosen(i) = all_alpha_90_omitted_unchosen(best);
    fit.alpha_70_omitted_unchosen(i) = all_alpha_70_omitted_unchosen(best);
    fit.alpha_50_omitted_unchosen(i) = all_alpha_50_omitted_unchosen(best);
    fit.beta(i) = all_beta(best);
    
    % save best fit according to BIC
    [fit_BIC.bic(i),bestBIC] = min(all_bic);
    
    % save parameter values of this fit
    fit_BIC.alpha_90_presented(i) = all_alpha_90_presented(bestBIC);
    fit_BIC.alpha_70_presented(i) = all_alpha_70_presented(bestBIC);
    fit_BIC.alpha_50_presented(i) = all_alpha_50_presented(bestBIC);
    fit_BIC.alpha_90_omitted(i) = all_alpha_90_omitted(bestBIC);
    fit_BIC.alpha_70_omitted(i) = all_alpha_70_omitted(bestBIC);
    fit_BIC.alpha_50_omitted(i) = all_alpha_50_omitted(bestBIC);
    fit_BIC.alpha_90_presented_unchosen(i) = all_alpha_90_presented_unchosen(bestBIC);
    fit_BIC.alpha_70_presented_unchosen(i) = all_alpha_70_presented_unchosen(bestBIC);
    fit_BIC.alpha_50_presented_unchosen(i) = all_alpha_50_presented_unchosen(bestBIC);
    fit_BIC.alpha_90_omitted_unchosen(i) = all_alpha_90_omitted_unchosen(bestBIC);
    fit_BIC.alpha_70_omitted_unchosen(i) = all_alpha_70_omitted_unchosen(bestBIC);
    fit_BIC.alpha_50_omitted_unchosen(i) = all_alpha_50_omitted_unchosen(bestBIC);
    fit_BIC.beta(i) = all_beta(bestBIC);

end