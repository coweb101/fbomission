% Model 10c: Fit rl model with learning rates separately for reward 
% probabilities, and each possible feedback type (e.g. positive presented, 
% positive omitted % etc.)
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% Hybrid model) and scaling of impact of PE

addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep
load datastruc_feedback_type; % created with dataprep

nsubs = length(datastruc_choice(:,1)); % number of subjects
ntrials = length(datastruc_choice(1,:)); % number of trials

% lower and upper bound for fit (alphas (12),beta, kappa, eta)
LB = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % lower bound
UB = [1 1 1 1 1 1 1 1 1 1 1 1 100 1 1]; % upper bound

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
        alpha_90_presented_pos = rand; % 4 learning rates
        alpha_70_presented_pos = rand; % separately for the stimuli 
        alpha_50_presented_pos = rand; % reward probabilities
        alpha_90_omitted_pos = rand; 
        alpha_70_omitted_pos = rand;  
        alpha_50_omitted_pos = rand;
        alpha_90_presented_neg = rand; 
        alpha_70_presented_neg  = rand; 
        alpha_50_presented_neg  = rand;
        alpha_90_omitted_neg  = rand;
        alpha_70_omitted_neg  = rand;  
        alpha_50_omitted_neg  = rand;
        
        beta = rand*100; % one exploration parameter
        
        kappa = rand; % scaling parameter of impact of dynamic learning rate
        eta = rand; % volatility scaling of learning rate update
        
        params = [alpha_90_presented_pos,alpha_70_presented_pos,...
            alpha_50_presented_pos,alpha_90_omitted_pos,alpha_70_omitted_pos,...
            alpha_50_omitted_pos,alpha_90_presented_neg,alpha_70_presented_neg,...
            alpha_50_presented_neg,alpha_90_omitted_neg,alpha_70_omitted_neg,...
            alpha_50_omitted_neg,beta,kappa,eta];
        
        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 
        
        % call fmincon to minimize output of function delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model10c,params,[],[],[],[],LB,UB,[],options,......
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
        all_kappa(j) = params(14);
        all_eta(j) = params(15);
        
        % save fit indices
        all_ll(j) = LL;
        all_bic(j) = aicbic(-LL,  length(params), ntrials);
    end
    
    % save best fit according to -LL
    [fit.ll(i),best] = min(all_ll);
    
    % save parameter values of this fit
    fit.alpha_90_presented_pos(i) = all_alpha_90_presented_pos(best);
    fit.alpha_70_presented_pos(i) = all_alpha_70_presented_pos(best);
    fit.alpha_50_presented_pos(i) = all_alpha_50_presented_pos(best);
    fit.alpha_90_omitted_pos(i) = all_alpha_90_omitted_pos(best);
    fit.alpha_70_omitted_pos(i) = all_alpha_70_omitted_pos(best);
    fit.alpha_50_omitted_pos(i) = all_alpha_50_omitted_pos(best);
    fit.alpha_90_presented_neg(i) = all_alpha_90_presented_neg(best);
    fit.alpha_70_presented_neg(i) = all_alpha_70_presented_neg(best);
    fit.alpha_50_presented_neg(i) = all_alpha_50_presented_neg(best);
    fit.alpha_90_omitted_neg(i) = all_alpha_90_omitted_neg(best);
    fit.alpha_70_omitted_neg(i) = all_alpha_70_omitted_neg(best);
    fit.alpha_50_omitted_neg(i) = all_alpha_50_omitted_neg(best);
    fit.beta(i) = all_beta(best);
    fit.kappa(i) = all_kappa(best);
    fit.eta(i) = all_eta(best);
    
    % save best fit according to BIC
    [fit_BIC.bic(i),bestBIC] = min(all_bic);
    
    % save parameter values of this fit

    fit_BIC.alpha_90_presented_pos(i) = all_alpha_90_presented_pos(bestBIC);
    fit_BIC.alpha_70_presented_pos(i) = all_alpha_70_presented_pos(bestBIC);
    fit_BIC.alpha_50_presented_pos(i) = all_alpha_50_presented_pos(bestBIC);
    fit_BIC.alpha_90_omitted_pos(i) = all_alpha_90_omitted_pos(bestBIC);
    fit_BIC.alpha_70_omitted_pos(i) = all_alpha_70_omitted_pos(bestBIC);
    fit_BIC.alpha_50_omitted_pos(i) = all_alpha_50_omitted_pos(bestBIC);
    fit_BIC.alpha_90_presented_neg(i) = all_alpha_90_presented_neg(bestBIC);
    fit_BIC.alpha_70_presented_neg(i) = all_alpha_70_presented_neg(bestBIC);
    fit_BIC.alpha_50_presented_neg(i) = all_alpha_50_presented_neg(bestBIC);
    fit_BIC.alpha_90_omitted_neg(i) = all_alpha_90_omitted_neg(bestBIC);
    fit_BIC.alpha_70_omitted_neg(i) = all_alpha_70_omitted_neg(bestBIC);
    fit_BIC.alpha_50_omitted_neg(i) = all_alpha_50_omitted_neg(bestBIC);
    
    fit_BIC.beta(i) = all_beta(bestBIC);
    fit_BIC.kappa(i) = all_kappa(bestBIC);
    fit_BIC.eta(i) = all_eta(bestBIC);
    
    
    % print to check whether bic and ll get the same fit
    best;
    bestBIC;
end

%save('probfb_funfit_results_model1','fit');