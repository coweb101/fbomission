% Model 10aac: Fit rl model with separate learning rates for valence + reward 
% probability + appearance + choice; plus update of unchosen

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
LB = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % lower bound
UB = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 100]; % upper bound

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
        alpha_90_presented_pos_unchosen = rand; % and separately for chosen
        alpha_70_presented_pos_unchosen = rand; % and unchosen action 
        alpha_50_presented_pos_unchosen = rand;
        alpha_90_omitted_pos_unchosen = rand; 
        alpha_70_omitted_pos_unchosen = rand;  
        alpha_50_omitted_pos_unchosen = rand;
        alpha_90_presented_neg_unchosen = rand; 
        alpha_70_presented_neg_unchosen  = rand; 
        alpha_50_presented_neg_unchosen  = rand;
        alpha_90_omitted_neg_unchosen  = rand;
        alpha_70_omitted_neg_unchosen  = rand;  
        alpha_50_omitted_neg_unchosen  = rand;
        
        beta = rand*100; % one exploration parameter
        
        params = [alpha_90_presented_pos,alpha_70_presented_pos,...
            alpha_50_presented_pos,alpha_90_omitted_pos,alpha_70_omitted_pos,...
            alpha_50_omitted_pos,alpha_90_presented_neg,alpha_70_presented_neg,...
            alpha_50_presented_neg,alpha_90_omitted_neg,alpha_70_omitted_neg,...
            alpha_50_omitted_neg,...
            alpha_90_presented_pos_unchosen,alpha_70_presented_pos_unchosen,...
            alpha_50_presented_pos_unchosen,alpha_90_omitted_pos_unchosen,alpha_70_omitted_pos_unchosen,...
            alpha_50_omitted_pos_unchosen,alpha_90_presented_neg_unchosen,alpha_70_presented_neg_unchosen,...
            alpha_50_presented_neg_unchosen,alpha_90_omitted_neg_unchosen,alpha_70_omitted_neg_unchosen,...
            alpha_50_omitted_neg_unchosen,...
            beta];
        
        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 
        
        % call fmincon to minimize output of function delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model10aac,params,[],[],[],[],LB,UB,[],options,......
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
        all_alpha_90_presented_pos_unchosen(j) = params(13);
        all_alpha_70_presented_pos_unchosen(j) = params(14);
        all_alpha_50_presented_pos_unchosen(j) = params(15);
        all_alpha_90_omitted_pos_unchosen(j) = params(16);
        all_alpha_70_omitted_pos_unchosen(j) = params(17);
        all_alpha_50_omitted_pos_unchosen(j) = params(18);
        all_alpha_90_presented_neg_unchosen(j) = params(19);
        all_alpha_70_presented_neg_unchosen(j) = params(20);
        all_alpha_50_presented_neg_unchosen(j) = params(21);
        all_alpha_90_omitted_neg_unchosen(j) = params(22);
        all_alpha_70_omitted_neg_unchosen(j) = params(23);
        all_alpha_50_omitted_neg_unchosen(j) = params(24);  
        all_beta(j) = params(25);
        
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
    fit.alpha_90_presented_pos_unchosen(i) = all_alpha_90_presented_pos_unchosen(best);
    fit.alpha_70_presented_pos_unchosen(i) = all_alpha_70_presented_pos_unchosen(best);
    fit.alpha_50_presented_pos_unchosen(i) = all_alpha_50_presented_pos_unchosen(best);
    fit.alpha_90_omitted_pos_unchosen(i) = all_alpha_90_omitted_pos_unchosen(best);
    fit.alpha_70_omitted_pos_unchosen(i) = all_alpha_70_omitted_pos_unchosen(best);
    fit.alpha_50_omitted_pos_unchosen(i) = all_alpha_50_omitted_pos_unchosen(best);
    fit.alpha_90_presented_neg_unchosen(i) = all_alpha_90_presented_neg_unchosen(best);
    fit.alpha_70_presented_neg_unchosen(i) = all_alpha_70_presented_neg_unchosen(best);
    fit.alpha_50_presented_neg_unchosen(i) = all_alpha_50_presented_neg_unchosen(best);
    fit.alpha_90_omitted_neg_unchosen(i) = all_alpha_90_omitted_neg_unchosen(best);
    fit.alpha_70_omitted_neg_unchosen(i) = all_alpha_70_omitted_neg_unchosen(best);
    fit.alpha_50_omitted_neg_unchosen(i) = all_alpha_50_omitted_neg_unchosen(best);
    fit.beta(i) = all_beta(best);
    
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
    fit_BIC.alpha_90_presented_pos_unchosen(i) = all_alpha_90_presented_pos_unchosen(bestBIC);
    fit_BIC.alpha_70_presented_pos_unchosen(i) = all_alpha_70_presented_pos_unchosen(bestBIC);
    fit_BIC.alpha_50_presented_pos_unchosen(i) = all_alpha_50_presented_pos_unchosen(bestBIC);
    fit_BIC.alpha_90_omitted_pos_unchosen(i) = all_alpha_90_omitted_pos_unchosen(bestBIC);
    fit_BIC.alpha_70_omitted_pos_unchosen(i) = all_alpha_70_omitted_pos_unchosen(bestBIC);
    fit_BIC.alpha_50_omitted_pos_unchosen(i) = all_alpha_50_omitted_pos_unchosen(bestBIC);
    fit_BIC.alpha_90_presented_neg_unchosen(i) = all_alpha_90_presented_neg_unchosen(bestBIC);
    fit_BIC.alpha_70_presented_neg_unchosen(i) = all_alpha_70_presented_neg_unchosen(bestBIC);
    fit_BIC.alpha_50_presented_neg_unchosen(i) = all_alpha_50_presented_neg_unchosen(bestBIC);
    fit_BIC.alpha_90_omitted_neg_unchosen(i) = all_alpha_90_omitted_neg_unchosen(bestBIC);
    fit_BIC.alpha_70_omitted_neg_unchosen(i) = all_alpha_70_omitted_neg_unchosen(bestBIC);
    fit_BIC.alpha_50_omitted_neg_unchosen(i) = all_alpha_50_omitted_neg_unchosen(bestBIC);
    
    fit_BIC.beta(i) = all_beta(bestBIC);
    
end
