% Model 7c: Fit model with two learning rates per stimulus reward 
% probabilities (i.e. a total of 6 learning rates) separately for valence
% and separately for chosen and unchosen action (i.e. a total of 12
% learning rates) + update of unchosen action 

addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep

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
    
    sub_choice = cell2mat(sub_choice);
    sub_not_choice = cell2mat(sub_not_choice);
    sub_outcome = cell2mat(sub_outcome);
    sub_stimuli = cell2mat(sub_stimuli);
    
    % loop through iterations
    for j = 1:niter
        
        % selects random start value for alpha and beta    
        alpha_90_pos = rand; % two learning rates
        alpha_70_pos = rand; % separately for the stimuli 
        alpha_50_pos = rand; % reward probabilities
        alpha_90_neg = rand; % and feedback valence
        alpha_70_neg = rand;  
        alpha_50_neg = rand;
        alpha_90_pos_unchosen = rand; % and separately for unchosen action
        alpha_70_pos_unchosen = rand; 
        alpha_50_pos_unchosen = rand; 
        alpha_90_neg_unchosen = rand;
        alpha_70_neg_unchosen = rand;  
        alpha_50_neg_unchosen = rand;
        
        beta = rand*100; % one exploration parameter
        params = [alpha_90_pos,alpha_70_pos,alpha_50_pos,alpha_90_neg,alpha_70_neg,alpha_50_neg,...
        alpha_90_pos_unchosen,alpha_70_pos_unchosen,alpha_50_pos_unchosen,alpha_90_neg_unchosen,alpha_70_neg_unchosen,alpha_50_neg_unchosen, beta];
        
        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 
        
        % call fmincon to minimize output of delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model7c,params,[],[],[],[],LB,UB,[],options,......
            sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli);
        
        % save optimized parameter of this iteration of fmincon fit
        all_alpha_90_pos(j) = params(1);
        all_alpha_70_pos(j) = params(2);
        all_alpha_50_pos(j) = params(3);
        all_alpha_90_neg(j) = params(4);
        all_alpha_70_neg(j) = params(5);
        all_alpha_50_neg(j) = params(6);  
        all_alpha_90_pos_unchosen(j) = params(7);
        all_alpha_70_pos_unchosen(j) = params(8);
        all_alpha_50_pos_unchosen(j) = params(9);
        all_alpha_90_neg_unchosen(j) = params(10);
        all_alpha_70_neg_unchosen(j) = params(11);
        all_alpha_50_neg_unchosen(j) = params(12);  
        all_beta(j) = params(13);
        
        % save fit indices
        all_ll(j) = LL;
        all_bic(j) = aicbic(-LL,  length(params), ntrials);
    end
    
    % save best fit according to -LL
    [fit.ll(i),best] = min(all_ll);
    
    % save parameter values of this fit
    fit.alpha_90_pos(i) = all_alpha_90_pos(best);
    fit.alpha_70_pos(i) = all_alpha_70_pos(best);
    fit.alpha_50_pos(i) = all_alpha_50_pos(best);
    fit.alpha_90_neg(i) = all_alpha_90_neg(best);
    fit.alpha_70_neg(i) = all_alpha_70_neg(best);
    fit.alpha_50_neg(i) = all_alpha_50_neg(best);
    fit.alpha_90_pos_unchosen(i) = all_alpha_90_pos_unchosen(best);
    fit.alpha_70_pos_unchosen(i) = all_alpha_70_pos_unchosen(best);
    fit.alpha_50_pos_unchosen(i) = all_alpha_50_pos_unchosen(best);
    fit.alpha_90_neg_unchosen(i) = all_alpha_90_neg_unchosen(best);
    fit.alpha_70_neg_unchosen(i) = all_alpha_70_neg_unchosen(best);
    fit.alpha_50_neg_unchosen(i) = all_alpha_50_neg_unchosen(best);
    fit.beta(i) = all_beta(best);
    
    % save best fit according to BIC
    [fit_BIC.bic(i),bestBIC] = min(all_bic);
    
    % save parameter values of this fit
    fit_BIC.alpha_90_pos(i) = all_alpha_90_pos(bestBIC);
    fit_BIC.alpha_70_pos(i) = all_alpha_70_pos(bestBIC);
    fit_BIC.alpha_50_pos(i) = all_alpha_50_pos(bestBIC);
    fit_BIC.alpha_90_neg(i) = all_alpha_90_neg(bestBIC);
    fit_BIC.alpha_70_neg(i) = all_alpha_70_neg(bestBIC);
    fit_BIC.alpha_50_neg(i) = all_alpha_50_neg(bestBIC);
    fit_BIC.alpha_90_pos_unchosen(i) = all_alpha_90_pos_unchosen(bestBIC);
    fit_BIC.alpha_70_pos_unchosen(i) = all_alpha_70_pos_unchosen(bestBIC);
    fit_BIC.alpha_50_pos_unchosen(i) = all_alpha_50_pos_unchosen(bestBIC);
    fit_BIC.alpha_90_neg_unchosen(i) = all_alpha_90_neg_unchosen(bestBIC);
    fit_BIC.alpha_70_neg_unchosen(i) = all_alpha_70_neg_unchosen(bestBIC);
    fit_BIC.alpha_50_neg_unchosen(i) = all_alpha_50_neg_unchosen(bestBIC);
    fit_BIC.beta(i) = all_beta(bestBIC);
    
end
