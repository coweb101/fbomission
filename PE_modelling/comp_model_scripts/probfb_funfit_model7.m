% Model 7: Fit rl model with two learning rates per stimulus reward 
% probabilities (i.e. a total of 6 learning rates) separately for valence
% but across feedback type/context

addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep

nsubs = length(datastruc_choice(:,1)); % number of subjects
ntrials = length(datastruc_choice(1,:)); % number of trials

% niter = 100; % number of iterations of fmincon(now globally set in batch file for all functions)

% lower and upper bound for fit (alpha,beta)
LB = [0 0 0 0 0 0 0]; % lower bound
UB = [1 1 1 1 1 1 100]; % upper bound

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
        
        beta = rand*100; % one exploration parameter
        params = [alpha_90_pos,alpha_70_pos,alpha_50_pos,alpha_90_neg,alpha_70_neg,alpha_50_neg,beta];
        
        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16); 
        
        % call fmincon to minimize output of delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model7,params,[],[],[],[],LB,UB,[],options,......
            sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli);
        
        % save optimized parameter of this iteration of fmincon fit
        all_alpha_90_pos(j) = params(1);
        all_alpha_70_pos(j) = params(2);
        all_alpha_50_pos(j) = params(3);
        all_alpha_90_neg(j) = params(4);
        all_alpha_70_neg(j) = params(5);
        all_alpha_50_neg(j) = params(6);    
        all_beta(j) = params(7);
        
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
    fit_BIC.beta(i) = all_beta(bestBIC);
    
    % print to check whether bic and ll get the same fit
    best;
    bestBIC;
end

%save('probfb_funfit_results_model1','fit');