% Model 2xb: Fit rl model with two learning rates (sep. for pos/neg FB)
% Only chosen action is updated (but not unchosen)
% PLUS learning rate updates in each trial (Rescorla Wagner - Pearce Hall
% hybrid model) and scaling of impact of PE (kappa)

addpath(genpath('./interim_datasets/'));

% load data
load datastruc_choice; % created with dataprep
load datastruc_not_choice; % created with dataprep
load datastruc_feedback; % created with dataprep
load datastruc_stimuli; % created with dataprep

nsubs = length(datastruc_choice(:,1)); % number of subjects
ntrials = length(datastruc_choice(1,:)); % number of trials

%niter = 100; % number of iterations of fmincon (now globally set in batch file for all functions)

% lower and upper bound for fit (alphas (2), beta, kappa, eta)
LB = [0 0 0 0 0]; % lower bound
UB = [1 1 100 1 1]; % upper bound

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

    n_valid = sum(~isnan(sub_choice) & ~isnan(sub_outcome));

    % loop through iterations
    for j = 1:niter

        % selects random start values
        alpha_pos = rand; % separate learning rates for
        alpha_neg = rand; % positive and negative outcomes
        
        beta = rand*100; % one exploration parameter
        kappa = rand; % scaling parameter of impact of dynamic learning rate
        eta = rand; % volatility scaling of learning rate update
        params = [alpha_pos, alpha_neg, beta, kappa, eta];

        % set options for fmincon
        options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16);

        % call fmincon to minimize output of delivered function (first argument):
        [params, LL] = fmincon(@probfb_fun_model2xb,params,[],[],[],[],LB,UB,[],options,......
            sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli);

        % save optimized parameter of this iteration of fmincon fit
        all_alpha_pos(j) = params(1);
        all_alpha_neg(j) = params(2);
        all_beta(j) = params(3);
        all_kappa(j) = params(4);
        all_eta(j) = params(5);

        % save fit indices
        all_ll(j) = LL;
        all_bic(j) = 2*LL + length(params)*log(n_valid);
    end

    % save best fit according to -LL
    [fit.ll(i),best] = min(all_ll);

    % save parameter values of this fit
    fit.alpha_pos(i) = all_alpha_pos(best);
    fit.alpha_neg(i) = all_alpha_neg(best);
    fit.beta(i) = all_beta(best);
    fit.kappa(i) = all_kappa(best);
    fit.eta(i) = all_eta(best);

    % save best fit according to BIC
    [fit_BIC.bic(i),bestBIC] = min(all_bic);

    % save parameter values of this fit
    fit_BIC.alpha_pos(i) = all_alpha_pos(bestBIC);
    fit_BIC.alpha_neg(i) = all_alpha_neg(bestBIC);
    fit_BIC.beta(i) = all_beta(bestBIC);
    fit_BIC.kappa(i) = all_kappa(bestBIC);
    fit_BIC.eta(i) = all_eta(bestBIC);

    % print to check whether bic and ll get the same fit
    best;
    bestBIC;
end
