% Parameter recovery: fit model 10ab (lowest AIC) to the simulated data
% Model 10ab: learning rates per reward probability x valence x appearance, exponential decay over presentations of the current stimulus (half-life lambda) (+ update of chosen and unchosen option)
% (start values, bounds and fmincon call copied from probfb_funfit_model10ab.m)

addpath(genpath('./interim_datasets/'));

% remove fits of the empirical data from the workspace (same variable names)
clear fit fit_BIC

% load simulated data (created with probfb_simulation_data_parameter_recovery)
load simulated_data_AIC;

nsubs = length(fieldnames(simulated_data)); % number of participants
ntrials = height(simulated_data.sim_data_01.sim_01); % number of trials

% names of the free parameters in the order of params(1), params(2), ...
param_names = {'alpha_90_presented_pos', 'alpha_70_presented_pos', 'alpha_50_presented_pos', ...
    'alpha_90_omitted_pos', 'alpha_70_omitted_pos', 'alpha_50_omitted_pos', ...
    'alpha_90_presented_neg', 'alpha_70_presented_neg', 'alpha_50_presented_neg', ...
    'alpha_90_omitted_neg', 'alpha_70_omitted_neg', 'alpha_50_omitted_neg', ...
    'beta', 'lambda'};

% lower and upper bound for fit
LB = [0 0 0 0 0 0 0 0 0 0 0 0 0 1]; % lower bound
UB = [1 1 1 1 1 1 1 1 1 1 1 1 100 1000]; % upper bound

for i = 1:nsubs % loop through participants

    i
    data = simulated_data.(sprintf('sim_data_%02d', i));

    for k = 1:nsim % loop through simulated data sets of this participant

        data_nsim = data.(sprintf('sim_%02d', k));
        sub_choice = data_nsim.chosen;
        sub_not_choice = data_nsim.unchosen;
        sub_outcome = data_nsim.feedback;
        sub_stimuli = data_nsim.stim;
        sub_feedback_type = data_nsim.feedback_type;

        n_valid = sum(~isnan(sub_choice) & ~isnan(sub_outcome));

        all_params = NaN(niter, numel(param_names));
        all_ll = NaN(1, niter);
        all_bic = NaN(1, niter);

        for j = 1:niter % iterations of fmincon

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

            % half-life in presentations, log-uniform in [1, 1000]
            lambda = exp(rand*log(1000));

            params = [alpha_90_presented_pos,alpha_70_presented_pos,...
                alpha_50_presented_pos,alpha_90_omitted_pos,alpha_70_omitted_pos,...
                alpha_50_omitted_pos,alpha_90_presented_neg,alpha_70_presented_neg,...
                alpha_50_presented_neg,alpha_90_omitted_neg,alpha_70_omitted_neg,...
                alpha_50_omitted_neg,beta, lambda];

            options = optimset('display','off','MaxFunEvals',100000,'TolFun',1e-16,'TolX',1e-16);
            [params, LL] = fmincon(@probfb_fun_model10ab,params,[],[],[],[],LB,UB,[],options,......
                sub_choice,sub_outcome,sub_not_choice,ntrials, sub_stimuli, sub_feedback_type);

            all_params(j,:) = params;
            all_ll(j) = LL;
            all_bic(j) = 2*LL + length(params)*log(n_valid);
        end

        % save best fit according to -LL and according to BIC
        [fit.ll(k,i), best] = min(all_ll);
        [fit_BIC.bic(k,i), bestBIC] = min(all_bic);
        for p = 1:numel(param_names)
            fit.(param_names{p})(k,i) = all_params(best,p);
            fit_BIC.(param_names{p})(k,i) = all_params(bestBIC,p);
        end
    end
end
