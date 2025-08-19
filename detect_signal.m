function det_signal_idx = detect_signal(measurement, GM)
% Function performs Maximum Likelihood (ML) detection under a Gaussian mixture model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
% measurement:      vector of single cell measurements of length equal to the number of measured pSTATs
% GM:               cell array of Gaussiam-Mixture models corresponding to each possible input signal/cytokine 
%                   in the communication model
%% Outputs:
% det_signal_idx:   index correponding to the input signal that is detected under the ML criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_dims = numel(measurement);
num_signals = numel(GM);
num_componets = numel(GM{1});

prob = zeros(num_componets,num_signals);
log_prob = zeros(num_signals,1);
        
% Perform ML detection
% Compute the probabilities for each possible input signal under the GM
% model
% Select the signal that maximizes the probability

for signal_idx = 1:num_signals
    gmd = GM{signal_idx};
    rho = gmd.ComponentProportion;
    
    for comp_idx = 1:num_componets
        mu = gmd.mu(comp_idx,:); 
        Sigma = gmd.Sigma(:,:,comp_idx); 
        
        S = det(Sigma);

        prob(comp_idx,signal_idx) = rho(comp_idx)/(2*pi)^(num_dims/2)/sqrt(S) * ...
                     exp(-0.5*(measurement - mu)/(Sigma)*(measurement - mu)');
    end
    
    log_prob(signal_idx) = log(sum(prob(:,signal_idx)));
end

% Find index corresponding to the maximum log probability

[~,det_signal_idx] = max(log_prob);

end
