function [out_data, GM, SNRdB] = preprocess_data(in_data, avals, cytokines)
% Function reads and pre-processes FACS data using the log transformation. It outputs the processed 
% data, the Gaussian-Mixture models and the average SNR in dB   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
% in_data:          cell array of all cell measurements (for a cell subtype) for each input signal 
%                   signal/cytokine. The length of the cell array is equal to the number of cytokines. 
%                   Each entry in the cell array is a matrix with number of columns equal to the number 
%                   of measured pSTAT/pSMADs and the number of rows equal to the number of cells  
% avals:            matrix of alpha values used in the pre-processing step.The number of columns is equal 
%                   to the number of measured pSTAT/pSMADs and the number of rows is equal to the number of
%                   available cytokines
% cytokines:        cell array of available cytokines      
%% Outputs:
% out_data:         cell arary of processed data
% GM:               cell array of Gaussian Mixture models. Each entry in the cell araryis a GM model with two 
%                   components for a single signal/cytokine
% SNRdB:            estimated average signal-to-noise ratio in dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_signals = numel(cytokines);
GM = cell(1, num_signals);
out_data = cell(1, num_signals);

options = statset('MaxIter',1000);
num_components = 2; 
Ps = zeros(num_components,1);
Pn = zeros(num_components,1);
avgPs = zeros(num_signals,1);
avgPn = zeros(num_signals,1);

for signal_idx = 1:num_signals
    signal = cytokines{signal_idx};
    switch signal
        case 'IL2'
            data = in_data.response.IL2;
        case 'IL4'
            data = in_data.response.IL4;
        case 'IL10'
            data = in_data.response.IL10;
        case 'IFNg'
            data = in_data.response.IFNg;
        case 'IL6'
            data = in_data.response.IL6;
        case 'Untreated'
            data = in_data.response.Untreated;
    end

    out_data{signal_idx} = log(data - avals(signal_idx,:) + 1);

    Sigma = cov(out_data{signal_idx});
    RegularizationValue = min(eig(Sigma));

    GM{signal_idx} = fitgmdist(out_data{signal_idx},num_components,'CovarianceType','full',...
                'SharedCovariance',false,...
                'RegularizationValue',RegularizationValue,...
                'Options',options);

    for comp_idx = 1:num_components       
        mu = GM{signal_idx}.mu(comp_idx,:); 
        Sigma = GM{signal_idx}.Sigma(:,:,comp_idx);
        rho = GM{signal_idx}.ComponentProportion(comp_idx);        

        Ps(comp_idx)=rho*mu*mu';
        Pn(comp_idx)=rho*trace(Sigma);
    end
    avgPs(signal_idx) = sum(Ps);
    avgPn(signal_idx) = sum(Pn);

end

SNRdB = 10*log10(sum(avgPs)/sum(avgPn));

end