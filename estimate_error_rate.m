function [Pe, CM, total_cells] = estimate_error_rate(all_data, GM)
% Function reads and processes FACS data, and outputs signal (cytokine) detection error rates 
% and confusion matrix (CM) using Gaussian-Mixture model   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
% all_data:         cell array of all cell measurements (for a cell subtype) for each input signal 
%                   signal/cytokine. The length of the cell array is equal to the number of cytokines. 
%                   Each entry in the cell array is a matrix with number of columns equal to the number 
%                   of measured pSTAT/pSMADs and the number of rows equal to the number of cells.  
% GM:               cell array of Gaussiam-Mixture models corresponding to each input signal/cytokine 
%                   in the communication model        
%% Outputs:
% Pe:               average probability of error
% CM:               confusion matrix
% total_cells:      number of total cells over which average Pe is estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_signals = numel(all_data);

total_cells = 0;
CM = zeros(num_signals);

for signal_idx = 1:num_signals

    signal_data = all_data{signal_idx};
    num_cells = size(signal_data,1);
    total_cells = total_cells + num_cells;

    for cell_idx = 1:num_cells % loop over trials
        response = signal_data(cell_idx,:);

        det_signal_idx = detect_signal(response, GM);

        CM(signal_idx, det_signal_idx) = CM(signal_idx, det_signal_idx) + 1;

    end
end

% average probability of error as one minus average probability of correct
% signal detection
Pe = 1 - sum(diag(CM))/total_cells; 

% confusion matrix
CM = CM/total_cells;

end

