function [Pe, SNRdB, CM, num_total_cells] = run_main_GM(cell_label, donor_label, fname)
% Function runs the communication model for a single cell subtype and a single donor and returns
% the average probability of error, SNR in dB, and confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
% cell_label:           label used to identify the cell subtype
% donor_label:          label to identify the donor (BC or HD)
% fname:                output file name. If set to zero, automatically generates the file name

%% Outputs:
% Pe:                   average probability of error
% SNRdB:                estimated average signal-to-noise ratio in dB
% CM:                   confusion matrix
% num_total_cells:      number of total cells over which average Pe is estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define cytokines
cytokines =  {'IFNg', 'IL6', 'IL10', 'IL2', 'IL4', 'Untreated'}; 
Save_ON = true; %false;

% Read donor data
if contains(donor_label, 'BC')
    data_fname = ['Data/BC Samples/' donor_label '_' cell_label '.mat'];
elseif contains(donor_label, 'HD')
    data_fname = ['Data/HD Samples/' donor_label '_' cell_label '.mat'];
end
response = load(data_fname);

% Load alpha parameters
if contains(donor_label, 'BC')
    if contains(cell_label, 'B')
        param_fname = 'AlphaVals/BC/BC_CD20+ B_alpha.mat';
    elseif contains(cell_label, 'CD4+')
        param_fname = 'AlphaVals/BC/BC_CD4+ T_alpha.mat';
    elseif contains(cell_label, 'CD8+')
        param_fname = 'AlphaVals/BC/BC_CD8+ T_alpha.mat';
    else
        param_fname = ['AlphaVals/BC/BC_' cell_label '_alpha.mat'];
    end
elseif contains(donor_label, 'HD')
     if contains(cell_label, 'B')
        param_fname = 'AlphaVals/HD/HD_CD20+ B_alpha.mat';
    elseif contains(cell_label, 'CD4+')
        param_fname = 'AlphaVals/HD/HD_CD4+ T_alpha.mat';
    elseif contains(cell_label, 'CD8+')
        param_fname = 'AlphaVals/HD/HD_CD8+ T_alpha.mat';
    else
        param_fname = ['AlphaVals/HD/HD_' cell_label '_alpha.mat'];
     end
end

avals = load(param_fname);

% Set output filename
out_fname = set_output_filename(cell_label, donor_label, fname);


% Run model
[out_data, GM, SNRdB] = preprocess_data(response, avals.alpha, cytokines);
[Pe, CM, num_total_cells] = estimate_error_rate(out_data, GM);

if Save_ON
    save(out_fname);
end

% Display outputs
disp('------------------ Results -----------------'); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename = set_output_filename(cell_label, donor_label, fname)

a = clock; a = a(2:5); % get time as year|month|day|hour|minute
d_str = num2str(a);
timestamp = d_str(~isspace(d_str)); % create timestamp

% Directory where results are saved
if ~exist('Results','dir')
    mkdir('Results');
end
ResultPath = ['Results' filesep];

% Determine whether to automatically generate a filename based on
% parameters or use a user-specified filename
switch fname
    case 0 % auto generate filename
        filename = auto_gen_filename(cell_label, donor_label, timestamp);
        filename = [ResultPath filename];
    otherwise %user-specified filename
        filename = [ResultPath fname];
end

end

function filename = auto_gen_filename(cell_label, donor_label, timestamp)

filename_cell{1} = 'RunSignalingModel';
filename_cell{2} = ['_' cell_label];
filename_cell{3} = ['_' donor_label];
filename_cell{4} = ['_' timestamp];
% Concatenate all entries to form filename
filename = strcat(filename_cell{:});
end
