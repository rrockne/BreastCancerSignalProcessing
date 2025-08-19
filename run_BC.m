% Script that sets up the simulation for all breast cancer (BC) donnors and all cell subtypes. Results are 
% saved in a summary output file as well as individual files for each donnor and cell subtype 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unzip('Data.zip')
addpath("AlphaVals/")
addpath("AlphaVals/BC")
addpath("AlphaVals/HD")
addpath("Data")
addpath("Data/BC Samples/")
addpath("Data/HD Samples/")

cell_types = { 'CD4+ TCM', 'CD4+ TEM', 'CD4+ Naive', 'CD4+ T', 'CD8+ TCM', 'CD8+ TEM', 'CD8+ Naive', 'CD8+ T', ...
               'Naive B', 'Memory B', 'CD20+ B', 'Classical monocytes', 'CD16+ NK'};
outfile = 'Results_BC_All.mat'; 

num_cell_types = numel(cell_types);
num_samples = 19; % number of BC donors

snr_table = zeros(num_samples,num_cell_types);
Perror_table = zeros(num_samples,num_cell_types);
num_cells = zeros(1,num_samples);
CM_all{num_samples}{num_cell_types} = [];
tic
for ii = 1:num_cell_types
    filename = ['casefile_BC_' cell_types{ii} '.csv'];
    Q = zeros(6);
    for jj = 1:num_samples
        donor_label = ['BC' num2str(jj)];
	    cell_label = cell_types{ii};
	    [Pe, SNRdB, CM, num_cells(jj)] = run_main_GM(cell_label, donor_label, 0);
        
        disp(donor_label)
        disp(SNRdB)
        disp(Pe)

        Q = Q + CM*num_cells(jj);
        snr_table(jj,ii) = SNRdB;
        Perror_table(jj,ii) = Pe;
        CM_all{jj}{ii} = CM;        
    end
    Q = Q/sum(num_cells); %average confusion matrix
    toc
end

save(outfile, 'Perror_table', 'snr_table', 'CM_all', 'Q');
