% Loops over (part of) drosophila genome and estimates F for each locus

clear all
close all
clc


% load trajectory data that was extracted with freebayes and the
% data_conversion script
dros_data = readtable("linvilla_1_10_allelefreq.tsv",'FileType', 'text', 'Delimiter', '\t');

chr_column = string(dros_data.chr); % Convert cell array to string array
pos_column = string(dros_data.pos); % Convert numeric values to strings
unique_combinations = unique([chr_column, pos_column], 'rows');
n_tot_wins = size(unique_combinations, 1);


chrom_vector = unique(dros_data.chr)';
n_chroms = length(chrom_vector);

% Initialize empty struct array to save results in
% Preallocate struct array for results
temp_struct = struct('chr', [], ...
    'pos', [], ...
    'F', [], ...
    'cost', [], ...
    'F_fminunc', [], ...
    'cost_fminunc', [], ...
    'comment', []);
F_results = repmat(temp_struct, n_tot_wins, 1);

current_win = 0;

% loop over all chromosomes
for j = 1:n_chroms

    chr = string(chrom_vector{j});


    % Create a sub-table for current chromosome
    chrom_data = dros_data(dros_data.chr == chr, :);
    
    win_vector = unique(chrom_data.pos)';
    n_wins = length(win_vector);
    
    % loop over all windows
    for i=1:n_wins
        current_win = current_win +1;
        position = win_vector(i);
        disp(chr + ": position " + position);
        win_data = chrom_data(chrom_data.pos == position, :);
    
        n = unique(chrom_data{chrom_data.pos == position, "n_alleles"});
        % generate names for allele frequencies depending on the number of
        % alleles
        af_columns = cell(1, n); % cell to save names of allele frequencies
        for k = 1:n
            af_columns{k} = sprintf('af_%d', k); % generate 'af_1', 'af_2', ..., 'af_n'
        end
        
        trajectory_data = win_data{:, af_columns}';
    
        generation_data = win_data{:, 'gen'}';
    
        try
            [F, Cost_opt, F_fminunc, Cost_opt_fminunc] = fit_F_sparse_func(trajectory_data, generation_data, n);
            
            F_results(current_win).chr = chr;
            F_results(current_win).pos = position;
            F_results(current_win).F = F;
            F_results(current_win).cost = Cost_opt;
            F_results(current_win).F_fminunc = F_fminunc;
            F_results(current_win).cost_fminunc = Cost_opt_fminunc;
        
        catch ME
            F_results(current_win).chr = chr;
            F_results(current_win).pos = position;
            F_results(current_win).F = NaN;
            F_results(current_win).cost = NaN;
            F_results(current_win).F_fminunc = NaN;
            F_results(current_win).cost_fminunc = NaN;
            F_results(current_win).comment = ME.message;
            disp("Error in " + chr + ", position " + position ...
                + ": " + ME.message);
        end
    end
end 

% filename
file_str = "F_results" + ".mat";

% save matrix as .mat file
save(file_str, 'F_results');