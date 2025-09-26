% Loops over (part of) drosophila genome and analyzes F for each window

clear all
close all
clc

% load F matrixes
F_results_unfiltered = load("F_results_allc_1.mat").F_results;

% filter:
% Exclude windows where optimization failed
keep_indices = ~isnan([F_results_unfiltered.cost]);
F_results = F_results_unfiltered(keep_indices);
opt_failed = isnan([F_results_unfiltered.cost]);
n_opt_failed = length(F_results_unfiltered(opt_failed));
p_opt_failed = 100 * (n_opt_failed / length(F_results_unfiltered));

% Exclude windows where freqs where all 0.25, F only zeros and cost = 0
F_results = F_results([F_results.cost] ~= 0);

% extract and save number of windows
n_wins = length(F_results);
disp(n_wins);

% initialise with 0

sum_cost = 0;
n_three = 0;

%loop through all windows
for i = 1:n_wins
    F = F_results(i).F; % extract F
    if size(F, 2) == 3 % only if it is a locus with 3 alleles
        n_three = n_three + 1; % Count the number of loci with 3 alleles
        sum_cost = sum_cost + F_results(i).cost; % Accumulate the cost for each window
    end
end

average_cost_three = sum_cost/n_three;
disp(n_three);
disp(average_cost_three);
  