clear all
close all
clc

% load F matrixes
F_analysis_results = load("F_analysis_results.mat");

data_struct = F_analysis_results.F_analysis_results;

data_table = struct2table(data_struct);



% Convert the 'type_top_genotype' column to categorical
chroms = categorical(data_table.chr);

% Extract the 'pos' column (numerical data)
positions = data_table.pos;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exract hav1
hav1_values = data_table.hav1;

% keep only valid values (nan and inf come if you try dividing by zero)
valid_indices = ~isnan(hav1_values) & ~isinf(hav1_values);

% filter positions, chromosomes and hav1
positions_cleaned = positions(valid_indices);
chroms_cleaned = chroms(valid_indices);
hav1_cleaned = hav1_values(valid_indices);


%%%%%%%%%%%%%%%%%%%%
%quantiles
% Calculate quantiles for hav1_log
quantiles_hav1 = quantile(hav1_cleaned, [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]);

disp(quantiles_hav1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manhattan plot for hav1

% generate a color for each chromosomes
unique_chromosomes = categories(chroms_cleaned); 
colors = lines(length(unique_chromosomes)); 

% plot
figure;
hold on;

num_chromosomes = length(unique_chromosomes);
chromosome_width = 1; % width of each chromosome

% loop through all chromosomes
for i = 1:num_chromosomes
    current_chromosome = unique_chromosomes(i); %current chromosome
    
    % find indexes for current chromosome
    idx = chroms_cleaned == current_chromosome;
    
    % calculate x-position with offset for current chromosome
    % dots are spread within [i - category_width/2, i + category_width/2]
    adjusted_positions = positions_cleaned(idx) * chromosome_width / max(positions_cleaned) + (i - 0.5) * chromosome_width;
    
    % generate scatter-plot for current chromosome
    scatter(adjusted_positions, hav1_cleaned(idx), 5 , colors(i,:), 'filled');
end

% add horizontal line at y=0
yline(0, 'r--', [], 'LineWidth', 1);

% add horizontal lines for each quantile
for i = 1:length(quantiles_hav1)
    yline(quantiles_hav1(i), 'g--', [], 'LineWidth', 1);
end

% add labels and title
xlabel('Chromosome');
ylabel('hav1');

% set x ticks for each chromosome
xticks(1:num_chromosomes); 
xticklabels(unique_chromosomes);

hold off;

% save plot
saveas(gcf, "hav1_manhattan.fig");
print("hav1_manhattan.png", '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% manhattan plot for hav1 only within first and last quantile
% (without outliers for better visualization)

% keep only valid values (nan and inf come if you try dividing by zero)
valid_indices = ~isnan(hav1_values) & ~isinf(hav1_values) & (hav1_values > quantiles_hav1(1)) & (hav1_values < quantiles_hav1(length(quantiles_hav1)));

% filter positions, chromosomes and hav1
positions_cleaned = positions(valid_indices);
chroms_cleaned = chroms(valid_indices);
hav1_cleaned = hav1_values(valid_indices);


% generate a color for each chromosomes
unique_chromosomes = categories(chroms_cleaned); 
colors = lines(length(unique_chromosomes)); 

% plot
figure;
hold on;

num_chromosomes = length(unique_chromosomes);
chromosome_width = 1; % width of each chromosome

% loop through all chromosomes
for i = 1:num_chromosomes
    current_chromosome = unique_chromosomes(i); %current chromosome
    
    % find indexes for current chromosome
    idx = chroms_cleaned == current_chromosome;
    
    % calculate x-position with offset for current chromosome
    % dots are spread within [i - category_width/2, i + category_width/2]
    adjusted_positions = positions_cleaned(idx) * chromosome_width / max(positions_cleaned) + (i - 0.5) * chromosome_width;
    
    % generate scatter-plot for current chromosome
    scatter(adjusted_positions, hav1_cleaned(idx), 5, colors(i,:), 'filled');
end

% add horizontal line at y=0
yline(0, 'r--', [], 'LineWidth', 1);

% add horizontal lines for each quantile
for i = 1:length(quantiles_hav1)
    yline(quantiles_hav1(i), 'g--', [], 'LineWidth', 1);
end


% add labels and title
xlabel('Chromosome');
ylabel('hav1');

% set x ticks for each chromosome
xticks(1:num_chromosomes); 
xticklabels(unique_chromosomes);

hold off;

% save plot
saveas(gcf, "hav1_woexc_manhattan.fig");
print("hav1_woexc_manhattan.png", '-dpng', '-r300');
