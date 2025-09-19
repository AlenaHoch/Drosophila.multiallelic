% load here the files that you get from the python script sample_for_plots
% and the F-Matrix you got from scan_dros_genome_F.m

clear all
close all
clc

F_results = load("F_results_allc_1.mat").F_results;


% read filenames from txt file
txt_file = "filelist.txt"; % name of txt file
file_names = readlines(txt_file);

% loop through all filenames of csv files
for i = 1:length(file_names)
    file_str = file_names(i); % current filename of csv file
    % check if csv file exists
    if isfile(file_str)
        
        disp(['file: ', file_str]);
       
        % load csv file as table
        trajectory_data = readtable(file_str);
        
        % save and display number of alleles at current position (in this
        % csv file)
        n = trajectory_data{1, "n_alleles"};
        disp('number of alleles at this position:');
        disp(n);
        
        % generate names for allele frequencies depending on the number of
        % alleles
        af_columns = cell(1, n); % cell to save allele frequency names
        for i = 1:n
            af_columns{i} = sprintf('af_%d', i);
        end
        
        % ! save in matrix
         X = trajectory_data{:, af_columns}'; %read in the allele frequencies
        
        
        % extract generation data
        g = trajectory_data{:, "gen"}';
        
     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulating trajectories for estimated F
        % Parameters
                                                
        T=max(g);                                      % Final time
         
        start_freqs = X(:, 1);               % Starting freqs from first drosophila sample
        
        % Initialisation
        X_sim=zeros(n,T+1);                           % Stores deterministic trajectories
        X_sim(:,1) = start_freqs;              % Assign starting freqs
        
        chr = trajectory_data{:, "chr"}';
        pos = trajectory_data{:, "pos"}';
    
        % find matching F matrix
        matching_index = find(([F_results.chr] == trajectory_data{1, "chr"}) & ([F_results.pos] == trajectory_data{1, "pos"}));
    
        if ~isempty(matching_index)
            F = F_results(matching_index).F;
        end
    
        % Time loop
        for k=1:T
            x=X_sim(:,k);
            V=diag(x)-x*x';
            D=V*F*x/(1+x'*F*x);
            xp=x+D;
            xp=xp/sum(xp);              % Forces normalisation (just to be sure!)
            X_sim(:,k+1)=xp;
        end
        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plots
        
        % generate matrix of n random colors (RGB values)
        colors = rand(n, 3);
        
        % cell to save legend entries
        legend_entries = cell(1, n);
        
        
        g_sim = 0:1:T;
       
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Combined plot
        
        font_size = 18;
        font_size_ticks = 16;
        
        figure;
        % ! Empirical trajectories plot
        subplot(1,2,1);
        hold on;
        % loop through all alleles and plot their frequency
        for i = 1:n
            plot(g, X(i,:), 'color', colors(i,:), 'LineWidth', 2);
            
            % add allele to legend
            legend_entries{i} = sprintf('af_%d', i);
        end
        
        % position legend
        legend(legend_entries(1:n), 'Location', 'northwest');
        
        title('Drosophila observed evolution', 'FontSize', font_size);
        xlabel('Generations', 'FontSize', font_size);
        ylabel('Allele frequency', 'FontSize', font_size);
        ax = gca;
        ax.FontSize = font_size_ticks;
        ylim([0.0,1.0]);
        xlim([0, T]);
        
        %%%%
        % ! Modelled trajectories plot
        subplot(1,2,2);
        hold on;
      
        %  loop through all alleles and plot their modeled frequency
        for i = 1:n
            plot(g_sim, X_sim(i,:), 'color', colors(i,:), 'LineWidth', 2);
            
            % add allele to legend
            legend_entries{i} = sprintf('af_%d', i);
        end
        
        % position legend
        legend(legend_entries(1:n), 'Location', 'northwest');
        
        title('Model trajectories', 'FontSize', font_size);
        xlabel('Generations', 'FontSize', font_size);
        ylabel('Allele frequency', 'FontSize', font_size);
        ax = gca;
        ax.FontSize = font_size_ticks;
        ylim([0.0,1.0]);
        xlim([0, T]);
        
        hold off;
        
        % Adjust figure size and appearance
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 5]); % Set size to width=6in &; height=5in
        print("trajectories_emp_model_combined_" + file_str + ".png",'-dpng','-r300');
           
        
    else
        disp(['file not found: ', file_str]);
    end
    
end



