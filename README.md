# Drosophila.multiallelic
Scripts from project "Selection of multiallelic loci in *D. melanogaster*"

## Project description
In this repository I collected the scripts used for my internship project about selection at multiallelic loci in Drosophila melanogaster. My project is based on the model developed by Nikolas Vellnow, Toni Goosmann and David Waxman. Their preprint can be found on bioRxiv: https://doi.org/10.1101/2024.11.08.622587
Many of the scripts I used for my project are originally from Vellnow et al. I kept some of them as they are, changed some of them slightly, some of them more and created new scripts as well. You can find the orignal scripts used by Vellnow et al. in his Github repository: https://github.com/NikolasVellnow/selection_multiallelic
We have used Matlab and Python to apply the model to empirical data from *D. melanogaster*.

### Extracting trajectories from pooled sequencing data
We downloaded `.bam`-files for different sampling time points from the DEST webpage https://dest.bio/data-files/sync-bam-bed. We had to modify the bam files so that the read groups had unique identifiers, which we did with the shell script `fix_read_groups.sh`.
The shell script `run_extraction.sh` takes these modified `.bam`-files as input and uses freebayes to call alleles and outputs these to a `.vcf`-file. It then filters the sites so that only sites with at least 3 alleles are included in the output `filtered_sites.tsv`.
A textfile with a list of the bam files and the `.fasta`-file of the *D. melanogaster* reference genome is necessary for this script.
The file `run_extraction_coveragefilter.sh` does the same, but it uses customized filter settings based on the coverage of your sample instead of the FreeBayes default filter settings.
The Python script `data_conversion.py` takes the textfile `filtered_sites.tsv` as input to calculate the number of sites and the number of alleles at those sites and save them in the output `pop_a_b_allelefreq.tsv`, where pop stands for your population, a for the threshold for hard filtering with the `min_reads` constant and b for the threshold for filtering with max_n_alleles. Those information, as well as the estimated generations for your population samples can be changed in the first lines of the script in the section "your settings".

### Genes
The script `makebed.py` uses the file `pop_a_b_allelefreq.tsv` to generate a `.bed`-file `pop.bed`. The script `genefinder.py`can then take the `.bed`-file `pop.bed` together with a `.fasta`-file of the *D. melanogaster* reference genome and the corresponding gene annotation `.gff`-file to save the names of the genes that our loci are in in the `.txt`-file `genes_pop.txt`.
The script `count_words.py` takes the `.txt`-file `genes_pop.txt` to save a `.tsv`-file `genes_pop.tsv`, where you have every gene name only once in the first column and its frequency in `genes_pop.txt` in the second column. It also outputs the proprtion of your loci that was found in genes and needs to load the file `pop_a_b_allelefreq.tsv` to do so.
With `proportions.py` you can calculate the proportion of DNA in genes in the chromosomes you are investigating as a reference. The script loads the `.fasta`-file of the *D. melanogaster* reference genome and the corresponding gene annotation `.gff`-file for that.

### Estimation of F
The script `scan_dros_genome_F.m` takes a `.tsv`-file with trajectory data across many genomic windows, e.g. `pop_a_b_allelefreq.tsv`, estimates F for each window and saves the information in a `.mat`-file `F_results.mat`. The script calls the functions `fit_F_sparse_func.m`, with its dependencies. Also, the "Optimization Toolbox" and the "Global Optimization Toolbox" for Matlab have to be installed to use this script.
The script `scan_dros_genome_analysis.m` then takes the estimated F matrices `F_results.mat` analyses them and outputs the results in a `.mat`-file `F_analysis_results.mat`. The script calls function `convert_F_to_w.m`

#### Manhattan Plots
The script `manhattan_plot.m` takes `F_analysis_results.mat` to create manhattan plots for the strength of heterozygous advantage based on the value hav1. It also displays the quantiles of hav1.

#### Trajectory plots for subgroups
The script `subgroups.m` also takes `F_analysis_results.mat` as input and displays the quantiles of hav1. Additionally, the script defines subgroups of the loci based on the analysis done before and saves them in `.csv`-files. `sample_for_plots.py`now takes such a subgroup `.csv`-file as input to pick randomly some of the included loci and save a `.csv`-file `pop_chr_pos.csv` for that, where pop stands again for the population and can be changed and chr and pos describe the Chromosome and Position of the specific locus. The `sample_size` can also be changed. The script generates a `.txt`-file `filelist.txt`, which stores the names of all the `.csv`-files created by the same script.
You can then use the script `trajectory_plots.m`, which takes `filelist.txt`, the corresponding `.csv`-files `pop_chr_pos.csv` and the estimated F matrices `F_results.mat` to create trajectory plots based on the observed data from the `.csv`-files, as well as modeled trajectory plots based on the estimated F matrices in `F_results.mat`.





