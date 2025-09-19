echo "Running freebayes"
freebayes -f /vol/storage/alenas_stuff/bamfiles/dmel-all-chromosome-r6.12.fasta\
 --pooled-continuous\
 -r 2R -r 2L -r 3R -r 3L \
 --min-coverage 611\
 -g 3663\
 -L /vol/storage/alenas_stuff/bamfiles/linvilla/list.txt > variants_allc_pooled_test.vcf
  



echo "Compressing VCF-file"
bgzip variants_allc_pooled_test.vcf
echo "Creating index for compressed VCF-file"
tabix variants_allc_pooled_test.vcf.gz
echo "Extracting allele counts with bcftools"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' variants_allc_pooled_test.vcf.gz > allele_counts_snps_allc_test.tsv

echo "Filtering with awk"
awk 'BEGIN { FS="\t"; OFS="\t" }
{
    # Split ALT field (column 4) on commas.
    n = split($4, altAlleles, ",");
    # Only consider sites with at least 2 alternate alleles (i.e. REF + >=2 ALTs = 3 alleles)
    if(n >= 2) {
        valid = 1;
        # Check that all sample fields (columns 5 to NF) are not missing (i.e. not ".")
        for(i = 5; i <= NF; i++){
            if($i == ".") { valid = 0; break; }
        }
        if(valid) print $0;
    }
}' allele_counts_snps_allc_test.tsv > filtered_sites_allc_test.tsv

echo "Done!"

