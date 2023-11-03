For step 1.1 

plink --vcf gwas_data/genotypes.vcf --make-bed --pca 10 'header' --out genotypes

For step 2.1

plink --vcf gwas_data/genotypes.vcf --freq

For step 3.1

plink --vcf gwas_data/genotypes.vcf --linear --pheno gwas_data/CB1908_IC50.txt --covar genotypes.eigenvec --allow-no-sex --out CB1908_IC50_gwas_results

plink --vcf gwas_data/genotypes.vcf --linear --pheno gwas_data/GS451_IC50.txt --covar genotypes.eigenvec --allow-no-sex --out GS451_IC50_gwas_results

Step 3.4

For CB1908 the top loci was within the TUBA1A gene which encodes tubulin alpha 1a. This could have made individuals more sensitive to CB1908 disrupting microtubule function.

For GS451 the top loci was within the pseudogene ZNF826P on chromosome 19. This pseudogene has been previously identified as being a predictior for lung cancer. The mutation here could modify regulatory functions of any produced transcripts, the change in regulatory function could increase the sensitivity to GS451 within immune cells. 