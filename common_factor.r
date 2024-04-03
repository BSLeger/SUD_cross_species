library(GenomicSEM)
library(here)

ldsc_filepath <- here("../nicsa_LD_score_regression/initial_run/LDSCoutput.RData")
sumstats_filepath <- here("../SUMSTATSoutput.RData")

load(ldsc_filepath)
load(sumstats_filepath)

#run the multivariate GWAS using parallel processing
pfactor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = SUMSTATSoutput, 
                            estimation = "DWLS", cores = NULL, parallel = TRUE)

#output_path <- paste(results_dir, "/results_subset_", subset_number, ".csv", sep = "")
output_path='pfactor_out_allfiles.csv'
write.csv(pfactor, file = output_path, row.names = FALSE)
#print(paste("Analysis for subset", subset_number, "completed and saved to", output_path))
