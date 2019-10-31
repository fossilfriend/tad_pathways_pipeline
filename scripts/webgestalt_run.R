# 2017 Gregory Way
# scripts/webgestalt_run.R

# Description: 
# Runs WebGestalt Pathway Analysis on Input genelist

# Usage:
# The script is run in the command line

#           Rscript --vanilla scripts/webgestalt_run.R

# and takes two positional arguments as input:

#         --tad_genelist_file       <LOCATION_OF_GWAS_FILE.csv>
#         --output_name             <FILENAME_ABBREVIATION>

# Output:
# Significantly overrepresented pathways from a WebGestalt Analysis

library(WebGestaltR)
library(tidyr)
library(dplyr)

# Load in command arguments
option_list <- list(optparse::make_option(c("-t", "--tad_genelist_file"),
                                          type = "character",
                                          help = "Location of TAD genes"),
                    optparse::make_option(c("-o", "--output_name"),
                                          type = "character",
                                          help = "Abbrev. to save results"),
                    optparse::make_option(c("-e", "--enrich_database"),
                                          type = "character",
                                          help = "Abbrev. to save results"))

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

# Load arguments
tad_gene_file <- opt$tad_genelist_file
output_name <- opt$output_name
enrich_db  <- c("geneontology_Biological_Process",  "geneontology_Molecular_Function",  "geneontology_Cellular_Component")
if (opt$enrich_database == "pathway") {
    enrich_db  <- c("pathway_KEGG", "pathway_Reactome", "pathway_Wikipathway")
}
if (opt$enrich_database == "tfbs") {
    enrich_db  <- c("network_Transcription_Factor_target")
}
print(paste("EnrichDB", enrich_db))

output_pval_file <- file.path("gestalt", paste0(output_name, "_pvals.tsv"))
output_path_file <- file.path("gestalt", paste0(output_name, "_gestalt.tsv"))

gene_df <- readr::read_tsv(tad_gene_file)
genes <- gene_df$gene_name

webgestalt_output <- WebGestaltR(enrichMethod = "ORA",

                                 enrichDatabase = enrich_db,
                                 organism = "hsapiens",
                                 interestGene = genes,
                                 interestGeneType = "genesymbol",
                                 minNum = 4,
                                 fdrMethod = "BH",
                                 isOutput = TRUE,
                                 outputDirectory = "gestalt",
                                 referenceSet = "genome",
                                 projectName = output_name)

# print(str(webgestalt_output))
# Process output files
webgestalt_output <- webgestalt_output %>%
  tidyr::separate_rows(overlapId, userId, sep = ";")
p_val <- webgestalt_output %>% dplyr::select(description, pValue) 
p_val <- p_val[!duplicated(p_val), ]

colnames(webgestalt_output) <- c("go_id", "go_name", "link", "count",
                                 "observed", "expected", "R", "pval",
                                 "adjP", "overlapGene", "enrichdb", "symbol")
colnames(p_val) <- c("go_name", "adjP")

write.table(p_val, output_pval_file, sep = "\t", row.names = FALSE)
write.table(webgestalt_output, output_path_file, sep = "\t", row.names = FALSE)
