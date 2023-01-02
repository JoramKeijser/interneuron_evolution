# Utilities/helper functions for preprocessing

# Gene names in "Elfn1" format (mouse convention)
recode_genes <- function(gene_names){
  sapply(gene_names, function(name){
    paste0(substr(name, 1, 1), tolower(substr(name, 2, nchar(name))))
  })
}
