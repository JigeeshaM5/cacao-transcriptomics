#ORA function (example usage in cacao leaf)

ORA_leaf_geno2 <- lapply(seq_along(leaf_sets_geno), function(i) {
  fora(pathways = sl, 
       genes = leaf_sets_geno[[i]][, geneid], # genes in go_term i
       universe = leaf_reference, # all genes expressed in leaves
       minSize = 15, 
       maxSize = 500) %>% 
    mutate(cluster = names(leaf_sets_geno)[[i]]) # add intersection/complement names column
}) %>% 
  data.table::rbindlist() %>% # combine tables
  filter(padj < 0.05) %>% 
  arrange(cluster, padj) %>% 
  # Add additional columns from BP_db
  left_join(distinct(go_terms, ids, go_names, category),
            by = c("pathway" = "ids")) 

