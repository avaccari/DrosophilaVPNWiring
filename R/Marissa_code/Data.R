### Read in the output matrices for each LC. 
# Contains the individual LCs (columns) and outputs with >10 synapses total (rows). 
print(getwd())
lc17out = read.csv("../../output/LC_outputs_v12/LC17outliersmeta.csv")
outrm = as.character(lc17out$bodyid)

loom_list = c("LC6", "LC4", "LPLC2", "LPLC1", "LC12", "LC17")
nonloom_list = c("LC21", "LC25", "LC15", "LC11", "LC18")

mark_lcs = c("LC17", "LC12", "LC4", "LC6", "LPLC1", "LPLC2", 
             "LC9", "LC10", "LC11", "LC13", "LC15", "LC16", 
             "LC18", "LC20", "LC21", "LC22", "LC24", "LC25", 
             "LC26",  "LPLC4", "LC28b")

LCs_interest = mark_lcs#c(loom_list, nonloom_list)
mat_list = list()

# Iterate over every VPN
for (lc in LCs_interest){#c(loom_list, nonloom_list)){
  # Read connectivity data for that VPN
  mat_list[[lc]] = readRDS(file=sprintf("../../output/LC_outputs_v12/%s_LC_outputs_mat.rds", lc))
  tmp = colnames(mat_list[[lc]][-c(1:3)]) %>% str_split_fixed(pattern = "_", n=2)
  
  #drop _weight from column names corresponding to LCs and set rownames to output body ids
  colnames(mat_list[[lc]])[-c(1:3)] = tmp[,1]
  rownames(mat_list[[lc]]) = mat_list[[lc]]$output
  
  # remove NAs
  mat_list[[lc]][is.na(mat_list[[lc]])] <- 0
  
  # remove outliers for LC17
  if (lc == "LC17"){ #removing outliers
    lcmat = mat_list[[lc]]
    lcmat = lcmat[!(rownames(lcmat) %in% outrm), !(colnames(lcmat) %in% outrm)]
    v = apply(lcmat[-c(1:3)], 1, var)
    lcmat = lcmat[!(rownames(lcmat) %in% c(names(v[v==0]))),]
    mat_list[[lc]] = lcmat
  }
}
names(mat_list) = LCs_interest

