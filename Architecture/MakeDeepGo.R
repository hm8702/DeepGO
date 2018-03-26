#################################################################
#
#	Make GO architecure 
#	2018-3-18
#	ccma
###################################################

library(GO.db)
library(org.Hs.eg.db)

tmp = org.Hs.egGO2ALLEGS
mapped_genes = mappedkeys(tmp)
g2e = as.list(tmp[mapped_genes])

tmp = org.Hs.egSYMBOL
mapped_genes = mappedkeys(tmp)
e2s = as.list(tmp[mapped_genes])

children <- as.list(GOBPCHILDREN)
offspring <- as.list(GOBPOFFSPRING)
parents <- as.list(GOBPPARENTS)
ancestor <- as.list(GOBPANCESTOR)

### function get additional genes of the tier
get_genes = function(go_list){
	upper = unique(unlist(sapply(go_list, FUN=function(x){parents[x]})))
	upper_genes = unique(unlist(sapply(upper, FUN=function(x){g2e[x]})))
	this_genes = unique(unlist(sapply(go_list, FUN=function(x){g2e[x]})))
	add_genes = setdiff(upper_genes, this_genes)
}
##########################################################
#	Make tier 1 GO modules only have genes no offspring
#########################################################
check_length = function(x){var_name = paste("t", as.character(x), "_go", sep=""); return(length(get(var_name)))}
t1_go = intersect(names(which(is.na(offspring))), names(g2e))
t0_gene = unique(unlist(g2e[t1_go]))

t1_add_genes = get_genes(t1_go)

############################################
#	Make GO tiers. GO tiers may corss tiers
###########################################
i = 2
ancestor_go = unique(unlist(parents)) ## 17492
go_tiers = list()	
go_tiers[["t1_go"]] = t1_go
while(T){
	pre_go = go_tiers[[paste("t", (i-1), "_go", sep="")]]
	this_go_name = paste("t", i, "_go", sep="") 
	tmp_go = unique(unlist(parents[pre_go]))
	this_go = setdiff(tmp_go, unique(unlist(ancestor[tmp_go])))
	#this_go = tmp_go
	if(length(this_go) >=1){
		go_tiers[[this_go_name]] = this_go
		i = i+1
	}else{
		break
	}
}
length(unique(unlist(go_tiers)))    ### 13297
length(unlist(go_tiers))		### 13297


######################################################################
#	add gene to each tier if it has not been considered. Genes won't cross tiers
#####################################################################

all_genes = unique(unlist(g2e))   ## 18906
rest_genes = all_genes

gene_tiers = list()
gene_tiers[["t0_genes"]] = t0_gene
for(i in 1:(length(go_tiers)-2)){
	needed_genes = unique(unlist(g2e[go_tiers[[i+1]]]))
	considered_genes = unique(unlist(g2e[unique(unlist(offspring[go_tiers[[i+1]]]))]))
	this_gene_name = paste("t", i, "_genes", sep="")
	this_genes = setdiff(needed_genes, considered_genes)
	this_genes  = intersect(this_genes, rest_genes)
	gene_tiers[[this_gene_name]] = this_genes
	rest_genes = setdiff(rest_genes, this_genes)
}
length(unique(unlist(gene_tiers)))    ### 10419
length(unlist(gene_tiers))         ### 10419


