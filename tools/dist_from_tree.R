library(ape)
args <- commandArgs(TRUE)
#tree_file_path="16S.phylip_phyml_tree.txt"
tree_file_path=args[1]
if (file.exists(tree_file_path))
{
	newick_tree <- read.tree(tree_file_path)
	number_of_tips <- length(newick_tree$tip)
	d_n <- dist.nodes(newick_tree)[1:number_of_tips,1:number_of_tips]
	rownames(d_n) <- newick_tree$tip
	colnames(d_n) <- newick_tree$tip

	# print distannes to console
	cat(number_of_tips, "\n")
	for (i in 1:number_of_tips )
	{
		cat(newick_tree$tip[i],d_n[i,], "\n", sep="\t")
	}
}
