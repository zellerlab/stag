import sys
import csv

class Taxonomy:
    TREE_ROOT = "tree_root"
    def __init__(self, file_name):
        self.file_name = file_name
        self.child_nodes = {Taxonomy.TREE_ROOT: set()}
        self.tree_root = Taxonomy.TREE_ROOT
        self.last_level_to_genes = dict()
        self.all_gene_ids = list()
        self.number_of_taxonomic_levels = 0
        self.annotation_per_gene = dict()

    def _add_child(self, parent, child):
        self.child_nodes.setdefault(parent, set()).add(child)

    # load taxonomy from the defined file --------------------------------------
    def load_from_file(self):
        with open(self.file_name,"r") as tax_in:
            for line_i, row in enumerate(csv.reader(tax_in, delimiter="\t")):
                # what is the purpose of this? removing file path from genes?
                # vals = line.rstrip().replace("/","-").split("\t")
                try:
                    gene, lineage = row
                except:
                    raise ValueError(f"line {line_i} is not properly formatted, expecting gene\tlineage.. :\n{row}")
                lineage = lineage.split(";")
                if line_i == 0:
                    self.number_of_taxonomic_levels = len(lineage)
                elif self.number_of_taxonomic_levels != len(lineage):
                    raise ValueError(f"line {line_i}'s lineage=depth {len(lineage)} does not match previous lineage-depth ({self.number_of_taxonomic_levels})\n{lineage}")
 
                self.annotation_per_gene[gene] = list(lineage)
                self._add_child(self.tree_root, lineage[0])
                for tax_i, child in enumerate(lineage):
                    parent = self.tree_root if tax_i == 0 else lineage[tax_i - 1]
                    self._add_child(parent, child)
                self.last_level_to_genes.setdefault(lineage[-1], set()).add(gene)
                self.all_gene_ids.append(gene)
    
            self.all_gene_ids.sort() # we sort the list, so that search should be faster

    # make a copy of this taxonomy ---------------------------------------------
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

    # return the classification of one gene
    def extract_full_tax_from_gene(self, gene_id):
        return self.annotation_per_gene[gene_id]

    # return number of levels --------------------------------------------------
    def get_n_levels(self):
        return self.number_of_taxonomic_levels

    # return the root id -------------------------------------------------------
    def get_root(self):
        return self.tree_root

    # find children of a node --------------------------------------------------
    def find_children_node(self, node):
        return list(self.child_nodes.get(node, list()))

    # return the last level to genes -------------------------------------------
    def get_last_level_to_genes(self):
        return {node: set(genes) for node, genes in self.last_level_to_genes.items()}

    # check if it is the last node before the genes ----------------------------
    def is_last_node(self, node):
        return self.last_level_to_genes.get(node) is not None 

    # find all genes under a given node ----------------------------------------
    # return a list of all genes
    def find_gene_ids(self, node):
        all_leaves = list()
        self.find_leaves_recoursive(node, all_leaves)
        return all_leaves
    def find_leaves_recoursive(self, node, all_leaves):
        genes = self.last_level_to_genes.get(node)
        if genes:
            all_leaves.extend(genes)
        else:
            for c in self.child_nodes[node]:
                self.find_leaves_recoursive(c, all_leaves)

    # function to remove nodes (and genes underneath), given a list of nodes ---
    # it returns the gene ids that were removed
    def remove_clades(self, node_list):
        # remove all clades under
        list_removed_genes = list()
        for n in node_list:
            self.remove_clade_iter(n, list_removed_genes)
        # remove all clades on top
        for n in node_list:
            # now need to remove from the higher level in child_nodes
            for i in self.child_nodes:
                self.child_nodes[i].discard(n) # discard does not raise a KeyError.
        # if it was the only child, then we should remove also at higher level
        self.remove_unused_branches()
        return list(list_removed_genes)

    def remove_clade_iter(self, node, list_removed_genes):
        if node in self.last_level_to_genes:
            # we arrived at the end of the tree, we remove the genes, but first:
            # add to the set of removed genes
            list_removed_genes.extend(self.last_level_to_genes[node])
            # remove the genes from the gene list
            self.all_gene_ids = [e for e in self.all_gene_ids if e not in self.last_level_to_genes[node]]
            # and, finally, remove the node from the last_level_to_genes dict
            self.last_level_to_genes.pop(node,None)
        else:
            try:
                for n in self.child_nodes[node]:
                    self.remove_clade_iter(n, list_removed_genes)
                # remove from child_nodes
                self.child_nodes.pop(node,None)
            except KeyError as e:
                sys.stderr.write("WARNING: key not present when removing a clade ["+str(e)+"]\n\n")


    def remove_unused_branches(self):
        removed_any = False # this becomes True if we remove any node from
                            # child_nodes, in which case we re-run remove_unused_branches
        list_to_remove = list()
        for i in self.child_nodes:
            if len(self.child_nodes[i]) == 0:
                removed_any = True
                list_to_remove.append(i)
                # remove from taxonomy at higher level
                for j in self.child_nodes:
                    self.child_nodes[j].discard(i)
        # remove nodes that are empty from child_nodes
        for n in list_to_remove:
            self.child_nodes.pop(n,None)
        # call remove_unused_branches again if necessary
        if removed_any:
            self.remove_unused_branches()

    # function to remove genes from a list -------------------------------------
    def remove_genes(self, gene_list):
        # remove the genes from the gene list
        self.all_gene_ids = [e for e in self.all_gene_ids if e not in gene_list]
        # remove the genes from last_level_to_genes
        for g in gene_list:
            for node in self.last_level_to_genes:
                self.last_level_to_genes[node].discard(g)
        # Check if all the genes from one clade are removed, and hence we should
        # remove that clade
        list_to_remove = list()
        for node in self.last_level_to_genes:
            if len(self.last_level_to_genes[node]) == 0:
                list_to_remove.append(node)
        self.remove_clades(list_to_remove)

    # function that returns all nodes at one level, ----------------------------
    # as a dictionary of the parent nodes
    def find_tax_level_iter(self, current_node, current_level, result):
        if current_node in self.child_nodes:
            for n in self.child_nodes[current_node]:
                result[n] = current_level+1
                self.find_tax_level_iter(n, current_level+1, result)
    def find_node_level(self, tax_level_find):
        # find tax level for each node
        tax_level = dict()
        tax_level[self.tree_root] = 0
        self.find_tax_level_iter(self.tree_root,0,tax_level)
        # select only the one from the correct level
        res = dict()
        for n in self.child_nodes:
            if tax_level[n] == tax_level_find:
                res[n] = set(self.child_nodes[n])
        return res

    # print the values in the taxonomy class -----------------------------------
    def __str__(self):
        string = ["NODES:"]
        string.extend(f"   (N):{node}: {children}" for node, children in self.child_nodes.items())
        string.extend(("", "GENES:"))
        string.extend(f"   (G):{node}: {genes}" for node, genes in self.last_level_to_genes.items()) 
        string.extend(("", "LIST GENES:", str(self.all_gene_ids)))
        string.extend(("", f"N LEVELS: {self.number_of_taxonomic_levels}", ""))
        return "\n".join(string)
