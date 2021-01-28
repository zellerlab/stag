# Author: Alessio Milanese <milanese.alessio@gmail.com>

# Input:
#  - one multiple sequence alignment (MSA) per marker gene. The MSA is obtained
#    from the function stag align, like:
#       >gene1\t0001010001010000100101000...
#       >gene2\t0000110001010100100101001...
#  - a taxonomy file that describes the taxonomy of the genes:
#       gene1\tBacteria;Firmicutes;...
#
# Output:
#  - a database file (hdf5) that can be used by stag classify

import sys

#===============================================================================
#                          CLASS FOR THE TAXONOMY
#===============================================================================
class Taxonomy:
    # create class -------------------------------------------------------------
    def __init__(self, file_name):
        self.file_name = file_name
        self.child_nodes = dict()
        self.tree_root = "tree_root"
        self.child_nodes[self.tree_root] = set()
        self.last_level_to_genes = dict()
        self.all_gene_ids = list()
        self.number_of_taxonomic_levels = 0
        self.annotation_per_gene = dict()

    # load taxonomy from the defined file --------------------------------------
    def load_from_file(self):
        with open(self.file_name) as tax_in:
            #first_line = next(tax_in).rstrip().split("\t")
            #self.number_of_taxonomic_levels = len(first_line[1].split(";"))
            #tax_in.seek(0)

            for line in tax_in:
                # expected line: gene1\tBacteria\tFirmicutes\t...
                gene_id, tax_levels = line.rstrip().replace("/", "-").split("\t")
                tax_levels = tax_levels.split(";")
                if self.number_of_taxonomic_levels != len(tax_levels) or not tax_levels:
                    # check if number of taxonomic levels matches expectations
                    if self.number_of_taxonomic_levels or not tax_levels:
                        sys.stderr.write("Error: taxonomy record does not have the expected number of taxonomic levels in:\n")
                        sys.stderr.write("  "+line+"\n")
                        sys.exit(1)
                    # if expected numbers were not set, take them from input
                    self.number_of_taxonomic_levels = len(tax_levels)

                # we add the annotation_per_gene:
                self.annotation_per_gene[gene_id] = list(tax_levels)

                # we enter the first level, to the root:
                self.child_nodes[self.tree_root].add(tax_levels[0])
                self.child_nodes.setdefault(tax_levels[0], set())

                # we enter all remaining levels
                for i in range(1, len(tax_levels) - 1):
                    # first we enter that this is a child
                    self.child_nodes[tax_levels[i - 1]].add(tax_levels[i])
                    # and second, we create a node if there is not already
                    self.child_nodes.setdefault(tax_levels[i], set())

                # We add the last level
                self.child_nodes[tax_levels[-2]].add(tax_levels[-1])
                # Finally we add from the last level to the genes ids
                self.last_level_to_genes.setdefault(tax_levels[-1], set()).add(gene_id)
                # and we add it to the list of gene ids
                self.all_gene_ids.append(gene_id)

            self.all_gene_ids.sort() # we sort the list, so that search should be faster


    # make a copy of this taxonomy ---------------------------------------------
    def copy(self):
        temp = Taxonomy(self.file_name)
        temp.child_nodes = dict()
        temp.tree_root = "tree_root"
        for i in self.child_nodes:
            temp.child_nodes[i] = set(self.child_nodes[i])
        temp.last_level_to_genes = dict()
        for i in self.last_level_to_genes:
            temp.last_level_to_genes[i] = set(self.last_level_to_genes[i])
        temp.all_gene_ids = list(self.all_gene_ids)
        temp.number_of_taxonomic_levels = self.number_of_taxonomic_levels
        for i in self.annotation_per_gene:
            temp.annotation_per_gene[i] = list(self.annotation_per_gene[i])
        return temp

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
        if node in self.child_nodes:
            return list(self.child_nodes[node])
        else:
            return None
    # return the last level to genes -------------------------------------------
    def get_last_level_to_genes(self):
        last_level_to_genes_temp = dict()
        for i in self.last_level_to_genes:
            last_level_to_genes_temp[i] = set(self.last_level_to_genes[i])
        return last_level_to_genes_temp

    # check if it is the last node before the genes ----------------------------
    def is_last_node(self, node):
        return node in self.last_level_to_genes
    # find all genes under a given node ----------------------------------------
    # return a list of all genes
    def find_gene_ids(self, node):
        all_leaves = list()
        self.find_leaves_recoursive(node, all_leaves)
        return all_leaves
    def find_leaves_recoursive(self, node, all_leaves):
        if node in self.last_level_to_genes:
            all_leaves.extend(self.last_level_to_genes[node])
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
        to_print = "NODES:\n"
        for i in self.child_nodes:
            to_print = to_print + "   (N):" + i + ": " + str(self.child_nodes[i]) + "\n"
        to_print = to_print + "\nGENES:\n"
        for i in self.last_level_to_genes:
            to_print = to_print + "   (G):" + i + ": " + str(self.last_level_to_genes[i]) + "\n"
        to_print = to_print + "\nLIST GENES:\n" + str(self.all_gene_ids) + "\n"
        to_print = to_print + "\nN LEVELS: " + str(self.number_of_taxonomic_levels) + "\n"
        return to_print
