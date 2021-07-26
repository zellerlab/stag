import csv
import logging

class Taxon:
    def __init__(self, level=0, parent=None, label=None):
        self.level = level
        self.label = label if label else Taxonomy.TREE_ROOT
        self.children = dict()
        self.genes = set()
        self.species_nodes = set()
        self.parent = parent
    def add_child(self, child):
        self.children.setdefault(child.label, child)
    def add_species_descendant(self, node):
        self.species_nodes.add(node)
    def add_gene(self, gene):
        self.genes.add(gene)
    def is_leaf(self):
        return not self.children


class Taxonomy(dict):
    TREE_ROOT = "tree_root"
    def __init__(self, fn=None):
        self[self.TREE_ROOT] = Taxon()
        self.n_taxlevels = 0
        self.gene_lineages = dict()
        self.fn = fn
        if self.fn:
            self._load_taxonomy()

    def _check_lineage_depth(self, lineage, line_no):
        lineage = lineage.replace("/", "-").split(";") # issue10
        if len(lineage) < self.n_taxlevels:
            raise ValueError(f"Line {line_no}: Taxonomy record does not have the expected number of taxonomic levels\n{lineage}")
        self.n_taxlevels = len(lineage)
        return lineage

    def _load_taxonomy(self):
        for line_no, (gene, lineage) in enumerate(csv.reader(open(self.fn), delimiter="\t"), start=1):
            parent = self[self.TREE_ROOT]
            lineage = self._check_lineage_depth(lineage, line_no)
            last_level = len(lineage) - 1
            for level, taxon in enumerate(lineage):
                if level > 0:
                    parent = node
                node = self.setdefault(taxon, Taxon(level=level, parent=parent, label=taxon))
                parent.add_child(node)
                if level < last_level:
                    node.add_species_descendant(lineage[-1])
            node.add_gene(gene)
            self.gene_lineages[gene] = lineage

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

    def extract_full_tax_from_gene(self, gene):
        return self.gene_lineages.get(gene)

    def get_n_levels(self):
        return self.n_taxlevels

    def get_root(self):
        return self.TREE_ROOT

    def find_children_node(self, node):
        return list(self.get(node, Taxon()).children.keys())

    def get_last_level_to_genes(self):
        return {node: set(node.genes) for node in self.values() if node.genes}

    def is_last_node(self, node):
        return self.get(node, Taxon()).is_leaf()

    def find_gene_ids(self, node=None):
        nodes = list(self[Taxonomy.TREE_ROOT].children.keys()) if (node is None or node == self.get_root()) else [node]

        if any(self.get(node) is None for node in nodes):
            raise ValueError(f"Nodelist invalid: {nodes}")

        genes = set()
        for node in nodes:
            if self[node].level == self.n_taxlevels - 1:
                genes.update(self[node].genes)
            else:
                for descendant in self[node].species_nodes:
                    genes.update(self.get(descendant, Taxon()).genes)

        return list(genes)

    def find_gene_ids_old(self, node=None):
        genes = set()
        nodes = [self[node if node else self.TREE_ROOT]]
        while nodes:
            node = nodes.pop(0)
            nodes.extend(node.children.values())
            genes.update(node.genes)
        return list(genes)

    def remove_clades(self, nodes):
        removed_genes = set()   
        for node in nodes:
            stack = [node]
            while stack:
                node2 = self[stack.pop()]
                removed_genes.update(node2.genes)
                stack.extend(node2.children)
                if node2.parent:
                    node2.parent.children.pop(node2.label, None)
                    self._clean_branch(node2.parent)
                self.pop(node2.label, None)
        for gene in removed_genes:
            self.gene_lineages.pop(gene, None)
        return list(removed_genes)

    def _clean_branch(self, node):
        while True:
            if node.children or not node.parent:
                break
            try:
                self.pop(node.label)
                node.parent.children.pop(node.label)
            except:
                pass
            node = node.parent

    def remove_genes(self, genes):
        empty_nodes = set()
        for gene in genes:
            node = self[self.gene_lineages[gene][-1]]
            node.genes.discard(gene)
            if not node.genes:
                empty_nodes.add(node.label)
        self.remove_clades(empty_nodes)

    def find_node_level(self, tax_level):
        nodes = dict()
        queue = [(self[self.TREE_ROOT], -1)]
        while queue:
            node, level = queue.pop(0)
            if level + 1 == tax_level:
                for child in node.children.values():
                    nodes[child.label] = set(child.children)
            else:
                queue.extend((child, level + 1) for child in node.children.values())
        return nodes

    def get_all_nodes(self, mode=None, get_root=False):
        queue = [(self[self.TREE_ROOT], set())]
        while queue:
            node, siblings = queue.pop(0)
            if node.label != self.get_root() or get_root:
                yield node.label, siblings

            children = set(node.children)
            for child in children:
                siblings = children.difference({child})
                queue.append((self[child], siblings))

    def ensure_geneset_consistency(self, genes):
        genes_in_tree = set(self.find_gene_ids())
        logging.info(f"   CHECK: genes in geneset: {len(genes)}")
        logging.info(f"   CHECK: genes in taxonomy: {len(genes_in_tree)}")

        # check that all genes in the geneset are in the taxonomy
        missing_genes = set(genes).difference(genes_in_tree)
        if missing_genes:
            logging.info(" Error: some genes in the alignment have no taxonomy.")
            for gene in missing_genes:
                logging.info(f"    {gene}")
            raise ValueError("Some genes in the alignment have no taxonomy.\n"
                             "Use the command 'check_input' to find more information.\n")
        else:
            logging.info("   CHECK: check all genes in the alignment have a taxonomy: correct")

        # the taxonomy can have more genes than the geneset, but these need to be removed
        # since the selection of the genes for training and testing is done at taxonomy level
        drop_genes = genes_in_tree.difference(genes)
        if drop_genes:
            n_drop_genes = len(drop_genes)
            self.remove_genes(drop_genes)
        else:
            n_drop_genes = None
        logging.info(f"   CHECK: check genes that we need to remove from the taxonomy: {n_drop_genes}")

        # verify number of genes is consistent between set and taxonomy tree
        genes_in_tree = self.find_gene_ids()
        if len(genes_in_tree) != len(genes):
            msg = "Even after correction, the genes in the taxonomy and the alignment do not agree."
            logging.info(f" Error: {msg.lower()}")
            raise ValueError(msg) 
