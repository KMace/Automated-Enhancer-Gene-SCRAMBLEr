import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import pyranges
import sys
import sklearn as sk
from sklearn.tree import DecisionTreeClassifier

chromosomes_of_interest = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

non_housekeeping_weight = (1 / 4)
enhancer_count_weight = (1 / 250)
enhancer_proportion_weight = (1 / 0.4)
HAP1_expression_weight = (1 / 14)
gene_size_weight = (16788 / 1)

def main():
    read_command_line()
    read_genes(gene_annotation_reference)
    clean_genes()
    read_expression(cell_lines_expression_reference)
    clean_expression()
    read_regulatory_elements(regulatory_elements_reference)
    clean_regulatory_elements()
    find_gene_size()
    merge_annotation_expression()
    proximal_enhancer_search()
    data_exploration()
    gene_scoring()

def read_command_line():
    print("Reading command line...")
    global gene_annotation_reference
    global cell_lines_expression_reference
    global regulatory_elements_reference
    global upstream_search
    global downstream_search
    global cell_line_of_interest
    global results_directory
    
    gene_annotation_reference = sys.argv[1]
    cell_lines_expression_reference = sys.argv[2]
    regulatory_elements_reference = sys.argv[3]
    upstream_search = int(sys.argv[4])
    downstream_search = int(sys.argv[5])
    cell_line_of_interest = sys.argv[6]
    results_directory = sys.argv[7]
    
    

def read_genes(genes_file):
    print("Reading gene annotations file...")
    global genes
    genes = pd.read_csv(genes_file, sep = "\t", names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"], skiprows = 5, dtype = {"Chromosome" : str, "Source" : str, "Type" : str, "Start" : int, "End" : int, "Score" : str, "Strand" : str, "Phase" : str, "Attributes" : str})
    
def read_expression(expression_file):
    print("Reading gene expression file...")
    global expression
    expression = pd.read_csv(expression_file).transpose()

def read_regulatory_elements(regulatory_file):
    print("Reading regulatory elements file...")
    global regulatory_elements
    regulatory_elements = pd.read_csv(regulatory_file, sep = "\t", names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"])

def clean_genes():
    print("Cleaning gene data...")
    global genes
    genes = genes[genes["Chromosome"].isin(chromosomes_of_interest)]
    genes = genes.drop(genes[genes["Type"] != "gene"].index)
    genes["Gene_biotype"] = genes["Attributes"].apply(lambda x : re.findall("gene_biotype \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    genes = genes.drop(genes[genes["Gene_biotype"] != "protein_coding"].index)
    genes["Gene_name"] = genes["Attributes"].apply(lambda x : re.findall("gene_name \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    genes = genes.drop(["Source", "Type", "Score", "Phase", "Attributes", "Gene_biotype"], axis = 1)

def clean_expression():
    print("Cleaning expression data...")
    global expression
    expression = expression.drop(["depmapID", "primary_disease"], axis = 0)
    expression.columns = expression.iloc[0]
    expression = expression.iloc[1:]
    expression = expression.reset_index().rename(columns = {"index" : "Gene_name"})
    expression["Mean"] = expression.loc[:, expression.columns != "Gene_name"].mean(axis = 1)
    expression["Std"] = expression.loc[:, expression.columns != "Gene_name"].std(axis = 1)
    expression = expression[["Gene_name", cell_line_of_interest, "Mean", "Std"]]

def clean_regulatory_elements():
    print("Cleaning regulatory elements data...")
    global regulatory_elements
    regulatory_elements = regulatory_elements.drop(regulatory_elements[regulatory_elements["Type"] != "enhancer"].index)
    regulatory_elements["Enhancer_ID"] = regulatory_elements["Attributes"].apply(lambda x : re.findall("ID=enhancer:(.*?);", x)[0])
    regulatory_elements = regulatory_elements.drop(["Source", "Type", "Score", "Strand", "Phase", "Attributes"], axis = 1)

def merge_annotation_expression():
    print("Merging annotations and expression data...")
    global genes, expression
    genes = pd.merge(genes, expression, on = "Gene_name", how = "inner")

def find_gene_size():
    global genes
    genes["Gene_size"] = genes["End"] - genes["Start"]

def proximal_enhancer_search():
    print("Searching for proximal enhancers...")
    global genes
    global regulatory_elements
    global overlaps
    genes_search = genes
    #genes_search["Start"] = np.where(genes_search["Strand"] == "+", genes_search["Start"] - upstream_search, genes_search["Start"] - downstream_search)
    #genes_search["End"] = np.where(genes_search["Strand"] == "+", genes_search["End"] + downstream_search, genes_search["End"] + upstream_search)
    genes_search["Start"] = genes_search["Start"] - upstream_search
    genes_search["End"] = genes_search["End"] + downstream_search
    genes_search["Search_size"] = genes_search["End"] - genes_search["Start"]
    print(upstream_search)
    print(downstream_search)
    print(genes)
    print(genes_search)
    gene_pr = pyranges.PyRanges(genes)
    print(gene_pr)
    search_pr = pyranges.PyRanges(genes_search)
    print(search_pr)
    #search_pr = search_pr.intersect(gene_pr, strandedness = False, invert = True)
    print(search_pr)
    regulatory_elements_pr = pyranges.PyRanges(regulatory_elements)
    overlaps = search_pr.intersect(regulatory_elements_pr, strandedness = False)
    overlaps_count = overlaps.df.groupby("Gene_name").size().reset_index(name = "Enhancer_count")
    overlaps_size = overlaps.df
    overlaps_size["Enhancer_content"] = overlaps_size["End"] - overlaps_size["Start"]
    overlaps_size = overlaps_size.groupby(["Gene_name", "Search_size"], as_index = False)["Enhancer_content"].sum().reset_index()
    overlaps_size["Enhancer_proportion"] = (overlaps_size["Enhancer_content"] / overlaps_size["Search_size"])
    overlaps_count = pd.merge(overlaps_count, overlaps_size, on = "Gene_name", how = "inner")
    genes = pd.merge(genes, overlaps_count, on = "Gene_name", how = "inner")

def data_exploration():
    print("Visualising data...")
    global relative_expression
    global overlaps
    global genes
    relative_expression = expression.sort_values(by = ["Std"], ascending = False)
    #relative_expression = relative_expression[relative_expression["mean"] < relative_expression["HAP1"]]
    relative_expression["Difference"] = relative_expression[cell_line_of_interest] - relative_expression["Mean"]

    figure, axis = plt.subplots(1, 3, figsize = (18.5, 10.5))
    axis[0].scatter(relative_expression["Std"], relative_expression[cell_line_of_interest], s = 0.3, c = "red")
    axis[0].scatter(relative_expression["Std"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[0].set(xlabel = "Standard deviation of expression within cell lines", ylabel = "Normalised expression")

    axis[1].scatter(relative_expression["Difference"], relative_expression[cell_line_of_interest], s = 0.3, c = "red")
    axis[1].scatter(relative_expression["Difference"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[1].set(xlabel = "Difference between expression within " + cell_line_of_interest + " and mean expression", ylabel = "Normalised expression")

    axis[2].scatter(relative_expression["Std"], relative_expression["Difference"], s = 0.3, c = "green")
    axis[2].set(xlabel = "Difference between expression within " + cell_line_of_interest + " and mean expression", ylabel = "Normalised expression")

    plt.savefig(results_directory + "Significant_" + cell_line_of_interest + "_expression.png")
    plt.close(figure)
    
    figure, axis = plt.subplots(2, 2, figsize = (18.5, 10.5))
    axis[0][0].scatter(genes["Enhancer_count"], genes["Std"], s = 0.3, c = "purple")
    axis[0][0].set(xlabel = "Number of enhancers within search area", ylabel = "Standard deviation of expression within cell lines")
    
    axis[1][0].scatter(genes["Enhancer_count"], genes[cell_line_of_interest], s = 0.3, c = "red")
    axis[1][0].scatter(genes["Enhancer_count"] + 0.5, genes["Mean"], s = 0.3, c = "blue")
    axis[1][0].set(xlabel = "Number of enhancers within search area", ylabel = "Normalised expression")
    
    axis[0][1].scatter(genes["Enhancer_proportion"], genes["Std"], s = 0.3, c = "purple")
    axis[0][1].set(xlabel = "Proportion of search area which is enhancer", ylabel = "Standard deviation of expression within cell lines")
    
    axis[1][1].scatter(genes["Enhancer_proportion"], genes[cell_line_of_interest], s = 0.3, c = "red")
    axis[1][1].scatter(genes["Enhancer_proportion"], genes["Mean"], s = 0.3, c = "blue")
    axis[1][1].set(xlabel = "Proportion of search area which is enhancer", ylabel = "Normalised expression")
    
    plt.savefig(results_directory + "Enhancers_near_genes.png")
    plt.close(figure)

def gene_scoring():
    print("Scoring genes...")
    global genes
    print(genes.columns)
    gene_scores = genes[["Gene_name", "Std", "HAP1", "Enhancer_count", "Enhancer_proportion", "Gene_size"]]
    gene_scores["Interest_score"] = gene_scores["Std"] * non_housekeeping_weight + \
        gene_scores["Enhancer_count"] * enhancer_count_weight + \
        gene_scores["Enhancer_proportion"] * enhancer_proportion_weight + \
        gene_scores["HAP1"] * HAP1_expression_weight + \
        (1 / gene_scores["Gene_size"]) * gene_size_weight
    gene_scores = gene_scores.sort_values("Interest_score", ascending = False)
    print(gene_scores.head(50))

if __name__ == "__main__":
    main()