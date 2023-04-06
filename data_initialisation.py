import sys
import re
import json
import pandas as pd
import numpy as np

configFileLocation = sys.argv[1]

# returnConfigFile opens and returns the configuration file containing important
# directories, weights and constants, etc.
def returnConfigFile(configFilePath):
    with open(configFilePath, "r") as file:
        config_file = json.load(file)
        
        return config_file

# returnConfigDirectoryPaths returns the configuration directory paths from the
# configuration file.
def returnConfigDirectoryPaths(configFilePath):
    file = returnConfigFile(configFilePath)
    
    results = file["results_directory"]
    gene_prioritisation_report = file["gene_prioritisation_report_directory"]
    
    return results, gene_prioritisation_report
    
# returnConfigReferencePaths returns gene reference paths from the configuration
# file.
def returnConfigReferencePaths(configFilePath):
    file = returnConfigFile(configFilePath)
    
    genome = file["reference_genome"]
    gene_annotation = file["gene_annotation_reference"]
    regulatory_elements = file["regulatory_elements_reference"]
    general_expression = file["general_expression_by_cell_line_reference_path"]
    specific_expression = file["specific_expression_by_cell_line_reference_path"]
    
    return genome, gene_annotation, regulatory_elements, general_expression, specific_expression

# returnConfigSearchMetrics returns the genome search metrics from the
# configuration file.
def returnConfigSearchMetrics(configFilePath):
    file = returnConfigFile(configFilePath)
    
    search_type = file["search_type"]
    search_within_gene = file["search_within_gene"] # rename this - is it an internal search window??
    upstream_search = file["upstream_search"]
    downstream_search = file["downstream_search"]
    
    return search_type, search_within_gene, upstream_search, downstream_search

# This is inappropriately named: you are creating a dataframe; more than just
# reading gene annotations.        

# read_gene_annotations creates a dataframe based on the gene annotations.
def read_gene_annotations(): 
    _, gene_annotation_path, _, _, _ = returnConfigReferencePaths(configFileLocation)
    
    gene_annotations = pd.read_csv(
        gene_annotation_path, sep = "\t",
        skiprows = 5,
        names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", 
                 "Attributes"],
        dtype = {"Chromosome" : str, "Source" : str, "Type" : str, "Start" : int, 
                 "End" : int, "Score" : str, "Strand" : str, "Phase" : str, "Attributes" : str}
    )
    
    return gene_annotations
