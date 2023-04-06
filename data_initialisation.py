import sys
import re
import json
import pandas as pd
import numpy as np

import find_metrics as fm

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
    
# read_config_file assigns global variables based on the config file.
def read_config_file():
    with open(sys.argv[1], "r") as config_file:
        settings = json.load(config_file)

    global STD_WEIGHT
    global ANOMALOUS_EXPRESSION_WEIGHT
    global ENHANCER_COUNT_WEIGHT
    global ENHANCER_PROPORTION_WEIGHT
    global CELL_LINE_EXPRESSION_WEIGHT
    global GENE_SIZE_WEIGHT

    global ENHANCER_KERNEL_SHAPE
    global ENHANCER_KERNEL_SIZE_TYPE
    global ABSOLUTE_ENHANCER_KERNEL_SIZE
    global RELATIVE_ENHANCER_KERNEL_SIZE
    global RELATIVE_ENHANCER_KERNEL_SIGMA
    global MIN_ABSOLUTE_ENHANCER_CLUSTER_WIDTH
    global MIN_ENHANCER_CLUSTER_PROMINENCE
    
    global QUIESCENT_KERNEL_SHAPE
    global QUIESCENT_KERNEL_SIZE_TYPE
    global ABSOLUTE_QUIESCENT_KERNEL_SIZE
    global RELATIVE_QUIESCENT_KERNEL_SIZE
    global RELATIVE_QUIESCENT_KERNEL_SIGMA
    global MIN_ABSOLUTE_QUIESCENT_CLUSTER_WIDTH
    global MIN_QUIESCENT_CLUSTER_PROMINENCE

    global SIGMOIDAL_SLOPE
    global SIGMOIDAL_MIDPOINT
    global CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD
    global INTERFERRING_GENE_OVERLAPS
    
    global ENHANCER_CONVOLUTION
    global QUIESCENT_CONVOLUTION
    global ENHANCER_CONVOLUTION_WEIGHT
    global QUIESCENT_CONVOLUTION_WEIGHT
    global PLATEAU_THRESHOLD

    STD_WEIGHT = settings["relative_std_weight"]
    ANOMALOUS_EXPRESSION_WEIGHT = settings["relative_anomalous_expression_weight"]
    ENHANCER_COUNT_WEIGHT = settings["relative_enhancer_count_weight"]
    ENHANCER_PROPORTION_WEIGHT = settings["relative_enhancer_proportion_weight"]
    CELL_LINE_EXPRESSION_WEIGHT = settings["relative_cell_line_expression_weight"]
    GENE_SIZE_WEIGHT = settings["relative_gene_size_weight"]

    ENHANCER_KERNEL_SHAPE = settings["enhancer_kernel_shape"]
    ENHANCER_KERNEL_SIZE_TYPE = settings["enhancer_kernel_size_type"]
    ABSOLUTE_ENHANCER_KERNEL_SIZE = settings["absolute_enhancer_kernel_size"]
    RELATIVE_ENHANCER_KERNEL_SIZE = settings["relative_enhancer_kernel_size"]
    RELATIVE_ENHANCER_KERNEL_SIGMA = settings["relative_enhancer_kernel_sigma"]
    MIN_ABSOLUTE_ENHANCER_CLUSTER_WIDTH = settings["min_absolute_enhancer_cluster_width"]
    MIN_ENHANCER_CLUSTER_PROMINENCE = settings["min_enhancer_cluster_prominence"]

    QUIESCENT_KERNEL_SHAPE = settings["quiescent_kernel_shape"]
    QUIESCENT_KERNEL_SIZE_TYPE = settings["quiescent_kernel_size_type"]
    ABSOLUTE_QUIESCENT_KERNEL_SIZE = settings["absolute_quiescent_kernel_size"]
    RELATIVE_QUIESCENT_KERNEL_SIZE = settings["relative_quiescent_kernel_size"]
    RELATIVE_QUIESCENT_KERNEL_SIGMA = settings["relative_quiescent_kernel_sigma"]
    MIN_ABSOLUTE_QUIESCENT_CLUSTER_WIDTH = settings["min_absolute_quiescent_cluster_width"]
    MIN_QUIESCENT_CLUSTER_PROMINENCE = settings["min_quiescent_cluster_prominence"]

    SIGMOIDAL_SLOPE = settings["sigmoidal_slope"]
    SIGMOIDAL_MIDPOINT = settings["sigmoidal_midpoint"]
    CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD = settings["cell_line_specific_expression_threshold"]
    INTERFERRING_GENE_OVERLAPS = settings["interferring_gene_overlaps"]

    ENHANCER_CONVOLUTION = settings["enhancer_convolution"]
    QUIESCENT_CONVOLUTION = settings["quiescent_convolution"]
    ENHANCER_CONVOLUTION_WEIGHT = settings["enhancer_convolution_weight"]
    QUIESCENT_CONVOLUTION_WEIGHT = settings["quiescent_convolution_weight"]
    PLATEAU_THRESHOLD = settings["plateau_threshold"]


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
