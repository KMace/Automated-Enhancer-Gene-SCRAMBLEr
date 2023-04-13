from sympy.utilities.iterables import multiset_permutations
import pandas as pd
import numpy as np
import data_initialisation as di
import find_metrics as fm
import region_convolutions as rc
import data_visualisation as dv


def main():
    steps = np.arange(0, 1, 0.5)

    for i in steps:
        print("step number ", i)

        a = np.array([i, i, i, i, i, i, i])
        for p in multiset_permutations(a):
            print(p)

            di.read_config_file(p[0], p[1], p[2], p[3], p[4], p[5], p[6])
        
            gene_annotations = di.read_gene_annotations()
            general_expression_data = di.read_general_expression_data()
            specific_expression_data = di.read_specific_expression_data()
            gene_data = pd.merge(gene_annotations, general_expression_data, on = "Gene_name", how = "inner")
            gene_data = pd.merge(gene_data, specific_expression_data, on = "Gene_name", how = "inner")
            
            del gene_annotations, general_expression_data, specific_expression_data
            
            gene_data = fm.find_gene_sizes(gene_data)
            gene_data = fm.find_interferring_genes(gene_data)
            gene_data = fm.find_search_windows(gene_data)
            
            regulatory_elements = di.read_regulatory_elements()
            enhancers, quiescent_regions = di.clean_regulatory_elements(regulatory_elements)
            del regulatory_elements
            
            enhancer_overlaps = fm.find_element_overlaps_within_search_window(enhancers, gene_data)
            del enhancers, quiescent_regions
            
            gene_data = fm.count_overlaps_per_gene(gene_data, enhancer_overlaps, "Enhancer")
            gene_data = fm.find_nearby_enhancer_densities(gene_data, enhancer_overlaps)
            gene_data = fm.calculate_interest_score(gene_data)

            fm.export_gene_scores_report(p[0], p[1], p[2], p[3], p[4], p[5], p[6])
            
            #gene_data = rc.convolution(gene_data, enhancer_overlaps, "Enhancer")

            #if di.QUIESCENT_CONVOLUTION == True:
            #    gene_data = rc.quiescent_convolution(gene_data, quiescent_overlaps)
            #    del quiescent_overlaps
            #del enhancer_overlaps
            
            # gene_data = rc.find_plateaus(gene_data)
            # rc.export_convolutions(gene_data)
            # rc.export_plateaus(gene_data)
            # dv.gene_report(gene_data)
    
    #enhancer_convolution, recombination_convolution = rc.overlay_convolutions(rc.enhancer_convolution(gene_data, enhancer_overlaps), rc.quiescent_convolution(gene_data, quiescent_overlaps)) 

# Python function to print permutations of a given list
def permutation(lst):
 
    # If lst is empty then there are no permutations
    if len(lst) == 0:
        return []
 
    # If there is only one element in lst then, only
    # one permutation is possible
    if len(lst) == 1:
        return [lst]
 
    # Find the permutations for lst if there are
    # more than 1 characters
 
    l = [] # empty list that will store current permutation
 
    # Iterate the input(lst) and calculate the permutation
    for i in range(len(lst)):
       m = lst[i]
 
       # Extract lst[i] or m from the list.  remLst is
       # remaining list
       remLst = lst[:i] + lst[i+1:]
 
       # Generating all permutations where m is first
       # element
       for p in permutation(remLst):
           l.append([m] + p)
    return l
 
 
if __name__ == "__main__":
    main()