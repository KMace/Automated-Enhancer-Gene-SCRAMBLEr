import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

import data_initialisation as di

def define_step_function_of_element_overlaps_within_search_window(gene_data, overlaps, element_type):
    
    #X and Y coordinates are generated for step functions representing elements within search window
    
    print("Generating step function for elements within search window...")
    
    gene_data["Search_window_start"] = gene_data["Search_window_start"].astype("int")
    gene_data["Search_window_end"] = gene_data["Search_window_end"].astype("int")

    gene_data[(element_type + "_step_function_x")] = [np.empty(0, dtype = float)] * len(gene_data)
    gene_data[(element_type + "_step_function_y")] = [np.empty(0, dtype = float)] * len(gene_data)

    gene_data = gene_data.sort_values("Interest_score", ascending = False).reset_index(drop = True)

    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
        
        print("Generating step function of " + element_type.lower() + "s for " + gene["Gene_name"] + " (" + str(index + 1) + " of " + str(di.CONVOLUTION_LIMIT) + ")...")

        gene_specific_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        
        step_function_x = np.arange(gene["Search_window_start"], gene["Search_window_end"])
        step_function_y = np.zeros((gene["Search_window_end"] - gene["Search_window_start"]), dtype = int)
        
        for _, overlap in gene_specific_overlaps.iterrows():
            
            single_step_x = np.where(np.logical_and(overlap["Start"] <= step_function_x, step_function_x <= overlap["End"]), 1, 0)
            step_function_y = np.where(single_step_x == 1, 1, step_function_y)
        
        gene_step_x = np.where(np.logical_and(gene["Gene_start"] <= step_function_x, step_function_x <= gene["Gene_end"]), 1, 0)
        step_function_y = np.where(gene_step_x == 1, 1, step_function_y)
        
        gene_data.at[index, (element_type + "_step_function_x")]  = step_function_x
        gene_data.at[index, (element_type + "_step_function_y")] = step_function_y
        
    return gene_data
    
def convolve_step_function_to_average_windowed_density(gene_data, element_type):

    #X and Y coordinates are generated for convolution of element step function with chosen kernel

    print("Converting step functions to convolved average windowed density signal...")

    gene_data[(element_type + "_convolution_x")] = [np.empty(0, dtype = float)] * len(gene_data)
    gene_data[(element_type + "_convolution_y")] = [np.empty(0, dtype = float)] * len(gene_data)
    
    gene_data = gene_data.sort_values("Interest_score", ascending = False).reset_index(drop = True)

    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
        
        kernel = get_kernel(di.ENHANCER_KERNEL_SHAPE, int((di.RELATIVE_ENHANCER_KERNEL_SIZE * (gene["Search_window_end"] - gene["Search_window_start"]))), int(di.RELATIVE_ENHANCER_KERNEL_SIGMA * (gene["Search_window_end"] - gene["Search_window_start"])))
        convolution_y = np.convolve(kernel, gene[(element_type + "_step_function_y")])
        convolution_x = (np.arange((gene["Search_window_start"] - (len(kernel) // 2)), (gene["Search_window_start"] - (len(kernel) // 2) + len(convolution_y))))

        gene_data.at[index, (element_type + "_convolution_x")] = convolution_x
        gene_data.at[index, (element_type + "_convolution_y")] = convolution_y
        
    return gene_data

def get_kernel(kernel_shape, size, sigma):
    
    #Kernel is generated as numpy array depending on desired shape and size
    
    if kernel_shape == "flat":
        
        kernel = np.ones(size)
        
    
    elif kernel_shape == "guassian":
        
        kernel = np.zeros(size)
        np.put(kernel, (size // 2), 10)
        kernel = gaussian_filter1d(kernel, sigma)
        
    else:
        
        raise Exception("Kernel shape is neither Flat nor Guassian")
        
    return kernel

def trim_convolutions(gene_data, element_type):
    
    #Removes parts of convolution that overlap the edges of the step function, horizontally
    
    print("Trimming convolutions...")
    
    

def combine_convolutions(enhancer_convolution, quiescent_convolution):
    
    #Eventually will be used to add convolutions together
    
    print("Merging convolutions...")

    enhancer_convolution = np.negative(enhancer_convolution)
    print(enhancer_convolution)
    combined_convolution = np.add((enhancer_convolution * di.ENHANCER_CONVOLUTION_WEIGHT), (quiescent_convolution * di.QUIESCENT_CONVOLUTION_WEIGHT))
    
    return combined_convolution
    
def export_convolutions(gene_data):
    
    #Coordinates of convolutions are exported to wig file, for each gene
    
    print("Exporting enhancer density convolutions to wig file...")
    
    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
    
        with open((di.RESULTS_DIRECTORY + gene["Gene_name"] + "_convolutions.wiggle"), "w") as f:
            
            f.write("fixedStep chrom=chr" + gene["Chromosome"] + " start=" + str(gene["Enhancer_convolution_x"][0]) +" step=1")
            f.write("\n")
    
        convolution_signal = pd.DataFrame({"Convolution_signal" : gene_data.loc[index, "Enhancer_convolution_y"]})
        convolution_signal.to_csv((di.RESULTS_DIRECTORY + gene["Gene_name"] + "_convolutions.wig"), sep = "\t", index = False, mode = "a", header = False)
    
def find_plateaus(gene_data):
    
    #find_plateaus takes convolved coordinates, and applies a threshold to
    #separate the search window into regions based on the y-value of each
    #convolved base.
    
    gene_data["Plateau_coordinates"] = ""
    gene_data["Plateau_starts"] = ""
    gene_data["Plateau_ends"] = ""
    
    gene_data = gene_data.sort_values("Interest_score", ascending = False)
    
    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
        
        print("Finding plateaus for gene " + gene["Gene_name"] + " (" + str(index + 1) + " of " + str(di.CONVOLUTION_LIMIT) + ")...")
        
        convolved_x = gene["Enhancer_convolution_x"]
        convolved_y = gene["Enhancer_convolution_y"]
        convolved_y = np.append(convolved_y, 0)
        boolean_below_threshold = convolved_y < di.PLATEAU_THRESHOLD
        boolean_below_threshold = np.concatenate((boolean_below_threshold[:(gene["Gene_start"] - convolved_x[0])], np.full((gene["Gene_end"] -  gene["Gene_start"]), False), boolean_below_threshold[(gene["Gene_end"] - convolved_x[0]):]))
        boolean_below_threshold = (boolean_below_threshold[:-1] != boolean_below_threshold[1:])
        plateau_coordinates = convolved_x[boolean_below_threshold]
        plateau_coordinates = np.concatenate([[convolved_x[0]], plateau_coordinates, [convolved_x[-1]]])
        
        gene_data.at[index, "Plateau_coordinates"] = plateau_coordinates
        gene_data.at[index, "Plateau_starts"] = plateau_coordinates[::2]
        gene_data.at[index, "Plateau_ends"] = plateau_coordinates[1::2]
        
        gene_data.drop(["Plateau_coordinates"], axis = 1)
        
    return gene_data
    
def export_plateaus(gene_data):
    #export_plateaus saves plateaus associated with each gene as a bed file
    
    print("Exporting plateaus to bed file...")
    
    with open((di.RESULTS_DIRECTORY + "plateaus.bed"), "w") as f:
            f.write("Chromosome Start	End	Gene_name")
            f.write("\n")

    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
        
        plateau_regions = pd.DataFrame({"Start" : gene_data.loc[index, "Plateau_starts"], "End" : gene_data.loc[index, "Plateau_ends"]})
        plateau_regions["Gene_name"] = gene["Gene_name"]
        plateau_regions["Chromosome"] = "chr" + gene["Chromosome"]
        plateau_regions["Strand"] = gene["Strand"]
        
        #plateau_regions = ss.find_fasta(plateau_regions)
        
        plateau_regions.to_csv((di.RESULTS_DIRECTORY + "plateaus.bed"), sep = "\t", index = False, columns = ["Chromosome", "Start", "End", "Gene_name"], mode = "a", header = False)