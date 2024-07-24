import pandas as pd
import numpy as np

# t function from test script
# First lets create a function to calculate the t-test
def calculate_t(data_1, data_2, n1=5, n2=5):
    mean1 = np.mean(data_1)
    mean2 = np.mean(data_2)
    var1 = np.std(data_1)
    var2 = np.std(data_2)
    
    # pooled standard deviation
    s = ((var1*(n1 - 1) + var2*(n2 - 1))/ (n1 + n2 - 2))**(1/2)
    
    return (mean1 - mean2)/ (s* (((1/n1) + (1/n2))**(1/2)))

# permutation test for a lean/obese version from test script
# ok now we are going to create a function that runs a permutation test for a single transcript
def permutation(data_lean, data_obese, num_perm=10000):
    # first calculate the baseline observed without 
    t_baseline = calculate_t(data_lean, data_obese)
    
    # lets create a counter that is incremented if success
    counter = 0
    
    # now lets loop through the permutations
    for i in range(num_perm):
        # first task: reshuffle data between lean and obese arrays
        order = np.random.permutation(np.concatenate((data_lean, data_obese)))
        new_data_lean = np.array(order[:5])
        new_data_obese = np.array(order[5:])
        
        trial_t = calculate_t(new_data_lean, new_data_obese)
        if (trial_t >= t_baseline and t_baseline > 0) or (trial_t <= t_baseline and t_baseline < 0):
            counter += 1
    
    return counter/num_perm, num_perm

# create a function to return data with the good columns
def good_columns(filename, tissue_name, time="4wk", strain="B6"):
    # columns we want to keep in pandas
    if (time, strain) not in [("4wk", "B6"), ("4wk", "BTBR"), ("10wk", "B6"), ("10wk", "BTBR")]:
        raise Exception("Not valid combo.")
                    
    # now lets get the lean/obese data for either 4wk, 10wk, b6, or btbr
    columns = ['a_substance_id', 'a_gene_id', 'GeneSymbol_1', 'gene_1']
    for obese in ["lean", "ob"]:
        for i in range(1, 6):
            columns.append(tissue_name + "." + strain + "." + obese + "." + time + "." + str(i))
            
    data = pd.read_csv(filename, usecols=columns)
    return data

def file_names(tissue_name):
    file_name = "TimeCourseData.mlratio.{}.csv".format(tissue_name)
    return tissue_name, file_name