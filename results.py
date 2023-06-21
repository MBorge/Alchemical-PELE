import os
import csv
import pandas as pd
import numpy as np

# Function to get the mean of the first quartile from a lists of numbers
def mean_first_quartile(numbers):
    q1 = np.quantile(numbers, 0.25)
    lower_values = [n for n in numbers if n <= q1]
    # print(np.mean(lower_values))
    return np.mean(lower_values)


# Function to get the binding energies from the report files in each subdirectory
def get_bind_energy():
    output_dir = "output"

    numbers_by_subdir = {}

    # Iterate through each subdirectory in the output directory
    for subdir in os.listdir(output_dir):
        subdir_path = os.path.join(output_dir, subdir)
        if os.path.isdir(subdir_path) and subdir.isdigit():
            numbers = []
            # Iterate through each file in the subdirectory
            for file in os.listdir(subdir_path):
                # Check if the file is a report file
                if file.startswith("report"):
                    report_path = os.path.join(subdir_path, file)
                    with open(report_path, "r") as f:
                        lines = f.readlines()
                        # Iterate through each line in the report file
                        for line in lines[1:]:
                            values = line.split()
                            numbers.append(float(values[4]))
            numbers_by_subdir[int(subdir)] = numbers

    # Sort the dictionary by keys
    values = dict(sorted(numbers_by_subdir.items()))
    # print("All the binding energies:")
    # print(values)
    # print(" ")
    return(values)

# Function to get the internal energies from the report files in each subdirectory
def get_internal_energy():
    output_dir = "output"

    numbers_by_subdir = {}

    # Iterate through each subdirectory in the output directory
    for subdir in os.listdir(output_dir):
        subdir_path = os.path.join(output_dir, subdir)
        if os.path.isdir(subdir_path) and subdir.isdigit():
            numbers = []
            # Iterate through each file in the subdirectory
            for file in os.listdir(subdir_path):
                # Check if the file is a report file
                if file.startswith("report"):
                    report_path = os.path.join(subdir_path, file)
                    with open(report_path, "r") as f:
                        lines = f.readlines()
                        # Iterate through each line in the report file
                        for line in lines[1:]:
                            values = line.split()
                            numbers.append(float(values[5]))
            numbers_by_subdir[int(subdir)] = numbers

    # Sort the dictionary by keys
    values = dict(sorted(numbers_by_subdir.items()))
    # print("All the internal energies:")
    # print(values)
    # print(" ")
    return(values)

# Function to get the strain energies with the reference in each epoch
def strain_energy_epoque(d):
    for lista in d.values():
        minim = min(lista)
        for i in range(len(lista)):
            if lista[i] > 0:
                lista[i] = abs(minim-lista[i])
            elif lista[i] < 0:
                lista[i] = abs(minim-abs(lista[i]))
    # print("All the strain energies:")
    # print(d)
    # print(" ")
    return d

# Function to get the strain energies with the global minimum
def strain_energy_global(d):
    minim = min(numero for lista in d.values() for numero in lista)
    for lista in d.values():
        for i in range(len(lista)):
            if lista[i] > 0:
                lista[i] = abs(minim-lista[i])
            elif lista[i] < 0:
                lista[i] = abs(minim-abs(lista[i]))
    # print("All the strain energies:")
    # print(d)
    # print(" ")
    return d

# Function to get the sum of the binding and the strain energies
def final_energy(dict1, dict2):
    result = {}
    for k in dict1.keys():
        result[k] = [a + b for a, b in zip(dict1[k], dict2[k])]
    # print("All the final energies:")
    # print(result)
    # print(" ")
    return result

# Function to get the minimum values from a dictionary
def get_min_values(numbers_by_subdir):
    min_values_by_subdir = {}
    for subdir, numbers in numbers_by_subdir.items():
        min_value = min(numbers)
        min_values_by_subdir[subdir] = min_value
    
    # print("Minimum values:")
    # print(min_values_by_subdir)
    # print(" ")
    return min_values_by_subdir

# Function to get the mean of the first quartile from a dictionary  
def get_mean_by_subdir(numbers_by_subdir):
    mean_by_subdir = {}
    for subdir, numbers in numbers_by_subdir.items():
        mean = mean_first_quartile(numbers)
        mean_by_subdir[subdir] = mean
    # print("Mean of the first quartile:")
    # print(mean_by_subdir)
    # print(" ")
    return mean_by_subdir

# Function to write the mean and minimum values to a CSV file
def write_csv_from_dict(bind_mean_dict, bind_min_dict, sum_epoc_mean_dict, sum_epoc_min_dict, sum_glob_mean_dict, sum_glob_min_dict):
    # Create the "results" directory if it doesn't exist
    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Write the data to a CSV file in the "results" directory
    filepath = os.path.join(results_dir, "lambdas_results.csv")
    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['lambda', 'binding_mean_low', 'binding_minimum', "sum_epoc_mean_low", "sum_epoc_minimum", "sum_glob_mean_low", "sum_glob_minimum"])
        subdirs = set(bind_mean_dict.keys()) | set(bind_min_dict.keys())
        for subdir in subdirs:
            low_mean = bind_mean_dict.get(subdir)
            minimum = bind_min_dict.get(subdir)
            sum_epoc_low_mean = sum_epoc_mean_dict.get(subdir)
            sum_epoc_minimum = sum_epoc_min_dict.get(subdir)
            sum_glob_low_mean = sum_glob_mean_dict.get(subdir)
            sum_glob_minimum = sum_glob_min_dict.get(subdir)
            writer.writerow([subdir, low_mean, minimum, sum_epoc_low_mean, sum_epoc_minimum, sum_glob_low_mean, sum_glob_minimum])

# Function to calculate the total difference of each column and write it to a txt file
def relative_energy():
    # Path to the lambda_results.csv file
    results_path = "./results/lambdas_results.csv"

    # Read the file as a Pandas dataframe
    df = pd.read_csv(results_path)

    # Ignore the first column
    df = df.iloc[:, 1:]

    # Initialize a string of comma-separated column names
    col_names = ','.join(df.columns) + '\n'

    # Initialize the list of differences
    total_diff = []

    # Iterate through each column and calculate the sum of differences
    for col in df.columns:
        values = df[col].values
        diff = values[-1] - values[0]
        total_diff.append(diff)

    # Write the string of column names and differences to a txt file
    with open("./results/total_diff.txt", "w") as f:
        f.write(col_names)
        f.write(','.join(map(str, total_diff)))


def main():
    bind_energy = get_bind_energy()
    strain_epoc = strain_energy_epoque(get_internal_energy())
    strain_global = strain_energy_global(get_internal_energy())
    sum_epoc_energy = final_energy(bind_energy, strain_epoc)
    sum_global_energy = final_energy(bind_energy, strain_global)
    bind_mean = get_mean_by_subdir(bind_energy)
    bind_min = get_min_values(bind_energy)
    sum_epoc_mean = get_mean_by_subdir(sum_epoc_energy)
    sum_epoc_min = get_min_values(sum_epoc_energy)
    sum_global_mean = get_mean_by_subdir(sum_global_energy)
    sum_global_min = get_min_values(sum_global_energy)
    write_csv_from_dict(bind_mean, 
                        bind_min,
                        sum_epoc_mean, 
                        sum_epoc_min,
                        sum_global_mean,
                        sum_global_min)
    relative_energy()



if __name__ == "__main__":
    main()