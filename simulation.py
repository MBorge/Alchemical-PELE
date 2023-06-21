import shutil
import os
from AdaptivePELE import adaptiveSampling
#import results
#import work
#import bar_analysis

# Function to copy the source files to the specified destination directories
def replace_files(num):
    # Define the source file paths using string formatting
    hybz_source_file = f"alchemical_files/hybz_{num}"
    ligandParams_source_file = f"alchemical_files/ligandParams_{num}.txt"
    HYB_rotamers_source_file = f"alchemical_files/HYB.rot.assign{num}"

    # Use shutil to copy the source files to the destination directories
    shutil.copy(hybz_source_file, "DataLocal/Templates/OpenFF/Parsley/hybz")
    shutil.copy(ligandParams_source_file, "DataLocal/OBC/ligandParams.txt")
    shutil.copy(HYB_rotamers_source_file,"DataLocal/LigandRotamerLibs/HYB.rot.assign")


# Function to count the number of hybz files in the alchemical_files directory
def count_files():
    # Change the current working directory to "alchemical_files"
    os.chdir("alchemical_files")
    count = 0
    # Count the number of files that start with "hybz"
    for file in os.listdir():
        if file.startswith("hybz"):
            count += 1
    # Change back to the parent directory
    os.chdir("..")
    # Return the count of hybz files
    return count


# Function to update the value of an iteration in the "adaptive_1.8.conf" file
def update_iterations(iteration):
    # Open the "adaptive_1.8.conf" file for reading
    with open("adaptive_1.8.conf", "r") as f:
        lines = f.readlines()
    # Update the 17th line with the desired iteration value
    lines[17] = '         "iterations": ' + str(iteration) + "," + "\n"
    # Open the "adaptive_1.8.conf" file for writing
    with open("adaptive_1.8.conf", "w") as f:
        f.writelines(lines)


# Function to update the value of the steps in the "adaptive_1.8.conf" file
def update_steps(steps):
    # Open the "adaptive_1.8.conf" file for reading
    with open("adaptive_1.8.conf", "r") as f:
        lines = f.readlines()
    # Update the 18th line with the desired steps value
    lines[18] = '         "peleSteps": ' + str(steps) + "," + "\n"
    # Open the "adaptive_1.8.conf" file for writing
    with open("adaptive_1.8.conf", "w") as f:
        f.writelines(lines)


# Main function
def main():
    # Get the number of iterations
    iterations = count_files()
    # Loop for each iteration
    for i in range(iterations):
        # Call the replace_files function
        replace_files(i)
        # Call the update_iterations function
        update_iterations(i+1)
        # When the iteration is 0 or the last iteration, change the steps to 50
        if i == 0 or i == iterations-1:
            update_steps(50)
        # When the iteration is 1 change the steps to 10
        elif i == 1:
            update_steps(10)
        # Execute the simulation
        adaptiveSampling.main("adaptive_1.8.conf")
    # Call results.main()
    #results.main()
    #work.main()
    #bar_analysis_2.main()

# Execute the script if it is the main program
if __name__ == "__main__":
    main()
