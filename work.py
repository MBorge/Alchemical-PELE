import os
#from simulation import replace_files
#from simulation import count_files
import subprocess
import shutil

def replace_files(num):
    # Define the source file paths using string formatting
    hybz_source_file = f"alchemical_files/hybz_{num}"
    ligandParams_source_file = f"alchemical_files/ligandParams_{num}.txt"
    HYB_rotamers_source_file = f"alchemical_files/HYB.rot.assign{num}"

    # Use shutil to copy the source files to the destination directories
    shutil.copy(hybz_source_file, "DataLocal/Templates/OpenFF/Parsley/hybz")
    shutil.copy(ligandParams_source_file, "DataLocal/OBC/ligandParams.txt")
    shutil.copy(HYB_rotamers_source_file,"DataLocal/LigandRotamerLibs/HYB.rot.assign")


def check__num_lambdas():
    os.chdir("alchemical_files")
    lambdas = len([name for name in os.listdir(".") if os.path.isfile(name) and name.startswith("hybz")])
    os.chdir("..")
    return lambdas


def check_num_trajectories():
    os.chdir("output/0/")
    trajectories = len([name for name in os.listdir(".") if os.path.isfile(name) and name.startswith("trajectory")])
    os.chdir("../..")
    return trajectories


def prepare_pdb():
    os.system("cp -r ../../alchemical_files/ .")
    os.system("cp ../../../../pele2.conf .")
    os.system("cp -r ../../../../DataLocal .")
    
    files = os.listdir(".")
    model = [file for file in files if file.startswith("model")]
    os.system("mv " + model[0] + " complex_1.pdb")


def recalculations(l, lambdas):
    lambdas_parameters = []
    if l == 0:
        lambdas_parameters = [l, l+1]
    elif l == lambdas-1:
        lambdas_parameters = [l-1, l]   
    else:
        lambdas_parameters = [l-1, l, l+1]
    for i in lambdas_parameters:
        replace_files(i)
        command = "/eb/x86_64/software/PELE/1.8b3/bin/PELE_serial pele2.conf"
        subprocess.call(command, shell=True)
        os.system("mv output output_" + str(i))


def model_files_creator(l, lambdas):
    files = os.listdir(".")
    report = [file for file in files if file.startswith("trajectory")]

    with open(report[0], "r") as file:
        lines = file.read()
        models = lines.split("MODEL")
        models = models[1:] 

    for i, model in enumerate(models, start=1):
        folder = f"model_{i}"
        os.makedirs(folder, exist_ok=True)
        out = os.path.join(folder, f"modelo_{i}.pdb")
        with open(out, 'w') as file:
            file.write("MODEL" + model)
        os.chdir(folder)
        prepare_pdb()
        recalculations(l, lambdas)
        os.chdir("..")

    print(f"{len(models)} PDB files in individual folders.")


def main():
    lambdas = check__num_lambdas()
    trajectories = check_num_trajectories()
    
    work_directory = "alchemical_work"
    if not os.path.exists(work_directory):
        os.mkdir(work_directory)
    os.chdir(work_directory)

    for l in range(lambdas):
        os.mkdir(str(l))
        os.chdir(str(l))

        os.mkdir("alchemical_files")
        os.chdir("alchemical_files")

        os.system("cp ../../../alchemical_files/hybz_" + str(l) + " .")
        os.system("cp ../../../alchemical_files/HYB.rot.assign" + str(l) + " .")
        os.system("cp ../../../alchemical_files/ligandParams_" + str(l) + ".txt" + " .")

        if l != lambdas-1:
            os.system("cp ../../../alchemical_files/ligandParams_" + str(l+1) + ".txt" + " .")
            os.system("cp ../../../alchemical_files/HYB.rot.assign" + str(l+1) + " .")
            os.system("cp ../../../alchemical_files/hybz_" + str(l+1) + " .")

        if l != 0:
            os.system("cp ../../../alchemical_files/ligandParams_" + str(l-1) + ".txt" + " .")
            os.system("cp ../../../alchemical_files/HYB.rot.assign" + str(l-1) + " .")
            os.system("cp ../../../alchemical_files/hybz_" + str(l-1) + " .")
        
        os.chdir("..")

        for t in range(trajectories):
            os.mkdir("trajectory_" + str(t+1))
            os.chdir("trajectory_" + str(t+1))            
            os.system("cp ../../../output/" + str(l) + "/trajectory_" + str(t+1) + ".pdb ." )
            os.system("cp ../../../output/" + str(l) + "/report_" + str(t+1) +  " .")

            model_files_creator(l, lambdas)

            os.chdir("..")
        os.chdir("..")
    os.chdir("..")

if __name__ == "__main__":
    main()
