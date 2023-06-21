import os
import pymbar
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def check_steps(l,t,lambdas):
    times = []

    with open("report_{}".format(t+1), "r") as f:
        lines = f.readlines()[1:]

    if l == 0 or l == lambdas-1:
        steps = 50
        accepted_steps = []
        for line in lines:
            accepted_steps.append(int(line.split()[1]))
        
        times = [accepted_steps[i+1]- accepted_steps[i] for i in range(len(accepted_steps)-1)]
        last = steps - accepted_steps[-1]
        times.append(last)

    else:
        steps = 10
        accepted_steps = []
        for line in lines:
            accepted_steps.append(int(line.split()[1]))
        
        times = [accepted_steps[i+1]- accepted_steps[i] for i in range(len(accepted_steps)-1)]
        last = steps - accepted_steps[-1]
        times.append(last)
        
    return times



def plot_histograms(csv_file):
    df = pd.read_csv(csv_file)
    for column in df:
        plt.hist(df[column], bins=50, edgecolor="black", density=True)
        plt.savefig("histogram_{}.png".format(column))
        plt.clf()
    


def obtain_energies(lambdas, trajectories):
    for l in range(lambdas):
        os.chdir("alchemical_work/{}".format(l))
        forward_global = []
        reverse_global = []
        energies_lambda_before_global = []
        energies_lambda_next_global = []
        energies_lambda_global = []

        for t in range(trajectories):
            os.chdir("trajectory_{}".format(t+1))
            models = len([name for name in os.listdir(".") if name.startswith("model")])
            forward_traj = []
            reverse_traj = []
            steps_for_conf = check_steps(l,t,lambdas)
            energies_lambda_before = []
            energies_lambda_next = []
            energies_lambda = []
            for m in range(models):
                if l == 0 and m == 0:
                    continue

                energies = []
                os.chdir("model_{}".format(m+1))

                if l == 0:
                    for i in range(2):
                        os.chdir("output_{}".format(i))
                        with open("report", "r") as f:
                            line = f.readlines()[1]
                            energy = float(line.split()[4])
                            energies.append(energy)
                        os.chdir("..")

                    number_of_steps = steps_for_conf[m]
                    repeated_number = list(itertools.repeat(energies[1] - energies[0], number_of_steps))
                    repeated_energy_lambda = list(itertools.repeat(energies[0], number_of_steps))
                    repeated_energy_lambda_next = list(itertools.repeat(energies[1], number_of_steps))

                    forward_global.extend(repeated_number)
                    forward_traj.extend(repeated_number)
                    energies_lambda.extend(repeated_energy_lambda)
                    energies_lambda_next.extend(repeated_energy_lambda_next)
                    energies_lambda_global.extend(repeated_energy_lambda)
                    energies_lambda_next_global.extend(repeated_energy_lambda_next)
                    #print("Forward: {}".format(forward_traj))
                        
                elif l == lambdas-1:
                    for i in [l-1,l]:
                        os.chdir("output_{}".format(i))
                        with open("report", "r") as f:
                            line = f.readlines()[1]
                            energy = float(line.split()[4])
                            energies.append(energy)
                        os.chdir("..")

                    number_of_steps = steps_for_conf[m]
                    repeated_number = list(itertools.repeat(energies[0] - energies[1], number_of_steps))
                    repeated_energy_lambda_before = list(itertools.repeat(energies[0], number_of_steps))
                    repeated_energy_lambda = list(itertools.repeat(energies[1], number_of_steps))

                    reverse_global.extend(repeated_number)
                    reverse_traj.extend(repeated_number)
                    energies_lambda_before.extend(repeated_energy_lambda_before)
                    energies_lambda.extend(repeated_energy_lambda)
                    energies_lambda_before_global.extend(repeated_energy_lambda_before)
                    energies_lambda_global.extend(repeated_energy_lambda)
                    #print("Reverse: {}".format(reverse_global))
            
                else:
                    for i in [l-1,l,l+1]:
                        os.chdir("output_{}".format(i))
                        with open("report", "r") as f:
                            line = f.readlines()[1]
                            energy = float(line.split()[4])
                            energies.append(energy)
                        os.chdir("..")

                    number_of_steps = steps_for_conf[m]
                    repeated_number_forward = list(itertools.repeat(energies[2] - energies[1], number_of_steps))
                    repeated_number_reverse = list(itertools.repeat(energies[0] - energies[1], number_of_steps))
                    repeated_energy_lambda_before = list(itertools.repeat(energies[0], number_of_steps))
                    repeated_energy_lambda = list(itertools.repeat(energies[1], number_of_steps))
                    repeated_energy_lambda_next = list(itertools.repeat(energies[2], number_of_steps))
                    
                    reverse_global.extend(repeated_number_reverse)
                    reverse_traj.extend(repeated_number_reverse)
                    forward_global.extend(repeated_number_forward)
                    forward_traj.extend(repeated_number_forward)
                    energies_lambda_before.extend(repeated_energy_lambda_before)
                    energies_lambda.extend(repeated_energy_lambda)
                    energies_lambda_next.extend(repeated_energy_lambda_next)
                    energies_lambda_before_global.extend(repeated_energy_lambda_before)
                    energies_lambda_global.extend(repeated_energy_lambda)
                    energies_lambda_next_global.extend(repeated_energy_lambda_next)
                    #print("Forward: {}".format(forward_global))
                    #print("Reverse: {}".format(reverse_global))

                #print("Lambda: {}, Trajectory: {}, Model: {}".format(l,t+1,m+1))
                #print(energies)
                os.chdir("..")

            if l == 0:
                data_t = {"forward": forward_traj}
                df_t = pd.DataFrame(data_t)
                df_t.to_csv("diff_traj.csv", index=False)
                #plot_histograms("diff_traj.csv")

                data_energy = {"lambda": energies_lambda, "lambda_next": energies_lambda_next}
                df_energy = pd.DataFrame(data_energy)
                df_energy.to_csv("energy_traj.csv", index=False)
                #plot_histograms("energy_traj.csv")

            elif l == lambdas-1:
                data_t = {"reverse": reverse_traj}
                df_t = pd.DataFrame(data_t)
                df_t.to_csv("diff_traj.csv", index=False)
                #plot_histograms("diff_traj.csv")

                data_energy = {"lambda_before": energies_lambda_before, "lambda": energies_lambda}
                df_energy = pd.DataFrame(data_energy)
                df_energy.to_csv("energy_traj.csv", index=False)
                #plot_histograms("energy_traj.csv")

            else:
                data_t = {"forward": forward_traj, "reverse": reverse_traj}
                df_t = pd.DataFrame(data_t)
                df_t.to_csv("diff_traj.csv", index=False)
                #plot_histograms("diff_traj.csv")

                data_energy = {"lambda_before": energies_lambda_before, "lambda": energies_lambda, "lambda_next": energies_lambda_next}
                df_energy = pd.DataFrame(data_energy)
                df_energy.to_csv("energy_traj.csv", index=False)
                #plot_histograms("energy_traj.csv")

            os.chdir("..")

        if l == 0:
            data = {"forward": forward_global}
            df = pd.DataFrame(data)
            df.to_csv("diff_lambda.csv", index=False)
            plot_histograms("diff_lambda.csv")

            data_t_energy = {"lambda": energies_lambda_global, "lambda_next": energies_lambda_next_global}
            df_t_energy = pd.DataFrame(data_t_energy)
            df_t_energy.to_csv("energy_lambda.csv", index=False)
            plot_histograms("energy_lambda.csv")
        
        elif l == lambdas-1:
            data = {"reverse": reverse_global}
            df = pd.DataFrame(data)
            df.to_csv("diff_lambda.csv", index=False)
            plot_histograms("diff_lambda.csv")

            data_t_energy = {"lambda_before": energies_lambda_before_global, "lambda": energies_lambda_global}
            df_t_energy = pd.DataFrame(data_t_energy)
            df_t_energy.to_csv("energy_lambda.csv", index=False)
            plot_histograms("energy_lambda.csv")

        else:
            data = {"forward": forward_global, "reverse": reverse_global}
            df = pd.DataFrame(data)
            df.to_csv("diff_lambda.csv", index=False)
            plot_histograms("diff_lambda.csv")

            data_t_energy = {"lambda_before": energies_lambda_before_global, "lambda": energies_lambda_global, "lambda_next": energies_lambda_next_global}
            df_t_energy = pd.DataFrame(data_t_energy)
            df_t_energy.to_csv("energy_lambda.csv", index=False)
            plot_histograms("energy_lambda.csv")

        os.chdir("../..")



def bar_calculation(lambdas):
    results = []
    for l in range(lambdas-1):
        forward = []
        reverse = []
        os.chdir("alchemical_work/{}".format(l))
        current_df = pd.read_csv("diff_lambda.csv")
        forward.append(current_df["forward"].tolist())
        next_folder = str(l+1)
        os.chdir("../{}".format(next_folder))
        next_df = pd.read_csv("diff_lambda.csv")
        reverse.append(next_df["reverse"].tolist())
        os.chdir("..")
        forward = np.array(forward)
        reverse = np.array(reverse)

        # Do boostrap with the vectors




        bar = pymbar.other_estimators.bar(forward, reverse)
        #print(bar)
        results.append(bar['Delta_f'])
        with open("bar_results.txt", "a") as f:
            f.write("Lambda: {} - {}, BAR: {}\n".format(l, l+1, bar))
        os.chdir("..")

    rbfe = np.sum(results)
    #print("RBFE: {}".format(rbfe))   
    with open("rbfe.txt", "a") as f:
        f.write("RBFE: {}\n".format(rbfe))

    transitions = range(len(results))
    plt.scatter(transitions, results)
    plt.savefig("transition.png")
    plt.clf()



def main():
    os.chdir("alchemical_work")
    lambdas =  len(os.listdir())
    os.chdir("0")
    trajectories = len(os.listdir())-1
    os.chdir("../..")
    obtain_energies(lambdas, trajectories)
    bar_calculation(lambdas)



if __name__ == "__main__":
    main()