"""
Iterate over analyses, and create summary statistics over replicates for each.
"""

import pyabc
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cloudpickle as pickle


# iterate over models and analysis settings

files = [f for f in os.listdir('results') if f.endswith('.db')]

print("Collect files ...")

files_dict = {}
for f in files:
    splits = f.split('__')
    files_dict.setdefault(splits[0][3:], {}).setdefault(splits[1], []).append(f)

if os.path.isfile("total_samples.dat"):
    with open("total_samples.dat", 'rb') as f:
        total_samples_dict = pickle.load(f)
else:
    # for each, get total samples number
    print("Collect sample numbers ...")
    total_samples_dict = {}
    # samples_dict = {}
    for model_id, model_dict in files_dict.items():
        total_samples_dict[model_id] = {}
        # samples_dict[model_id] = {}
        for analysis_id, analysis_list in model_dict.items():
            total_samples_dict[model_id][analysis_id] = []
            # samples_dict[model_id[analysis_id] = []
            for f in analysis_list:
                h = pyabc.History("sqlite:///results/" + f)
                h.id = 1  # to be sure
                samples = h.get_all_populations()['samples']
                # sample_dict[model_id][analysis_id].append(samples)
                total_samples_dict[model_id][analysis_id].append(sum(samples))
    with open("total_samples.dat", 'wb') as f:
        pickle.dump(total_samples_dict, f)

# visualize total sample numbers
print("Plot total sample numbers ...")
for model_id, model_dict in total_samples_dict.items():
    fig, ax = plt.subplots()
    means = []
    stds = []
    labels = []
    for analysis_id, analysis_list in model_dict.items():
        labels.append(analysis_id)
        means.append(np.mean(analysis_list))
        stds.append(np.std(analysis_list))
    sorted_indices = sorted(range(len(labels)), key=labels.__getitem__)
    labels = [labels[i] for i in sorted_indices]
    means = [means[i] for i in sorted_indices]
    stds = [stds[i] for i in sorted_indices]
    ax.bar(range(len(labels)), height=np.log10(means))
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_title(f"Total samples for model {model_id}")
    ax.set_xlabel("Acceptor algorithm")
    ax.set_ylabel("$\log_{10}$(Samples)")
    fig.tight_layout()
    plt.savefig(f"total_samples__{model_id}.png")

# visualize sample numbers
print("Plot sample numbers ...")
for model_id, model_dict in files_dict.items():
    for i_rep in range(10):
        histories = []
        labels = []
        for analysis_id, analysis_list in model_dict.items():
            histories.append(pyabc.History("sqlite:///results/" + analysis_list[i_rep]))
            labels.append(analysis_id)
        pyabc.visualization.plot_sample_numbers(histories, labels, rotation=90)
        plt.savefig(f"samples__{model_id}__{i_rep}.png")
