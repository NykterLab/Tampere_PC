#!/usr/bin/env python3

from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from bpnet.modisco.files import ModiscoFile, ModiscoFileGroup
import os
import sys
import argparse as ap

parser = ap.ArgumentParser(description='Run TF-MoDISCo for a given BPNET model')
parser.add_argument('--model_dir', dest='model_dir', help='BPNET model directory')
parser.add_argument('--task', dest='task', help='name of the BPNET task whose predictions are being analyzed; typically the name of a TF')
args = parser.parse_args()
model_dir = args.model_dir
task = [args.task]

model_dir_path = Path(model_dir)
modisco_dir = model_dir_path/'modisco'

## MultipleModiscoResult is a convenience wrapper around ModiscoResult
mf = ModiscoFileGroup({t: ModiscoFile(modisco_dir / t / 'modisco.h5') for t in task})

## create directories for plots: each dir corresponds to a different task/metacluster/
task_metacluster_dirs = []
for i in mf.pattern_names():
    dir_name_items = i.split('/')
    dir_name = '/'.join(dir_name_items[:len(dir_name_items) - 1])
    task_metacluster_dirs.append(dir_name)
task_metacluster_dirs = list(set(task_metacluster_dirs))
for t in task_metacluster_dirs:
    Path(f"{modisco_dir}/modisco_plots/{t}").mkdir(parents=True, exist_ok=True)

## Plot all the patterns
fig_names = []
for p in mf.patterns():
    fig = p.plot("seq_ic", title=f"{p.name} ({mf.n_seqlets(p.name)})")
    plt.ylim([0, 2])
    sns.despine(top=True, bottom=True, right=True)
    fig.savefig(f"{modisco_dir}/modisco_plots/{p.name}.png", bbox_inches = "tight", dpi=300)
    fig_names.append(f"{modisco_dir}/modisco_plots/{p.name}.png")

print(fig_names)


## Join the plots
fig, ax = plt.subplots(len(mf.patterns()), frameon=False)
for p in range(len(fig_names)):
    name = fig_names[p]
    img=mpimg.imread(name)
    ax[p].frameon=False
    ax[p].set_frame_on(False)
    ax[p].xaxis.set_visible(False)
    ax[p].yaxis.set_visible(False)
    ax[p].imshow(img)
fig.savefig(f"{modisco_dir}/modisco_plots/all_patterns.png", bbox_inches = "tight", dpi=1200)
