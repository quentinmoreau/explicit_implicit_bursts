import pickle
import sys
import numpy as np
import pandas as pd
from copy import copy
import pickle
import json
import time
from pathlib import Path
from extra.util import get_directories, get_files, make_directory, check_many
from mne import read_epochs
from extra.superlet import superlet, scale_from_period
from extra.burst_detection import extract_bursts, harmonics_detection
from datetime import datetime
import warnings


def extract_from_df(value, comp_column, targ_column):
    val_ix = comp_column.loc[comp_column == value].index.values.flatten()
    output = targ_column.iloc[val_ix].to_numpy().flatten()
    if output.shape[0] == 0:
        return np.nan
    else:
        return output

warnings.filterwarnings("ignore")
# parsing command line arguments
try:
    index = int(sys.argv[1])
except:
    print("incorrect arguments")
    sys.exit()

try:
    json_file = sys.argv[2]
    print("USING:", json_file)
except:
    json_file = "settings.json"
    print("USING:", json_file)

start_time = time.time()

# opening a json file
with open(json_file) as pipeline_file:
    parameters = json.load(pipeline_file)

max_freq = 120
foi = np.linspace(1, max_freq, 200)
scales = scale_from_period(1/foi)

search_range = np.where((foi >= 10) & (foi <= 33))[0]
beta_range = np.where((foi >= 13) & (foi <= 30))[0]
beta_lims = [13, 30]

path = Path(parameters["dataset_path"])
sfreq = parameters["downsample_dataset"]
hi_pass = parameters["hi_pass_filter"]
sub_path = path.joinpath("data")

output_path = Path("/home/common/mszul/explicit_implicit_beta/derivatives/PCA_results/beta_output/")

der_path = make_directory(path, "derivatives", check=True)
proc_path = make_directory(der_path, "processed", check=True)

epo_files = get_files(proc_path, "*-motor-epo.fif")
beh_files = get_files(proc_path, "*-beh.csv")
match_files = get_files(proc_path, "*-beh-match.json")

epo_file = epo_files[index]
subject = epo_file.parts[-2]
run = epo_file.stem.split("-")[-3]
beh_file = [i for i in beh_files if check_many([subject, run], i.stem, func="all")][0]
match_file = [i for i in match_files if subject in i.stem][0]
with open(match_file) as pipeline_file:
    beh_match = json.load(pipeline_file)
indices = beh_match[epo_file.name]

beta_output_file = output_path.joinpath(f"{epo_file.stem}-beta_power.pickle")

if beta_output_file.exists():
    print(f"{epo_file.stem} files are already processed, bye!")
    sys.exit()
else:
    pass

epochs = read_epochs(epo_file, verbose=False)
epochs = epochs.pick_types(meg=True, ref_meg=False, eeg=False, misc=False)
sens = ["MLC1", "MLC25", "MLC32", "MLC42", "MLC54", "MLC55", "MLC63"]
channels_used = [i for i in epochs.ch_names if check_many(["MLC"], i, func="any")]
channels_used = [i for i in channels_used if not check_many(sens, i, func="any")]
cu_ix = {cu: epochs.ch_names.index(cu) for cu in channels_used}


times = epochs.times
sfreq = epochs.info["sfreq"]
epoch_array = epochs.get_data()
trials, _, _ = epoch_array.shape

misc_dict = {}
for ixx, channel in enumerate(channels_used):
    print(f"{epo_file.stem} channel {ixx+1}/{len(channels_used)} {datetime.now()}")
    channel_ix = cu_ix[channel]
    tf_channel = [
        np.single(np.abs(superlet(
            epoch_array[trial, channel_ix], sfreq, 
            scales, 40, c_1=2, adaptive=True
        ))) for trial in range(trials)
    ]
    
    tf_channel = np.array([np.mean(i[beta_range,:], axis=0) for i in tf_channel])
    
    misc_dict[channel] = tf_channel

misc_dict["time"] = times
misc_dict["epo_to_beh_indices"] = indices
misc_dict["matched_beh_file"] = beh_file

with open(beta_output_file, "wb") as output_file:
    pickle.dump(misc_dict, output_file)

end_time = time.time()
print(f"{epo_file.stem} {(time.time() - start_time)/60} DONE")
