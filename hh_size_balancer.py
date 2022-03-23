# this program takes the first run HH size results and tries to shift household distribution upwards.
# So the average household size for 7+ could be reasonable.

# inputs:
# popsim control file with hh sizes and POPBASE, such as SEMCOG_2019_control_totals_blkgrp.csv
# sythesized summary from previous run, such as final_summary_BLKGRP.csv

# output:
# popsim control with updated HH size: SEMCOG_2019_control_totals_blkgrp_hhsize_adj.csv

# %%
import pandas as pd
import numpy as np
import os


# %%
def rebalance(x0, tar, w=None):
    # adjust household size distribution based on population targets
    # x0: # of hhs by hh size (list/vector)
    # tar: target population by geo unit(block group, tract, etc)
    # weight: for HHsize 7+, the average size in this block group (based on first synthesis run)
    x = x0.copy()
    hh = x.sum()
    rx = np.array(range(len(x0)))
    if w is None:
        w = rx + 1
    while x.dot(w) < tar and x.max() < hh:
        # print tar, x.dot(w), x
        i = np.random.choice(rx[:-1], p=(1.0 * x[:-1] / x[:-1].sum()))
        j = np.random.choice(rx[i:], p=(1.0 * x[i:] / x[i:].sum()))
        x[i] -= 1
        x[j] += 1
    return x


def update_hhsize_control(df_ctrl, df_output):
    cols = [[f'HHPERSONS{x}', f'hh_persons_{x}_control', f'hh_persons_{x}_result']
            for x in range(1, 8)]
    syn_ctrls, out_ctrls, out_reslts = [list(c) for c in zip(*cols)]

    # calculate "actual"/synthesized household size for 7+ person HH from first run
    df_output['weight'] = (df_output['persons_num_result']-(df_output[out_reslts[:-1]]
                                                            * [1, 2, 3, 4, 5, 6]).sum(axis=1))/df_output[out_reslts[-1]]
    df_output['weight'].fillna(0, inplace=True)
    df_output.drop('geography', axis=1, inplace=True)
    # before = (df_output[out_ctrls] * [1, 2, 3, 4, 5, 6, 7]).sum().sum()

    # rebalance HH size and create new HH size target DF
    new_targets = []
    for id, row in df_output.iterrows():
        ntarget = rebalance(row[out_ctrls].values, row['persons_num_control'],
                            [1, 2, 3, 4, 5, 6, max(row['weight'], 7)])
    # print(row['persons_num_control'], ntarget)
        new_targets.append([id] + list(ntarget))

    df_new_target = pd.DataFrame(
        new_targets, columns=['BLKGRPID'] + syn_ctrls).set_index('BLKGRPID').astype(int)
    # after = (df_new_target[ctrl_cols]*[1, 2, 3, 4, 5, 6, 7]).sum().sum()

    return df_new_target, syn_ctrls


# %%
if __name__ == "__main__":

    ctrl_file = 'data/SEMCOG_2019_control_totals_blkgrp.csv'
    output_file = 'output/run6_pass1/final_summary_BLKGRP.csv'

    df_ctrl = pd.read_csv(ctrl_file, index_col='BLKGRPID')
    df_output = pd.read_csv(output_file, index_col='id')

    df_new_target, ctrl_cols = update_hhsize_control(df_ctrl, df_output)
    df_ctrl[ctrl_cols] = df_new_target[ctrl_cols]

    # back up control file from previous run, change its name to xxx_run##.csv
    i = 0
    ctrl_backup = ctrl_file.replace('.csv', '_run%s.csv')
    while os.path.exists(ctrl_backup % i):
        i += 1
    os.rename(ctrl_file, ctrl_backup % i)


    # save 2 copies of new control files, one has original name, the other as a backup "hhsize_adj"
    df_ctrl.to_csv(ctrl_file)
    df_ctrl.to_csv(ctrl_backup.replace('.csv', '_hhsize_adj.csv') % i)

    print(f'original control file {ctrl_file} has been renamed to ', ctrl_backup % i)
    print(f'new control file "{ctrl_file}" were produced and replaced the original control')
   


# %%
