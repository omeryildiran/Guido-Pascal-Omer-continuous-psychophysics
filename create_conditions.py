# Create conditions csv file for gaussian_noise_v2.py
# conditions are the blob width withe different arcmin
import numpy as np
import pandas as pd
import random
def condition_creater():
    blob_widths=[11,13,17,21,25,29]
# repeat each condition 5 times
    repeats=5
# create a list of conditions
    trialN_per_condition=20
    trialN_per_block=5
    total_trialN=trialN_per_condition*len(blob_widths)
    conditions=[]
    for i in range(trialN_per_condition//trialN_per_block):
        random.shuffle(blob_widths)
        blob_widths_repeated=np.repeat(blob_widths,repeats)
        conditions.append(blob_widths_repeated)
    conditions=np.array(conditions).flatten()
    return conditions

#conditions = condition_creater()

# # create a dataframe
# conditions_df=pd.DataFrame(conditions, columns=['blob_width'])
# # add index
# conditions_df['index']=conditions_df.index
# # add block number
# conditions_df['block']=conditions_df['index']//trialN_per_block
# # add trial number
# conditions_df['trial']=conditions_df['index']%trialN_per_block
# # add trial number
# conditions_df['trialN']=conditions_df['index']+1
# # add total trial number
# conditions_df['total_trialN']=total_trialN

# # create a csv file
# conditions_df.to_csv('conditions.csv', index=False)