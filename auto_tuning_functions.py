# behavior figures

# looking at bulk tuning curves in one plot for ACx and IC

import pandas as pd
import numpy as np

import math
from scipy.optimize import curve_fit
#plots
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.signal import find_peaks

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from array import *
import sys
sys.path.append('/home/dawg/coding_stuff/spike_distance_tune/PySpike') # your path to Pyspike folder
sys.path.append('/home/dawg/coding_stuff/spike_distance_tune/PySpike/test') # path to Pyspike test of you need to use testing plots

import pyspike as spk

#freq_cols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] # need to set this to match your own data future version will be more general
                                                                                                        

freq_cols = [0] # need to set this to match your own data future version will be more general


                         

#  (`-')                <-. (`-')_  _     <-. (`-')_ (`-')  _ (`-')                   (`-')  
#  ( OO).->       .->      \( OO) )(_)       \( OO) )(OO ).-/ ( OO).->       .->   <-.(OO )  
#  /    '._  ,--.(,--.  ,--./ ,--/ ,-(`-'),--./ ,--/ / ,---.  /    '._  (`-')----. ,------,) 
#  |'--...__)|  | |(`-')|   \ |  | | ( OO)|   \ |  | | \ /`.\ |'--...__)( OO).-.  '|   /`. ' 
#  `--.  .--'|  | |(OO )|  . '|  |)|  |  )|  . '|  |)'-'|_.' |`--.  .--'( _) | |  ||  |_.' | 
#     |  |   |  | | |  \|  |\    |(|  |_/ |  |\    |(|  .-.  |   |  |    \|  |)|  ||  .   .' 
#     |  |   \  '-'(_ .'|  | \   | |  |'->|  | \   | |  | |  |   |  |     '  '-'  '|  |\  \  
#     `--'    `-----'   `--'  `--' `--'   `--'  `--' `--' `--'   `--'      `-----' `--' '--' 

# Written by Walker Gauthier 2024 UIUC @ Auerbach Lab 
# email: davidwg3@illinois.edu                                                                                       
                                                                                       
                                                                                    
def is_significantly_activated(current_cell,sanity,filename):
    #print(current_cell)
    current_cell = np.concatenate(current_cell).ravel()#.tolist()

    #print(current_cell)
    #prestim_spks = count_range_in_list(current_cell,-0.05,0)
    bin1count = 0
    bin2count = 0
    bin3count = 0
    bin4count = 0
    bin5count = 0
    bin6count = 0
    bin7count = 0
    bin8count = 0
    bin9count = 0
    bin10count = 0
    

    for spike in current_cell:
        if(spike > -0.05 and spike < -0.045):
            bin1count = bin1count + 1
        if(spike > -0.045 and spike < -0.04):
            bin2count = bin2count + 1
        if(spike > -0.040 and spike < -0.035):
            bin3count = bin3count + 1
        if(spike > -0.035 and spike < -0.03):
            bin4count = bin4count + 1
        if(spike > -0.03 and spike < -0.025):
            bin5count = bin5count + 1
        if(spike > -0.025 and spike < -0.02):
            bin6count = bin6count + 1
        if(spike > -0.02 and spike < -0.015):
            bin7count = bin7count + 1
        if(spike > -0.015 and spike < -0.01):
            bin8count = bin8count + 1
        if(spike > -0.01 and spike < -0.05):
            bin9count = bin9count + 1
        if(spike > -0.05 and spike < -0.00009):
            bin10count = bin10count + 1
    spk_ls = [bin1count,bin2count,bin3count,bin4count,bin5count
                          ,bin6count,bin7count,bin8count,bin9count]

    prestim_sd = np.array(spk_ls)
    prestim_sd = np.std(prestim_sd)
    
    prestim_avg = np.mean(prestim_sd)
    
    #print(spk_ls)
    #print(prestim_sd)
    
    
    
    prestim_sd4 = prestim_sd*4 # multiply by 4
    #print(prestim_sd4)
    
    
    # Define the bin edges with a set width
    bin_width = 0.005
    bin_edges = np.arange(min(current_cell), max(current_cell) + bin_width, bin_width)

    # Compute the histogram with the specified bin edges
    hist, _ = np.histogram(current_cell, bins=bin_edges)
    
    #print(hist)
    #print(_)

    
    
    # Find the index of the maximum value in the histogram
    peak_index = np.argmax(hist)

    # Get the peak value
    peak_latency = bin_edges[peak_index] # calculate post stim peak FR (NP histogram peak from 0 to 50 ms should work)
    
    peak_fr = hist[peak_index]
    bl_peak_fr = peak_fr - prestim_avg
    #peak_fr = peak_count
    #print(peak_latency)
    #print(peak_fr)
    time_0_ind = np.where(_ >= 0)[0][0]
    #print(time_0_ind)
    
    
    if peak_latency > 0 and peak_fr > prestim_sd4:
        # if peak firing rate exceeds that, then it is sound evoked
        evoked_Status = True
        thresh_ind = np.where(hist > prestim_sd4)[0][0]
        
        
        if thresh_ind >= time_0_ind:
            first_bin_exceeding_threshold = _[thresh_ind]
        else:
            first_bin_exceeding_threshold = -10
            
        #print('Thresh index: ',thresh_ind)
        #print('SD4 time: ',first_bin_exceeding_threshold)
            
        #print('Evoked')
        
    else:
        evoked_Status = False
        first_bin_exceeding_threshold = -10
        #print('Not evoked')

        
        
    # Find the time range of bins following the peak bin until the count goes below the standard deviation
    bins_after_peak = []
    for i in range(peak_index, len(hist)):
        if hist[i] >= prestim_sd4:
            bins_after_peak.append(bin_edges[i])
            
    #print('Bins after peak: ', bins_after_peak)
    if not bins_after_peak:
        #print('Empty Duration')
        duration = 0.005
        last_bin = 0.02
    
    else:   
        duration = (max(bins_after_peak) - min(bins_after_peak))
        last_bin = max(bins_after_peak)
    
    duration_rounded = round(duration, 8)

    start_bin = np.searchsorted(bin_edges,0.01, side = 'left') # get total spike count from 10 to 40ms (the reason we are doing this is becuase non-evoked cells would not have an evoked portion)
    end_bin = np.searchsorted(bin_edges,0.04,side='right') # can use this same function but use the last bin and peak latency times for the second argument (np.searchsorted(bin_edges,peak_latency,side='left))
    sum_counts = hist[start_bin:end_bin].sum()
    # still need spike count during evoked portion for evoked cells to compare genotypes.

    #print('Response Duration: ', duration_rounded)
    
    # sanity plots if needed
    if sanity == 'yes':
        plt.figure(figsize=(8, 6))
        plt.hist(current_cell, bins=bin_edges, color='blue', edgecolor='black', alpha=0.7)
        plt.axvline(x=peak_latency, color='red', linestyle='--', label=f'Peak Value: {peak_latency:.2f}')
        plt.axhline(y=prestim_sd4,color='black',linestyle='--')
        plt.axvspan(0, 0.1, facecolor='g',alpha=.1,zorder=-100) #new line
        plt.title(filename)
        
        plt.show()
    
    evoked_ls = [evoked_Status, peak_latency, peak_fr, duration_rounded,bl_peak_fr,first_bin_exceeding_threshold,last_bin,sum_counts]
    
    return evoked_ls


def make_evoked_df(current_df,freq_ls,db_ls,sanity):

    

    rows = len(db_ls)
    cols = len(freq_ls)

    sound_evoked_array = []
    non_evoked_array = []

    big_df = current_df

    new_df = pd.DataFrame(columns = big_df.columns.values) # placeholder dataframe to append to
                        # not efficient at all, but I don't want to wrestle with replacing by index either


    for i in big_df['file'].unique(): # loop over each tuning curve
        current_file = big_df.loc[big_df['file'] == i]
        
        for chan in current_file['channel'].unique(): 
            test = current_file.loc[current_file['channel'] == chan]
            
            current_geno = test['Genotype'][0] # grab the current genotype
            spiking_only_df = test[range(cols)] # grab only the columns that correspond to freqs
            #print(chan,x,current_geno)
            #print(spiking_only_df)

            current_output_df = spiking_only_df.map(lambda x: is_significantly_activated(x,sanity,i))  
            
            #print(current_output_df)
            current_output_df['file'] = i
            current_output_df['channel'] = chan
            current_output_df['Genotype'] = test['Genotype']
            
            new_df = pd.concat([new_df,current_output_df],ignore_index=False)
        
    return new_df

def plot_evoke_status_crude(evoked_df):

    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            #print(current_geno)
            print(test)
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'
            
            current_df = test[freq_cols].map(lambda x: x[0]) # get evoked True/False value
            current_df = current_df.astype(int) # converts true false to 1 0

            fr_df = test[freq_cols].map(lambda x: x[4]) # get baseline corrected firing rate

            fig1,ax1 = plt.subplots(2,1,figsize=(5,5))
            sns.heatmap(current_df,ax=ax1[0])
            sns.heatmap(fr_df,ax=ax1[1])
            ID_file = 'File: ' + str(i) +' Unit: ' + str(j)
            plt.title(ID_file)
            plt.tight_layout()
            plt.show()

    return current_df

def calc_dprime(evoked_ls,nonevoked_ls):

    #print(evoked_ls)

    #print(nonevoked_ls)

    mean_evoked = np.mean(evoked_ls)
    mean_nonevoked = np.mean(nonevoked_ls)

    sd_evoked = np.std(evoked_ls,ddof=1)
    sd_non_evoked = np.std(nonevoked_ls,ddof =1)

    avg_sd = (sd_evoked+sd_non_evoked)/2

    dprime = (mean_evoked-mean_nonevoked)/avg_sd


    return dprime

def make_dB_df(evoked_df):
    plot_df = pd.DataFrame(columns = ['latency','abs_peak_fr','rel_peak_fr','resp_duration','first_bin','last_bin','smoothed_data','spks_10_40ms','file','channel','genotype'])

    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            #print(current_geno)
            
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'
            
            
            for index, row in test.iterrows():
                
                current_df = test[freq_cols].apply(lambda x: list(x)[index]) # need to get only cells with True
                true_columns = current_df.columns[current_df.iloc[0] == True] # check where columns have a True string
                #print(true_columns)
                
                result_df = current_df[true_columns] # Select columns where 'True' is present in the first row
                
                
                # need to drop negative first_bin columns too

                true_columns = result_df.columns[result_df.iloc[5] >= 0] # check where columns have a True string
                result_df = result_df[true_columns]
                #print(result_df)
                
            
            
                latency = result_df.apply(lambda x: list(x)[1]) 
                #print(latency)
                lat_mean = latency.mean()
                # extract peak latency
                peak_fr = result_df.apply(lambda x: list(x)[2])
                
                peak_mean = peak_fr.mean()
                # exctract peak fr
                resp_dur = result_df.apply(lambda x: list(x)[3]) 
                dur_mean = resp_dur.mean()
                # extract response duration
                bl_peak_fr = result_df.apply(lambda x: list(x)[4])
                rel_peak_mean = bl_peak_fr.mean()
                # extract baseline relative peak fr
                # 3x3 smooth the array
                smoothed_data = median_filter(bl_peak_fr, size=3)
                smoothed_mean = smoothed_data.mean()
                
                first_sd4_bin = result_df.apply(lambda x: list(x)[5]) 
                first_bin = first_sd4_bin.mean()
                
                last_sd4_bin = result_df.apply(lambda x: list(x)[6]) 
                last_bin = last_sd4_bin.mean()

                spks_10_40ms = result_df.apply(lambda x: list(x)[7]) 
                spks_10_40ms_mean = spks_10_40ms.mean()


                plot_df = plot_df._append({'latency': lat_mean, 'abs_peak_fr': peak_mean,
                                                    'rel_peak_fr':rel_peak_mean,
                                                    'resp_duration': dur_mean,
                                                    'first_bin':first_bin,
                                                    'last_bin':last_bin,
                                                    'smoothed_data':smoothed_mean,
                                                    'spks_10_40ms':spks_10_40ms_mean,
                                                    'file': i,
                                                    'channel':j,
                                                    'genotype' : current_geno}, ignore_index=True)
    return plot_df           

def get_thresh(current_dbx1_array,db_ls,sanity):

# Smooth the line
    
    smoothed_y = savgol_filter(current_dbx1_array, window_length=5, polyorder=3)


    #thresh = 0.2*np.max(y)
    maxi=np.max(smoothed_y)
    addi = 0.2*maxi
    thresh_y = np.min(smoothed_y)+addi
    #print(maxi,thresh_y)

    thresh_rounded = round(thresh_y, 0)# round to nearestwhole number both thresh and smoothed Y
    smooth_rounded = np.round(smoothed_y, 0)
    #print(thresh_rounded,smooth_rounded)



     
    int_ind = np.where(smooth_rounded == thresh_rounded) # if there is an equal value then we are good
    

    #print(int_ind)
    try:
        thresh_ind = int_ind[0][0] # need to get the corresponding x value intensity
        
        #print('Not Empty')
    except:
        # if no equal value need to find the closest one to get the intensity index
        thresh_diff_array = smooth_rounded-thresh_rounded # need to get the closest smoothed value to our threshold fr
        smallest_diff = np.min(abs(thresh_diff_array))
        thresh_ind = np.where(abs(thresh_diff_array) == smallest_diff)
        #print('It was empty: ',smallest_diff,thresh_diff_array,thresh_ind)
        #print(thresh_ind[0][0])
        
        thresh_ind = thresh_ind[0][0]
    
    

        # add to get thresh for sanity plots
    thresh_db = db_ls[thresh_ind]
    '''
    if thresh_db == 0:
        # sometimes end up with 0 threshold due to the smoothing
        try:
            thresh_ind = int_ind[1][0] # so just take the next intercept
            thresh_db = db_ls[thresh_ind]
        except:
            thresh_ind = thresh_ind[1][0]
            thresh_db = db_ls[thresh_ind]
    '''


    if sanity == 'yes':
        plt.figure(figsize=(5,5))
        plt.plot(current_dbx1_array, label='Original Data')
        plt.plot(smoothed_y,label = 'Smooth')
        plt.axhline(y = thresh_y, color='k', label='Inflection',linestyle='--')
        sns.despine()
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.legend()
        plt.title('Smoothing with 5 point filter')
        plt.show()

    
    
    #print(thresh_db)

    return thresh_db

def add_thresh_col(dB_df,db_ls,sanity):

    
    thresh_df = pd.DataFrame(columns = ['threshold','file','channel','genotype'])
    intensity = db_ls

    for i in dB_df['file'].unique(): # need ti group by file and by channel
        current_file = dB_df.loc[dB_df['file'] == i]
        for j in current_file['channel'].unique():
            #print(j)
            current_unit = current_file.loc[current_file['channel'] == j]
            #print(current_unit)
            current_geno = current_unit['genotype'] # grab the current genotype
            #print(current_geno)
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT' # weird error were some geno labels were saved in an array

            curr_fr = current_unit['rel_peak_fr'] 
     
            thresh = get_thresh(curr_fr,intensity,sanity)

                # make a threshold dataframe columns being threshold, file, channel, genotype
            thresh_df = thresh_df._append({'threshold': thresh,
                                                    'file': i,
                                                    'channel':j,
                                                    'genotype' : current_geno}, ignore_index=True)
            
            

    return thresh_df
    



def make_freq_df(evoked_df,freq_ls):
    cf_df = pd.DataFrame(columns = ['file','channel','genotype'])
    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            #print(current_geno)
            #print(test)
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'

            current_df = test[freq_cols]
            #print(current_df)

            sum_all_freq = {}
            for series_name, series in current_df.items():
                current_freq = series

                                
                current_sum = np.sum(current_cell[4] for current_cell in current_freq) # sum the baseline corrected frs for this column
                current_freq = freq_ls[series_name]
                sum_all_freq[current_freq] = current_sum






                
            CF = max(sum_all_freq, key=sum_all_freq.get)
                

                # make a threshold dataframe columns being threshold, file, channel, genotype
            cf_df = cf_df._append({'CF': CF,
                                    'file': i,
                                    'channel':j,
                                    'genotype' : current_geno}, ignore_index=True)


                
    return cf_df




def get_dprime(evoked_df,cf_df,thresh_df,db_ls,freq_ls):
    # need to loop through each unit
    # need top get the threshold and BF from respective dataframes
    # the cells in the rows above threshold are all considered non-evoked
    # throw out threshold row for now (two peaks can mess this up)
    # take all sound evoked cells above the threshold row as evoked
    # do the d-prime formula
    # save into a dprime df of same format as the thresh and BF dataframes so we can merge if they want 
    dprime_df = pd.DataFrame(columns = ['file','channel','genotype'])
    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            #current_bf = cf_df.loc[(cf_df['file'] == i) & (cf_df['channel'] == j)]
            #current_bf = current_bf['CF']
            current_thresh = thresh_df.loc[(thresh_df['file'] == i) & (thresh_df['channel'] == j)]
            current_thresh = current_thresh['threshold']
            #print(current_thresh)
            current_thresh = current_thresh.iloc[0]
            #print(current_thresh)

            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'

            test = test[freq_cols]

            non_evoked_ls = []
            yes_evoked_ls = []
            for index, row in test.iterrows():
                for col in test.columns:
                    current_row = row


                    current_df = test[freq_cols].apply(lambda x: list(x)[index])


                    # the index will loop from 0 to 9 (which simply corresponds to each intensity)
                    current_db = db_ls[index]
                    #print(current_db)

                    if current_db < current_thresh:
                        # count every cell as non_evoked
                        #print(current_df)

                        no_evoke_bl_peak_fr = current_df.apply(lambda x: list(x)[7]) # just save them all, simple stuff
                        no_evoke_bl_peak_fr = no_evoke_bl_peak_fr.tolist()

                        for item in no_evoke_bl_peak_fr:
                            # save the bl corrected firing into a dataframe/ dictionary
                            non_evoked_ls.append(item)


                    if (current_db > current_thresh) and (current_db != np.max(db_ls)): # if above threshold and less than max intensity

                        #col = int(col)


                        current_cell = test.at[index,col]
                        current_below = test.at[(index+1),col]

                        max_freq_col = len(freq_ls)-1
                        

                        if col == 0: # if there is no column on the left then just reuse the center below
                            current_below_left = current_below
                        else:    
                            current_below_left = test.at[(index+1),(col-1)]

                        if col == max_freq_col: # if there is no column on the right then just reuse the center below
                            current_below_right = current_below
                        else:
                            current_below_right = test.at[(index+1),(col+1)]

                        #print(i,j)



                        try:
                            if current_cell[0] and (current_below[0] or current_below_left[0] or current_below_right[0]): # if current cell and one of its higher intensity neighbors are True evoked then add it to the evoked list
                                yes_evoke_bl_peak_fr = current_cell[7]
                                yes_evoked_ls.append(yes_evoke_bl_peak_fr)

                            else:
                                no_evoke_bl_peak_fr = current_cell[7]
                                non_evoked_ls.append(no_evoke_bl_peak_fr)
                                # if not add current cell to non-evoked list

                        except:
                            print('Problem with this unit/file')
                            


                    if current_db == np.max(db_ls):

                        true_columns = current_df.columns[current_df.iloc[0] == True]
                        current_df = current_df[true_columns]
                        yes_evoke_bl_peak_fr = current_df.apply(lambda x: list(x)[7])
                        yes_evoke_bl_peak_fr = yes_evoke_bl_peak_fr.tolist()

                        for item in yes_evoke_bl_peak_fr:
                            yes_evoked_ls.append(item)

                        
                        

                        # at the max intensity we can't check for contigiuos bins at a higher intensity
                        # so just take all cells that are evoked without checking


            #print(non_evoked_ls)    
            dprime = calc_dprime(yes_evoked_ls,non_evoked_ls)    

                # make a threshold dataframe columns being threshold, file, channel, genotype
            dprime_df = dprime_df._append({'dprime': dprime,
                                    'file': i,
                                    'channel':j,
                                    'genotype' : current_geno}, ignore_index=True)

    return dprime_df

# other idea is to generate an I/O function for each frequency to get threshold on a per column basis
# then take any cells above that threshold as evoked per column
# this is becuase the standard sound_evoked check leaves in too much noise


def savgol_dprime(evoked_df,db_ls,sanity):
    plot_df = pd.DataFrame(columns = ['dprime','file','channel','genotype'])

    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            #print(current_geno)
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'
            
            non_evoked_ls = [] # clear out the firing list
            yes_evoked_ls = []



            current_df = test[freq_cols] # need to get out the list values
            for col in current_df.columns:
                current_io = current_df[col]
                
                current_io = [sublist[7] for sublist in current_io]

                #print('File:',i,'Unit:',j,'Freq:',col)
                
                #print('Current IO',current_io)


                current_threshold = get_thresh(current_io,db_ls,sanity) # get threshold for current IO
                #print('Current Threshold',current_threshold)
                #thresh_ls = range(len(current_io)) # makes a list 0-X for x intensities
                #print(db_ls)
                #thresh_index = np.where(db_ls == current_threshold) # get index value of current threshold
                thresh_index = db_ls.index(current_threshold)
                #print(thresh_index) 
                non_evoked_fr = current_io[:thresh_index] # separate the firing values by evoked and non_evoked
                #print(non_evoked_fr)
                yes_evoked_fr = current_io[thresh_index+1:]# save into big list of evoked and non_evoked
                #print(yes_evoked_fr)

                if thresh_index > 0:                  
                    for item in yes_evoked_fr:
                                yes_evoked_ls.append(item)
                    for item in non_evoked_fr:
                                non_evoked_ls.append(item)

                
                

            dprime = calc_dprime(yes_evoked_ls,non_evoked_ls) # get dprime after done with all the columns for the unit



            plot_df = plot_df._append({'dprime':dprime,
                                                'file': i,
                                                'channel':j,
                                                'genotype' : current_geno}, ignore_index=True)
    return plot_df




def get_tuning_edges(current_khzx1_array,freq_ls,sanity):

# Smooth the line
    
    smoothed_y = savgol_filter(current_khzx1_array, window_length=5, polyorder=3)


    #thresh = 0.2*np.max(y)
    maxi=np.max(smoothed_y)
    mini = np.min(smoothed_y)
    addi = 0.5*(maxi-mini)
    thresh_y = np.min(smoothed_y)+addi
    #print(maxi,thresh_y)


    diff = smoothed_y - thresh_y

    sign_changes = np.where(np.diff(np.sign(diff)))[0] # find sign of the differences

    peaks = 'one'
    if len(sign_changes) > 2:
        peaks = 'multi'

    first_intercept = sign_changes[0]
    last_intercept = sign_changes[-1]


    thresh_khz_low = freq_ls[first_intercept]
    thresh_khz_high = freq_ls[last_intercept]
    

    
    


    all_evoked_cols = np.where(smoothed_y > thresh_y)[0]


    if sanity == 'yes':
        plt.figure(figsize=(5,5))
        plt.plot(current_khzx1_array, label='Original Data')
        plt.plot(smoothed_y,label = 'Smooth')
        plt.axhline(y = thresh_y, color='k', label='Inflection',linestyle='--')
        sns.despine()
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.legend()
        plt.title('Smoothing with 5 point filter')
        plt.show()

    
    #print(thresh_db)


    # use all_evoked_cells for bandwidth calculation for now
    return thresh_khz_low,thresh_khz_high,all_evoked_cols,peaks,smoothed_y




def calc_octaves(low,high):
    #low_khz = freq_ls[low] # convert index to freq in khz
	#high_khz = freq_ls[high]
	ratio = high/low # leave frequency in Hz if it is already, the q value code will convert it to kHz
	bandwidth = math.log2(ratio)
	return bandwidth

def get_qval(cf, bw):
     # calculates qvalue given characteristic freq index and bandwidth
      
	#cf_hz = freq_ls[cf] # convert cf index to cf in khz
	cf_khz =  cf/1000 # use if in Hz, comment out this line if in kHz
    
	qvalue = cf_khz/bw
     
    
	return qvalue


def get_bandwidth(evoked_df,thresh_df,freq_ls,db_ls,cf_df,sanity):
    # need to use the threshold value to get starting point
    # start at the row above threshold (10 dB in our case)
    # do the same sort of savgol filtering for the frequency curve
    # look by eye as to where the edges would be in the tuning curve heatmap and this savgol curve 
    # use the 20% (maybe more) cutoff to get the edges where each x axis point is a frequency
    # save upper and lower frequency numbers
    # calc octaves
    # save octave number

    # also add another measure of d-prime where we classify a frequency column for a given intesnity as evoked if it crosses a set threshold based on BF max firing
    # so then we can get an intensity dependent d-prime

    plot_df = pd.DataFrame(columns = ['dprime','low_freq','high_freq','BW','CF','db_above_threshold','peaks','file','channel','genotype'])

    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            


            cf_row = cf_df.loc[(cf_df['file'] == i) & (cf_df['channel'] == j)]
            
            current_cf = cf_row['CF']
            
            current_cf = current_cf.values[0]
            
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'

            current_df = test[freq_cols] # need to get out the list values

            for index, row in test.iterrows():  
                non_evoked_ls = [] # need to reset evoked lists because we are getting d-prime for every row
                yes_evoked_ls = []              
                current_df = test[freq_cols].apply(lambda x: list(x)[index])

                
                current_io = current_df.apply(lambda x: list(x)[7]) 
                #print(current_io)

                #print('File:',i,'Unit:',j)
                


                #print(index)

                current_edges = get_tuning_edges(current_khzx1_array = current_io,freq_ls = freq_ls,sanity=sanity)
                
                current_max = np.max(current_edges[2])
                current_min = np.min(current_edges[2])

                current_max_Hz = freq_ls[current_max]
                current_min_Hz = freq_ls[current_min]

                #print(current_min_Hz,current_max_Hz)

                current_bw = calc_octaves(current_min_Hz,current_max_Hz)
                #print(current_bw)




                current_threshold = thresh_df.loc[(thresh_df['file'] == i) & (thresh_df['channel'] == j)] # get threshold for current IO
                current_threshold = current_threshold['threshold']
                
                current_intensity = db_ls[index]
                

                relative_intensity = current_intensity - current_threshold
     
                relative_intensity = relative_intensity.values[0]
                

                

                all_evoked_indexes = current_edges[2] # gives every frequency value above the cutoff

                for a,item in enumerate(current_io): # if a frequency column for this intensity corresponds to an evoked index add to evoked list
                    if a in all_evoked_indexes:
                        yes_evoked_ls.append(item)
                    else:
                        non_evoked_ls.append(item)
   

                dprime = calc_dprime(yes_evoked_ls,non_evoked_ls) # get dprime after done with all the columns for the unit
                # can take the column values that are over the minimum to be evoked, and others to be not then calculate d-prime at each intensity above threshold


                
                
                
                #cf = float(cf)

                if current_bw == 0:
                     low = current_cf
                     high = current_cf + 1000
                     current_bw = calc_octaves(low,high)
                
                qval = get_qval(cf=current_cf,bw=current_bw)
                qval = qval/10 # correct for units


                plot_df = plot_df._append({'dprime':dprime,
                                        'low_freq':current_min_Hz,
                                        'high_freq':current_max_Hz,
                                        'BW':current_bw,
                                        'CF':current_cf,
                                        'Q-value':qval,
                                        'db_above_threshold':relative_intensity,
                                        'actual_dB': current_intensity,
                                        'peaks':current_edges[3],
                                        'file': i,
                                        'channel':j,
                                        'genotype' : current_geno}, ignore_index=True)

    return plot_df



def normalize(array,t_min,t_max):
    norm_arr = []
    diff = t_max - t_min
    diff_arr = max(array) - min(array)
    for item in array:
        temp = (((item - min(array))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return norm_arr


def area_under_curve(evoked_df,freq_ls,sanity):
    prec_df = pd.DataFrame(columns = ['file','channel','genotype'])
    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            #print(current_geno)
            #print(test)
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'

            current_df = test[freq_cols]
            #print(current_df)

            sum_all_freq = {}
            for series_name, series in current_df.items():
                current_freq = series

                


                
                current_sum = np.sum(current_cell[4] for current_cell in current_freq) # sum the baseline corrected frs for this column
                current_freq = freq_ls[series_name]
                sum_all_freq[current_freq] = current_sum
            y_edge_values = np.array(list(sum_all_freq.values()))
            
            edges = get_tuning_edges(y_edge_values,freq_ls,sanity)
            #print(edges)

            y_values = edges[4]
            
            y_values = normalize(y_values,0,1)


            x_values = np.arange(len(y_values))

            area = np.trapz(y_values, x_values)

            peaks, _ = find_peaks(y_values)

            #peaks = edges[4]
            peak_count = len(peaks)
            
            precision = (area*peak_count)         

                # make a threshold dataframe columns being threshold, file, channel, genotype
            prec_df = prec_df._append({'Precision': precision,
                                   'area':area,
                                   'peak_number':peaks,
                                    'file': i,
                                    'channel':j,
                                    'genotype' : current_geno}, ignore_index=True)


                
    return prec_df



# Define the Gaussian model
def gaussian(x, a, b, c):
    return a * np.exp(-((x - b) ** 2) / (2 * c ** 2))

# Define a function to calculate the slope and threshold at 20% of the range
def calculate_slope_and_threshold(params, x_data):
    a, b, c = params
    # Generate a dense set of x values for more accurate min/max detection
    x_fit = np.linspace(min(x_data), max(x_data), 1000)
    y_fit = gaussian(x_fit, *params)
    
    y_min = np.min(y_fit)
    y_max = np.max(y_fit)
    y_range = y_max - y_min
    y_threshold = y_min + 0.2 * y_range
    
    # Find the corresponding x value (intensity) where y_threshold is reached
    x_threshold = x_fit[np.argmin(np.abs(y_fit - y_threshold))]
    
    # Find the corresponding x value (intensity) where 80% of the range is reached
    y_asymp = y_min + 0.8 * y_range
    x_asymp = x_fit[np.argmin(np.abs(y_fit - y_asymp))]

    # Calculate slope based on 20% and 80% values of fitted line
    rise = y_asymp - y_threshold
    run = x_asymp - x_threshold
    gaussian_slope = rise / run

    return gaussian_slope, x_threshold, y_min, y_max

def get_gauss_params_fra(evoked_df,freq_ls,db_ls,sanity):
    # Example data
    data = evoked_df

    # Fit the model for each file and genotype
    files = data['file'].unique()
    results = []

    for file in files:
        subset_file = data[data['file'] == file]
        genotypes = subset_file['genotype'].unique()
        channels = subset_file['channel'].unique()

        for chan in channels:
            for genotype in genotypes:

                subset = subset_file[subset_file['genotype'] == genotype]
                subset = subset_file[subset_file['channel'] == chan]

                column_range_sum = subset.loc[:, freq_ls].sum(axis=1) # sum each column

                y_data = column_range_sum # y_data = sums
                x_data = freq_ls # x_data = freq_ls i.e. columns summed


                x_data = subset['dB'].values
                y_data = subset['spk_sec_abs'].values
                
                # Initial guesses for the parameters
                initial_guess = [max(y_data), np.mean(x_data), np.std(x_data)]
                
                try:
                    popt, pcov = curve_fit(gaussian, x_data, y_data, p0=initial_guess)
                    slope, threshold, y_min, y_max = calculate_slope_and_threshold(popt, x_data)
                    
                    result = {
                        'file': file,
                        'genotype': genotype,
                        'channel': chan,
                        'params': popt,
                        'slope': slope,
                        'threshold': threshold,
                        'y_min': y_min,
                        'y_max': y_max
                    }
                    results.append(result)
                    
                    # Plotting the fit
                    if sanity == 'yes':
                        plt.figure()
                        x_fit = np.linspace(min(x_data), max(x_data), 100)
                        y_fit = gaussian(x_fit, *popt)
                        plt.plot(x_fit, y_fit, label=f'{genotype} fit')
                        plt.axhline(y=y_min + 0.2 * (y_max - y_min), color='r', linestyle='--', label='20% threshold')
                        plt.axhline(y=y_min + 0.8 * (y_max - y_min), color='r', linestyle='--', label='80% threshold')
                        plt.xlabel('Sound Intensity')
                        plt.xticks(np.arange(0, 85, 5))
                        plt.ylabel('Firing Rate')
                        plt.title(f'File: {file}, Genotype: {genotype}')
                        plt.legend()
                        plt.show()
                    
                except RuntimeError:
                    print(f"Fit could not be performed for file {file}, genotype {genotype}")

    # Save results to a CSV file
    results_df = pd.DataFrame(results)
    results_df.to_csv('fit_results.csv', index=False)

    # Print the results
    for result in results:
        print(f"File: {result['file']}, Genotype: {result['genotype']}")
        print(f"  Parameters: {result['params']}")
        print(f"  Slope: {result['slope']}")
        print(f"  Threshold: {result['threshold']}")
        print(f"  Min Y: {result['y_min']}")
        print(f"  Max Y: {result['y_max']}")
        print()

# Define the Gaussian model
from scipy.optimize import curve_fit
   
def gaussian(x, a, b, c):
    return a * np.exp(-((x - b) ** 2) / (2 * c ** 2))

# Define a function to calculate the slope and threshold at 20% of the range
def calculate_slope_and_threshold(params, x_data):
    a, b, c = params
    # Generate a dense set of x values for more accurate min/max detection
    x_fit = np.linspace(min(x_data), max(x_data), 1000)
    y_fit = gaussian(x_fit, *params)
    
    y_min = np.min(y_fit)
    y_max = np.max(y_fit)
    y_range = y_max - y_min
    y_threshold = y_min + 0.5 * y_range
    print('Threshold Value: ',y_threshold)
    
    # Find the corresponding x value (intensity) where y_threshold is reached
    x_threshold = x_fit[np.argmin(np.abs(y_fit - y_threshold))]
    print(x_threshold)
    
    # Find the corresponding x value (intensity) where 80% of the range is reached
    y_asymp = y_min + 0.8 * y_range
    x_asymp = x_fit[np.argmin(np.abs(y_fit - y_asymp))]

    # Calculate slope based on 20% and 80% values of fitted line
    #rise = y_asymp - y_threshold
    #run = x_asymp - x_threshold
    #gaussian_slope = rise / run

    #gaussian_derivative = lambda x: -((x - b) / (c ** 2)) * gaussian(x, a, b, c)
    #gaussian_slope = gaussian_derivative(y_threshold)
    gaussian_slope = []
    h = 1e-5
    slope = np.abs((gaussian(x_threshold + h, a, b, c) - gaussian(x_threshold - h, a, b, c)) / (2 * h))
    gaussian_slope.append(slope)
    print('Slope:',gaussian_slope)
    

    return gaussian_slope, y_threshold, y_min, y_max

def make_gauss_df(evoked_df,freq_ls):
    gauss_df = pd.DataFrame(columns = ['file','channel','genotype'])
    for i in evoked_df['file'].unique(): # need ti group by file and by channel
        current_file = evoked_df.loc[evoked_df['file'] == i]
        for j in evoked_df['channel'].unique():
            test = current_file.loc[current_file['channel'] == j]
            current_geno = test['Genotype'][0] # grab the current genotype
            
            if isinstance(current_geno, str):
                current_geno = current_geno
            else:
                current_geno = 'WT'

            current_df = test[freq_cols]

            sum_all_freq = {}
            for series_name, series in current_df.items():
                current_freq = series

                                
                current_sum = np.sum(current_cell[4] for current_cell in current_freq) # sum the baseline corrected frs for this column
                current_freq = freq_ls[series_name]
                sum_all_freq[current_freq] = current_sum

                

                # make a threshold dataframe columns being threshold, file, channel, genotype
                gauss_df = gauss_df._append({'freq': current_freq,
                                            'sum': current_sum,
                                            'file': i,
                                            'channel':j,
                                            'genotype' : current_geno}, ignore_index=True)


                
    return gauss_df

def get_gauss_params_fra(gauss_df,sanity):
    # Example data
    data = gauss_df

    # Fit the model for each file and genotype
    files = data['file'].unique()
    results = []

    for file in files:
        subset_file = data[data['file'] == file]
        genotypes = subset_file['genotype'].unique()
        channels = subset_file['channel'].unique()

        for chan in channels:
            for genotype in genotypes:

                subset = subset_file[subset_file['genotype'] == genotype]
                subset = subset_file[subset_file['channel'] == chan]


                x_data = subset['freq'].values
                y_data = subset['sum'].values

                print(x_data)
                print(y_data)

                # Initial guesses for the parameters
                initial_guess = [max(y_data), np.mean(x_data), np.std(x_data)]
                
                try:
                    popt, pcov = curve_fit(gaussian, x_data, y_data, p0=initial_guess)
                    slope, threshold, y_min, y_max = calculate_slope_and_threshold(popt, x_data)
                    
                    result = {
                        'file': file,
                        'Genotype': genotype,
                        'channel': chan,
                        'params': popt,
                        'slope': slope,
                        'threshold': threshold,
                        'y_min': y_min,
                        'y_max': y_max
                    }
                    results.append(result)
                    
                    # Plotting the fit
                    if sanity == 'yes':
                        plt.figure()
                        x_fit = np.linspace(min(x_data), max(x_data), 100)
                        y_fit = gaussian(x_fit, *popt)
                        plt.plot(x_fit, y_fit, label=f'{genotype} fit')
                        plt.axhline(y=y_min + 0.5 * (y_max - y_min), color='r', linestyle='--', label='50% threshold')
                        plt.xlabel('Sound Intensity')
                        #plt.xticks(np.arange(0, 85, 5))
                        plt.ylabel('Firing Rate')
                        plt.title(f'File: {file}, Genotype: {genotype}')
                        plt.legend()
                        plt.show()
                    
                except RuntimeError:
                    print(f"Fit could not be performed for file {file}, genotype {genotype}")

    # Save results to a CSV file
    results_df = pd.DataFrame(results)
    results_df['slope'] = results_df['slope'].apply(lambda x: x[0]) 
    #results_df.to_csv('fit_results.csv', index=False)

    # Print the results
    for result in results:
        print(f"File: {result['file']}, Genotype: {result['Genotype']}")
        print(f"  Parameters: {result['params']}")
        print(f"  Slope: {result['slope']}")
        print(f"  Threshold: {result['threshold']}")
        print(f"  Min Y: {result['y_min']}")
        print(f"  Max Y: {result['y_max']}")
        print()

    return results_df



def get_distance_for_rlf(current_cell,edges,sanity):
    # time_range in form (-20,100)

    #spike_train = spk.SpikeTrain(np.array([0.1, 0.3, 0.45, 0.6, 0.9], [0.0, 1.0]))
    spike_trains = []

    for train in current_cell: # trains 31-60 are for the non_cf_cell
        #print(train)
        current_spike_train = spk.SpikeTrain(spike_times=train,edges=edges)
        spike_trains.append(current_spike_train) # need a list of spike trains

    
    
    ri_spike_dist = spk.spike_distance(spike_trains, RI=True)

    spike_distance_mat = spk.spike_distance_matrix(spike_trains)
    
    if sanity == 'yes':
        plt.figure(figsize=(5,5))

        plt.figure()
        
        plt.imshow(spike_distance_mat, interpolation='none')
        plt.title("SPIKE-distance")

        

    return ri_spike_dist,spike_distance_mat

def add_distance_col_and_df_rlf(raw_train_df,db_ls,freq_ls,edges,cf_df,sanity):

    
    form_spk_df = pd.DataFrame(columns = ['spike_train','spike_distance','spike_synch','file','channel','genotype'])
    intensity = db_ls

    for i in raw_train_df['file'].unique(): # need ti group by file and by channel
        current_file = raw_train_df.loc[raw_train_df['file'] == i]
        for j in current_file['channel'].unique():
            #print(j)
            current_unit = current_file.loc[current_file['channel'] == j]
            #print(current_unit)
            current_geno = current_unit['Genotype'] # grab the current genotype
            #print(current_geno)

            current_unit = current_unit[[0]]

            # compare evoked and non-evoked cells

            

            # loop over dB
            
            for index, row in current_unit.iterrows():
                current_int = current_unit.iloc[index]
                current_db = db_ls[index] # get current intensity]


                CF = cf_df.loc[cf_df['file'] == i] # get the CF cell from this data
                CF = CF.loc[CF['channel'] == j] 
                CF = CF['CF']
                CF = CF.iloc[0]
                cf_ind = freq_ls.index(CF)
                
                for series_name, series in current_int.items(): # loop over kHz
                    current_freq = freq_ls[series_name]  # get Hz freq  
                    current_cell =   current_int[series_name] # get currrent cell
                    current_cell = np.array(current_cell)

                    CF_octs = calc_octaves(CF,current_freq)

                    distance_cell = get_distance_for_rlf(current_cell = current_cell,edges=edges,sanity=sanity)


                        # make a threshold dataframe columns being threshold, file, channel, genotype
                    form_spk_df = form_spk_df._append({
                                                    'spike_train':current_cell,
                                                    'freq':current_freq,
                                                    'intensity':current_db,
                                                    'dist_matrix':distance_cell[1],
                                                    'ri_distance':distance_cell[0],
                                                            'file': i,
                                                            'channel':j,
                                                            'genotype' : current_geno}, ignore_index=True)
            
            

    return form_spk_df
