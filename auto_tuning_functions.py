# behavior figures

# looking at bulk tuning curves in one plot for ACx and IC

import pandas as pd
import numpy as np



#plots
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from array import *

print(pd.__version__)

 
                                                                                                        



                         

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
                                                                                       
                                                                                    
def is_significantly_activated(current_cell):
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
    #print('Response Duration: ', duration_rounded)
    
    # sanity plots if needed
    #plt.figure(figsize=(8, 6))
    #plt.hist(current_cell, bins=bin_edges, color='blue', edgecolor='black', alpha=0.7)
    #plt.axvline(x=peak_latency, color='red', linestyle='--', label=f'Peak Value: {peak_latency:.2f}')
    #plt.axhline(y=prestim_sd4,color='black',linestyle='--')
    #plt.axvspan(0, 0.1, facecolor='g',alpha=.1,zorder=-100) #new line
    #plt.title('SD4: ',prestim_sd4)
    
    #plt.show()
    
    evoked_ls = [evoked_Status, peak_latency, peak_fr, duration_rounded,bl_peak_fr,first_bin_exceeding_threshold,last_bin]
    
    return evoked_ls


def make_evoked_df(current_df,freq_ls,db_ls):

    

    rows = len(db_ls)
    cols = len(freq_ls)

    sound_evoked_array = []
    non_evoked_array = []

    big_df = current_df

    new_df = pd.DataFrame(columns = big_df.columns.values) # placeholder dataframe to append to
                        # not efficient at all, but I don't want to wrestle with replacing by index either


    for x in big_df['file'].unique(): # loop over each tuning curve
        current_file = big_df.loc[big_df['file'] == x]
        
        for chan in current_file['channel'].unique(): 
            test = current_file.loc[current_file['channel'] == chan]
            
            current_geno = test['Genotype'][0] # grab the current genotype
            spiking_only_df = test[range(cols)] # grab only the columns that correspond to freqs
            #print(chan,x,current_geno)
            #print(spiking_only_df)

            current_output_df = spiking_only_df.map(lambda x: is_significantly_activated(x))  
            
            #print(current_output_df)
            current_output_df['file'] = x
            current_output_df['channel'] = chan
            current_output_df['Genotype'] = test['Genotype']
            
            new_df = pd.concat([new_df,current_output_df],ignore_index=False)
        
    return new_df
            

def make_dB_df(evoked_df):
    plot_df = pd.DataFrame(columns = ['latency','abs_peak_fr','rel_peak_fr','resp_duration','first_bin','last_bin','smoothed_data','file','channel','genotype'])

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
            
            
            for index, row in test.iterrows():
                
                current_df = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]].apply(lambda x: list(x)[index])# need to get only cells with True
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
                print(peak_fr)
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


                plot_df = plot_df._append({'latency': lat_mean, 'abs_peak_fr': peak_mean,
                                                    'rel_peak_fr':rel_peak_mean,
                                                    'resp_duration': dur_mean,
                                                    'first_bin':first_bin,
                                                    'last_bin':last_bin,
                                                    'smoothed_data':smoothed_mean,
                                                    'file': i,
                                                    'channel':j,
                                                    'genotype' : current_geno}, ignore_index=True)
    return plot_df           

def get_thresh(current_dbx1_array,db_ls):

# Smooth the line
    
    smoothed_y = savgol_filter(current_dbx1_array, window_length=5, polyorder=3)


    #thresh = 0.2*np.max(y)
    maxi=np.max(smoothed_y)
    addi = 0.1*maxi
    thresh_y = np.min(smoothed_y)+addi
    #print(thresh_y)

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

    
    thresh_db = db_ls[thresh_ind]
    #print(thresh_db)

    return thresh_db

def add_thresh_col(dB_df,db_ls):

    
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
     
            thresh = get_thresh(curr_fr,intensity)

                # make a threshold dataframe columns being threshold, file, channel, genotype
            thresh_df = thresh_df._append({'threshold': thresh,
                                                    'file': i,
                                                    'channel':j,
                                                    'genotype' : current_geno}, ignore_index=True)

    return thresh_df
    

# add to get thresh for sanity plots
#plt.figure(figsize=(5,5))
#plt.plot(data, label='Original Data')
#plt.plot(smoothed_y,label = 'Smooth')
#plt.scatter(inflection_indices, smoothed_y[inflection_indices], color='red', label='Inflection Points')
#plt.axhline(y = thresh, color='k', label='Inflection',linestyle='--')
#sns.despine()
#plt.xlabel('Index')
#plt.ylabel('Value')
#plt.legend()
#plt.title('Smoothing with 5 point filter')
#plt.show()





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

            #current_df = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]].apply(lambda x: list(x)[0])
            #current_df = current_df.apply(lambda x: list(x)[4])
            current_df = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]]
            #print(current_df)

            sum_all_freq = {}
            for series_name, series in current_df.items():
                current_freq = series

                


                #current_df = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]].apply(lambda x: list(x)[0])
                
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

            #current_df = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]]
            #print(current_df)

            test = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]]

            non_evoked_ls = []
            yes_evoked_ls = []
            for index, row in test.iterrows():
                for col in test.columns:
                    current_row = row


                    current_df = test[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]].apply(lambda x: list(x)[index])


                    # the index will loop from 0 to 9 (which simply corresponds to each intensity)
                    current_db = db_ls[index]
                    #print(current_db)

                    if current_db < current_thresh:
                        # count every cell as non_evoked
                        #print(current_df)

                        no_evoke_bl_peak_fr = current_df.apply(lambda x: list(x)[4]) # just save them all, simple stuff
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
                                yes_evoke_bl_peak_fr = current_cell[4]
                                yes_evoked_ls.append(yes_evoke_bl_peak_fr)

                            else:
                                no_evoke_bl_peak_fr = current_cell[4]
                                non_evoked_ls.append(no_evoke_bl_peak_fr)
                                # if not add current cell to non-evoked list

                        except:
                            print('Problem with this unit/file')
                            


                    if current_db == np.max(db_ls):

                        true_columns = current_df.columns[current_df.iloc[0] == True]
                        current_df = current_df[true_columns]
                        yes_evoke_bl_peak_fr = current_df.apply(lambda x: list(x)[4])
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

def get_bandwidth():

    return


def get_qvals():

    return

