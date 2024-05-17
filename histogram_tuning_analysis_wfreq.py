#%% file import pre procesing
import pandas as pd
import numpy as np
import math

import io 
import requests 
import os
import glob
#plots
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as stats
from scipy.special import ndtri
from re import search
import tkinter
from tkinter import Tk
from tkinter.filedialog import askopenfilename
Tk().withdraw()

# frequency list
freq_ls = [1000.,  1500.,  2000.,  3000.,  4000.,  6000.,  8000., 12000.,
       16000., 20000., 24000., 30000., 35000., 40000., 45000., 50000.,
       55000., 60000., 65000.]

db_ls = [0,10,20,30,40,50,60,70,80,90]

def bandwidth_calc(low,high): # calculates distance in octaves given index input
	#low_khz = freq_ls[low] # convert index to freq in khz
	#high_khz = freq_ls[high]
	ratio = high/low
	bandwidth = math.log2(ratio)
	return bandwidth

def q_val(cf, bw): # calculates qvalue given characteristic freq index and bandwidth
	#cf_hz = freq_ls[cf] # convert cf index to cf in khz
	cf_khz =  cf/1000
	qvalue = cf_khz/bw
	return qvalue


tuning_df = pd.DataFrame(columns = ['threshold','CF','BW','Q10','Q20','Q30','low10','hi10','low20','hi20','lo30','hi30','file','channel','genotype','notes'])


# load in csvs with tuning info
filename_fr= askopenfilename(title="Select spks csv")


#ic_tuning = pd.read_csv("C:/Users/admin/Box/Behavior Lab/Shared/Manasi/Walker/fmr1_SD_analysis/ic_first_from_acx_p.csv")
ic_tuning_fr = pd.read_pickle(filename_fr)
#ic_tuning_fr = ic_tuning_fr.reset_index()
rows = len(db_ls)
cols = len(freq_ls)

col_names = ['{}'.format(col) for col in freq_ls]
row_names = ['{} dB'.format(row) for row in db_ls]
scale_ls = []
#fig, axes = plt.subplots(rows,cols,figsize=(20,30)) # need to adjust for how many conditions you are plotting


for x in ic_tuning_fr['file'].unique():

	current_file = ic_tuning_fr.loc[ic_tuning_fr['file'] == x]
	print(x)
	#geno_current = ic_tuning_fr['Genotype'][0]
	#print(geno_current[0])

	if 'WT' in x:
		geno_current = 'WT'
	if 'KO' in x:
		geno_current = 'KO'

	print(geno_current)

	for chan in current_file['channel'].unique():
		fig, axes = plt.subplots(rows,cols,figsize=(20,30)) # need to adjust for how many conditions you are plotting


		

		count = 0 # for selecting axis number in the loop
		i = 0
		j = 0
		for run in range(2):
		    for i in range(rows):
		        for j in range(cols):


		            test = current_file.loc[current_file['channel'] == chan]
		            test = test.iloc[:,:20]
		            testcell = test.at[i,j]
		            testcell_flat = np.hstack(testcell)

		            if run == 0:
		                #test = current_file.loc[current_file['channel'] == chan]
		                #test = test.iloc[:,:20]
		                #testcell = test.at[i,j]
		                #testcell_flat = np.hstack(testcell)
		                
		                #y, x, _ = axes[i,j].hist(testcell_flat,color='r',alpha = 0.5, bins = 150)
		                c,v = np.histogram(testcell_flat, bins = 150)
		                #print(y)
		               # scale_ls.append(y.max())
		                scale_ls.append(c.max())
		                


		                
            
		            if run == 1:
		                #test = current_file.loc[current_file['channel'] == chan]
		                #test = test.iloc[:,:20]
		                #testcell = test.at[i,j]
		                #testcell_flat = np.hstack(testcell)

		                y, a, _ = axes[i,j].hist(testcell_flat,alpha = 0.5, bins = 150)
		                axes[i,j].set(xticklabels=[])  # remove the tick labels
		                axes[i,j].tick_params(bottom=False)  # remove the ticks
		                #axes[i,j].set(yticklabels=[])  # remove the tick labels
		                axes[i,j].tick_params(labelleft=False, left=False)
		                axes[i,j].set_ylim(0,max(scale_ls))
		                axes[i,j].axvspan(0, 0.05, facecolor='g',alpha=.1,zorder=-100) #new line

		                axes[i,j].spines['top'].set_visible(False)
		                axes[i,j].spines['right'].set_visible(False)
		                axes[i,j].spines['bottom'].set_visible(False)
		                axes[i,j].spines['left'].set_visible(False)



		

		                for ax, col in zip(axes[0], col_names):
		                    ax.set_title(col)

		                for ax, row in zip(axes[:,0], row_names):
		                    ax.set_ylabel(row, rotation=0)

		plt.show(block=False) #need this to take input while plot is up
		                


		    

		while True:
			try:
				print(x)

				current_thresh = int(input('Threshold:'))
				current_cf = int(input('CF:'))

				# range of indexes that will correspond to frequency list
				current_bw10 = bandwidth_calc(int(input('BW10 Lower:')),int(input('BW10 Upper:')))
				current_bw20 = bandwidth_calc(int(input('BW20 Lower:')),int(input('BW20 Upper:')))
				current_bw30 = bandwidth_calc(int(input('BW30 Lower:')),int(input('BW30 Upper:')))

				# get frequency absolute values
				current_lo10
				current_hi10
				current_lo20
				current_hi20
				current_lo30
				current_hi30

				# any comments or notes to make
				notes = input('Notes: ')

				# convert to qvalues
				curr_q10 = q_val(current_cf,current_bw10)
				curr_q20 = q_val(current_cf,current_bw20)
				curr_q30 = q_val(current_cf,current_bw30)
				break
			except:
				print('Invalid input type. Please use ints')
		plt.close(fig)# clear current plot
		

		# add everything to dataframe

		dict = {'threshold':[current_thresh],
		        'CF':[current_cf],
		        'BW':[current_bw10],
		        #'BW20':[current_bw20],
		        #'BW30':[current_bw30],
		        'Q10':[curr_q10],
		        'Q20':[curr_q20],
		        'Q30':[curr_q30],
		        'low10':[],
		        'hi10':[],
		        'low20':[],
		        'hi20':[],
		        'lo30':[],
		        'hi30':[],
		        'file':[x],
		        'channel':[chan],
		        'genotype':[geno_current],
		        'notes':[notes]
		       }
		curr_df = pd.DataFrame(dict)
		tuning_df = tuning_df.append(curr_df)
		tuning_df.to_csv('your/path/tofile_name.csv')

