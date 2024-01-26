#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:15:20 2024

@author: mdouaihy, IGMM-CNRS / LPHI-U Montpellier
"""



import sys
sys.path.append('./utilities/')
from read_rawDataFiles import readDataInFile
from movies_combining import movies_combining_rawData
#sys.path.append('/home/mdouaihy/deconv_python/')
#sys.path.append('/home/mdouaihy/snail_crispr/code/')
#sys.path.append('/home/mdouaihy/bcpd_test/code/')
#sys.path.append('/home/mdouaihy/DeconvolutionPython+Notebook/utilities/')
import numpy as np
import os
import matplotlib.pyplot as plt
import math
import seaborn as sns
from MLE_gamma import MLE_gamma
from gregory_GaussianUnknownMeanAndVariance import GaussianUnknownMeanUnknownVariance
from bcpd_attribute import bocd
from retriving_time_points import retriving_time_100723
import pandas as pd

plt.close('all')

inputpath = '/home/mdouaihy/MLL_data/sna_crispr_paper/code/Data/'
data_type = '' # if we wan t to add a specific name to the output folder otherwise it's generate automatically
extension = '.xlsx' # extension type of the data (either xlsx or xls)
outputpath =  '/home/mdouaihy/MLL_data/sna_crispr_paper/code/Results/'
    


percentage = 0.6 # time point after percentage of max intensity
kernel = 2 #  kernel of smoothing the data
filtered = 1 # filter nuclei that goes into repression
p = 0.8 # probability of having a changepoint >=p

""" parameters of the signal and the movie """

FrameLen = 4.64 # frame length in seconds or time resolution
retention = 0 # in seconds
Polym_speed = 25 # Polymerase speed'
TaillePreMarq = 1649 # length of mRNA until we reach the beginning of the MS2 sequence in bp
TailleSeqMarq = 1292 # length of MS2 sequence in bp
TaillePostMarq = 1312 + Polym_speed*retention # length of mRNA from last loop of MS2 until the polymerase leaves the TS
EspaceInterPolyMin = 30 # minimal distance between polymerase (in bp)
Intensity_for_1_Polym = 1 # calibration factor ( 1 is the data is already divided by the calibration)    
DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed # duration of the signal
FreqEchImg = 1/FrameLen  
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed)
    



""" plotting parameters """
nbr_subplot_per_fig = 36
lw = 0.6 #line width of the plot



    
######################################################
if not os.path.exists(outputpath):
    os.mkdir(outputpath)    
    

fParam =  outputpath +  'HyperParameters.npz'
np.savez(fParam, 
          Polym_speed = Polym_speed,  
          TaillePreMarq = TaillePreMarq,
          TailleSeqMarq = TailleSeqMarq,
          TaillePostMarq = TaillePostMarq,
          EspaceInterPolyMin = EspaceInterPolyMin,
          FrameLen = FrameLen,
          Intensity_for_1_Polym = Intensity_for_1_Polym,
          FreqEchImg = FreqEchImg,
          DureeSignal = DureeSignal,
          FreqEchSimu = FreqEchSimu,    
          retention = retention
        )



readDataInFile(inputpath, outputpath, data_type + '/', fParam, extension) # transform the excel sheets to npz files

#npzFile = 


dirwrite = outputpath + '/CP_' + str(percentage) + 'xMax_K_' + str(kernel) + '/'
if filtered:
    data_type = '_bcpd_' + str(percentage) + 'xMax_K_' + str(kernel) + '_filtered' + '/'
    
else:
    data_type = '_bcpd_' + str(percentage) + 'xMax_K_' + str(kernel) + '/'

if not os.path.exists(dirwrite):
    os.mkdir(dirwrite)
    
    
          
dirwrite_full = dirwrite + '/full_signal/'
if not os.path.exists(dirwrite_full):
    os.mkdir(dirwrite_full)
    
dirwrite_activation = dirwrite + '/activation/'
if not os.path.exists(dirwrite_activation):
    os.mkdir(dirwrite_activation)

dirwrite_repressed = dirwrite + '/repressed/'
if not os.path.exists(dirwrite_repressed):
    os.mkdir(dirwrite_repressed)
    
hazard = 1/10000  # Constant prior on changepoint probability.
    


#######################################################


""" combining movies """
[FrameLen, DataExp, tmax_combined, tstart] = movies_combining_rawData(outputpath + 'npzFile/', inputpath, Intensity_for_1_Polym, extension)
n2 = DataExp.shape
nexp = n2[1]
frame_num = n2[0]


""" extracting files names """

file_name_list = np.array(os.listdir(inputpath )) # list of the data
nfiles = len(file_name_list) # length of the list
file_name_list = [file_name_list[i].replace(extension,'') for i in range(nfiles)]
        
nname = []        
n0 = 0        
for iii in range(len(file_name_list)):    
    fname = inputpath + file_name_list[iii] + extension
    
    ### loading time into mitosis    
    rawData = pd.read_excel(fname)
    rawData_np = rawData.iloc[1:, 4:].values
    time_into_mitosis = rawData['Time'].dropna(how='all').to_numpy()
    n2 = rawData_np.shape
    
    if len(time_into_mitosis) ==0:
        true_FrameLen =  3.86 
        print('!!!!! framelen pre-assigned to 3.86')
        tstart = 0
        tend = n2[0]*3.86
    else:
        true_FrameLen = np.unique(np.round(np.diff(time_into_mitosis),3))
        tstart = time_into_mitosis[0]
        tend = time_into_mitosis[-1]
        
    if tstart==0:
        DataExp_movie_i = DataExp[:n2[0],n0 :n2[1] + n0]
    else:
        DataExp_movie_i = DataExp[round(tstart/FrameLen) :round(tstart/FrameLen) + n2[0],n0 :n2[1] + n0]
    chck_pt = np.where(rawData_np[:DataExp.shape[0]-round(tstart/FrameLen), :] != DataExp_movie_i)[0]
    
    if len(chck_pt)!=0:
        print('!!!! names incoorect')
    
    else:
        splt_names = rawData.columns[4:].tolist()
        
        file_name = file_name_list[iii]
        file_name = file_name.replace('CalibratedTraces','')
        file_name = file_name.replace('__','_')
        file_name = file_name[:-1]
        nname = nname + [file_name + ' ' + s for s in splt_names]
        
    n0 = n0+n2[1]    

""" prior of normal gamma distribution parameters """
k0 = 1.5    # The prio
mean0 = np.mean(DataExp[int(5*60/FrameLen):int(10*60/FrameLen),:])
[alpha0, beta0] = MLE_gamma(np.var(DataExp[int(5*60/FrameLen):int(10*60/FrameLen),:],axis=0))


set_parameters_exp = [mean0, k0, alpha0, beta0]




    
######## Data needed for bcpd
DataExp_1st_part = np.full((DataExp.shape[0],DataExp.shape[1]),np.nan)
DataExp_2nd_part = np.full((DataExp.shape[0],DataExp.shape[1]),np.nan)
T0_repressed = []
R_exp_list = []




def smooth(y, box_pts): # the kernel is nomalized (box blur)
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

for data_i in range(nexp):
    data_exp = smooth(DataExp[:,data_i],kernel)
    set_parameters_exp = [mean0, alpha0, k0, beta0]
    
    model_exp          = GaussianUnknownMeanUnknownVariance(set_parameters_exp[0], set_parameters_exp[1], set_parameters_exp[2], set_parameters_exp[3])
    R_exp = bocd(data_exp, model_exp, hazard)
    [x_axis_lines, proba_line] = retriving_time_100723(R_exp,p)
    R_exp_list.append(R_exp)

    
    
    cond_for_cpd = np.where(data_exp>=np.nanmax(data_exp)*percentage)[0][-1] + 1
    cp_start_point = [x[0] for x in x_axis_lines]
    indx_bcpd_pt = np.where((cp_start_point-cond_for_cpd)>=0)[0]
    if len(indx_bcpd_pt) ==0:
        bcpd_pt = len(data_exp)
    else:
        indx_bcpd_pt = indx_bcpd_pt[0]
    
        bcpd_pt = x_axis_lines[indx_bcpd_pt][0]

    DataExp_1st_part[:bcpd_pt, data_i] = DataExp[ :bcpd_pt,data_i]
    DataExp_2nd_part[:DataExp.shape[0]-bcpd_pt, data_i] = DataExp[bcpd_pt:,data_i]
    
    
    T0_repressed.append(bcpd_pt)


    """ plotting """
    ifig= math.floor((data_i)/nbr_subplot_per_fig)+1
#    h = plt.figure( ifig )#, figsize=[10,12]   
    
    h1 = plt.figure(ifig, figsize=(13.69,10.27))  
    h2 = plt.figure((1+ifig)*100, figsize=(13.69,10.27))  
    h3 = plt.figure((1+ifig)*1000, figsize=(13.69,10.27))  
    

    
    subplot1 = h1.add_subplot(round(np.sqrt(nbr_subplot_per_fig)),round(np.sqrt(nbr_subplot_per_fig)),(data_i%nbr_subplot_per_fig+1))
    subplot2 = h2.add_subplot(round(np.sqrt(nbr_subplot_per_fig)),round(np.sqrt(nbr_subplot_per_fig)),(data_i%nbr_subplot_per_fig+1))
    subplot3 = h3.add_subplot(round(np.sqrt(nbr_subplot_per_fig)),round(np.sqrt(nbr_subplot_per_fig)),(data_i%nbr_subplot_per_fig+1))

    h1.subplots_adjust(hspace = 0.5, wspace=0.3)
    h2.subplots_adjust(hspace = 0.5, wspace=0.3)
    h3.subplots_adjust(hspace = 0.5, wspace=0.3)
 
    
    x = DataExp_1st_part[:, data_i]
    subplot1.plot(np.arange(len(x[~np.isnan(x)]))/FreqEchImg/60,x[~np.isnan(x)], color = 'k', linewidth = lw+0.1)
    subplot2.plot(np.arange(len(x[~np.isnan(x)]))/FreqEchImg/60,x[~np.isnan(x)], color = 'k', linewidth = lw)
    
    x = DataExp_2nd_part[:,data_i]
    subplot1.plot(np.arange(bcpd_pt,bcpd_pt + len(x[~np.isnan(x)]))/FreqEchImg/60,x[~np.isnan(x)], color = 'b', linewidth = lw+0.1)
    subplot3.plot(np.arange(bcpd_pt,bcpd_pt + len(x[~np.isnan(x)]))/FreqEchImg/60, x[~np.isnan(x)], color = 'b', linewidth = lw)
    
    subplot1.plot(np.arange(bcpd_pt)/FreqEchImg/60, data_exp[:bcpd_pt], color = 'r', linewidth = lw-0.4)



    subplot1.tick_params(axis='both', which='major', labelsize=5)
    
    subplot2.tick_params(axis='both', which='major', labelsize=10)
    
    subplot3.tick_params(axis='both', which='major', labelsize=5)
    
    sns.despine()
    
    h1.suptitle(str(percentage*100) + '% : full intensity')
    h2.suptitle(str(percentage*100) + '% : activation part')
    h3.suptitle(str(percentage*100) + '% : repressed part')
    
    if data_i%nbr_subplot_per_fig==nbr_subplot_per_fig-1 or (data_i+1)==nexp:
        figfile_repressed = dirwrite_repressed +'/figure_repressed_'+str(ifig)+'.pdf'
        h3.savefig(figfile_repressed)    
        plt.close(h3)
        
        figfile_activation = dirwrite_activation +'/figure_activation_'+str(ifig)+'.pdf'
        h2.savefig(figfile_activation)    
        plt.close(h2)
        
        figfile_full = dirwrite_full +'/figure_full_'+str(ifig)+'.pdf'
        h1.savefig(figfile_full)    
        plt.close(h1)        

nan_rows = np.where(np.isnan(np.sum(DataExp_1st_part, axis=1)))[0]

if len(nan_rows)==0:
    DataExp_0 = DataExp_1st_part
else:
    DataExp_0 = DataExp_1st_part[:nan_rows[0],:]

nan_rows = np.where(np.nansum(DataExp_2nd_part, axis=1)==0)[0]
        
if len(nan_rows)==0:
    DataExp_1 = DataExp_2nd_part
else:
    DataExp_1 = DataExp_2nd_part[:nan_rows[0],:]



T0_repressed = np.array(T0_repressed)



np.savez(dirwrite + '/result_BCPD_filtered_' + str(filtered) + '.npz', 
         DataExp_activation = DataExp_0,
         DataExp_repressed = DataExp_1,
         DataExp = DataExp,
         
         T0_repressed = (T0_repressed-1)/FreqEchImg)

np.savez(dirwrite + '/full_bcpd_results_p_' + str(p) + '.npz', 
         DataExp = DataExp,
         R_exp_list = R_exp_list)

if filtered ==1:

    """ check points to eliminate nuclei in two steps:
        1- remove 0 nuclei
        2- remove nuclei that has less then 10% of max number of pol II per zone """

    check_nan_columns = np.where(np.sum(DataExp_0==0, axis =0)+np.sum(np.isnan(DataExp_0), axis = 0)==DataExp_0.shape[0])[0] # this checks nuclei that are full with 0 and nan

    
    """ remove nuclei with less than 1 pol II """
    sd=DataExp_0.shape
    frame_num=sd[0] ### number of frames
    DureeSimu = frame_num*FrameLen  ### film duration in s
    DureeAnalysee = DureeSignal + DureeSimu
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)) # maximal number of polymerase positions
    area = FreqEchImg/Polym_speed*Intensity_for_1_Polym*(TaillePostMarq+TailleSeqMarq/2)
    DataExpSmooth = np.minimum(np.round(np.sum(DataExp_0,axis = 0) / area), num_possible_poly)
    check_less_than_pol_II = np.where(DataExpSmooth==0)[0]

    ### where to write figure files     
    delete_elemts_0 = np.unique(np.append(check_nan_columns,check_less_than_pol_II))
    
    for data_i in range(len(delete_elemts_0)):
        ifig= math.floor((data_i)/9)+1
        h = plt.figure(ifig)#, figsize=[10,12]   
#                plt.subplots_adjust(hspace=1,wspace=1)
        plt.subplot(3,3,(data_i%9+1))
        plt.plot(np.arange(0, frame_num)/FreqEchImg/60, DataExp_0[:,data_i].T, color = 'k', linewidth = 0.1)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        sns.despine()
        if data_i%9==8 or (data_i+1)==len(delete_elemts_0):
            figfile=dirwrite_activation +'/deleted_figure_'+str(ifig)+'.pdf'
            plt.tight_layout()
            h.savefig(figfile)
            plt.close()
            
    DataExp_0 = np.delete(DataExp_0, delete_elemts_0, 1)


    check_nan_columns = np.where(np.sum(DataExp_1==0, axis =0)+np.sum(np.isnan(DataExp_1), axis = 0)==DataExp_1.shape[0])[0] # this checks nuclei that are full with 0 and nan

    
    """ remove nuclei with less than 1 pol II """
    sd=DataExp_1.shape
    frame_num=sd[0] ### number of frames
    DureeSimu = frame_num*FrameLen  ### film duration in s
    DureeAnalysee = DureeSignal + DureeSimu
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)) # maximal number of polymerase positions
    area = FreqEchImg/Polym_speed*Intensity_for_1_Polym*(TaillePostMarq+TailleSeqMarq/2)
    DataExpSmooth = np.minimum(np.round(np.sum(DataExp_1,axis = 0) / area), num_possible_poly)
    check_less_than_pol_II = np.where(DataExpSmooth==0)[0]
    


    ### where to write figure files     
    delete_elemts_1 = np.unique(np.append(check_nan_columns,check_less_than_pol_II))
    
    for data_i in range(len(delete_elemts_1)):
        ifig= math.floor((data_i)/9)+1
        h = plt.figure(ifig)#, figsize=[10,12]   
#                plt.subplots_adjust(hspace=1,wspace=1)
        plt.subplot(3,3,(data_i%9+1))
        plt.plot(np.arange(0, frame_num)/FreqEchImg/60, DataExp_1[:,data_i].T, color = 'k', linewidth = 0.1)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        sns.despine()
        if data_i%9==8 or (data_i+1)==len(delete_elemts_1):
            figfile=dirwrite_repressed+'/deleted_figure_'+str(ifig)+'.pdf'
            plt.tight_layout()
            h.savefig(figfile)
            plt.close()
            
    DataExp_1 = np.delete(DataExp_1, delete_elemts_1, 1)
    
    f = open(dirwrite + '/removed_nuclei', "a")
    if sum(delete_elemts_0)!=0:
        f.write("\n "+ 'Activation : ' + str(delete_elemts_0))
    if sum(delete_elemts_1)!=0:
        f.write("\n "+ 'Repression : ' + str(delete_elemts_1))    
    f.close()
            
   
np.savez(dirwrite + '/result_BCPD_filtered_' + str(filtered) + '.npz', 
         DataExp_activation = DataExp_0,
         DataExp_repressed = DataExp_1,
         DataExp = DataExp,
         
         T0_repressed = (T0_repressed-1)/FreqEchImg)        

plt.close('all')    

bcpd = str(percentage) + '_K_' + str(kernel)
T0_0 = (T0_repressed-1)/FreqEchImg/60
delete_elemts_1 = np.unique(np.append(check_nan_columns,check_less_than_pol_II))
T0_1 = np.delete(T0_0, delete_elemts_1)

if len(T0_0)!= len(nname):
    nname_1 = ['nuclei ' + str(i) for i in range(len(T0_1))]
    nname = ['nuclei ' + str(i) for i in range(len(T0_0))]
else:
    nname_1 = np.delete(nname, delete_elemts_1)
            
    
h, ax = plt.subplots(2,1)
h.suptitle(('CP_' + bcpd))
ax[0].hist(T0_0, bins = 50)
ax[0].set_title('All nuclei')

ax[1].hist(T0_1, bins = 50)
ax[1].set_title('fitered nuclei')

h.tight_layout()
h.savefig(dirwrite +  '/T0_density_1.pdf')


h, ax = plt.subplots(2,1)
h.suptitle(('CP_' + bcpd).replace('_', ' '))
sorted_T0_0 = np.sort(T0_0)
cdf_T0 = np.arange(1, len(sorted_T0_0) + 1) / len(sorted_T0_0)
T0_0_half = sorted_T0_0[np.where(cdf_T0>= 0.5)[0][0]]
ax[0].plot(sorted_T0_0, cdf_T0, marker='.', linestyle='none')
ax[0].plot([T0_0_half,T0_0_half], [0,1], 'r')
ax[0].set_xlabel('Data')
ax[0].set_ylabel('CDF')
ax[0].set_title('All nuclei')
#        ax[0].set_xlim(0,37)

sorted_T0_1 = np.sort(T0_1)
cdf_T1 = np.arange(1, len(sorted_T0_1) + 1) / len(sorted_T0_1)
T0_1_half = sorted_T0_1[np.where(cdf_T1>= 0.5)[0][0]]
ax[1].plot(sorted_T0_1, cdf_T1, marker='.', linestyle='none')
ax[1].plot([T0_1_half,T0_1_half], [0,1], 'r')
ax[1].set_xlabel('Data')
ax[1].set_ylabel('CDF')
ax[1].set_title('fitered nuclei')
h.tight_layout()
h.savefig(dirwrite + '/T0_cdf_1.pdf')


df_0 = pd.DataFrame({'All nuclei: names': nname, 'All nuclei: T0': T0_0})
df_1 = pd.DataFrame({'Filtered: names': nname_1, 'Filtered: T0':T0_1})


df = pd.concat([df_0,df_1], ignore_index=True, axis=1)
df.columns = ['All nuclei: names', 'All nuclei: T0', 'Filtered: names', 'Filtered: T0']
df.insert(2,'','')
xlsfilename = dirwrite + '/T0_' + 'CP_' + bcpd +  '.xlsx'
writer = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
df.to_excel(writer, index=False)

writer.save()