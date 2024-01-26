#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:55:45 2023

@author: mdouaihy
"""


import numpy as np

import pandas as pd

import os


def movies_combining(xlsfile,respath, filesList, fParam, extension):
    ###############################################
        
    fparam = np.load(fParam)
    
    FreqEchSimu = fparam['FreqEchSimu']
    FreqEchImg = fparam['FreqEchImg']
    TailleSeqMarq = fparam['TailleSeqMarq']
    TaillePostMarq = fparam['TaillePostMarq']
    Polym_speed = fparam['Polym_speed']
    EspaceInterPolyMin = fparam['EspaceInterPolyMin']
    FrameLen = fparam['FrameLen']
    DureeSignal = fparam['DureeSignal']
        
    
    tend = np.zeros((len(filesList)))
    tstart = np.zeros((len(filesList)))
    true_FrameLen = np.zeros((len(filesList))) # some movies may have a different time resolution
    nfiles_l = np.ones((len(filesList)))
    
    
    ### this iteration is to fix the total number of nuclei and movie length
    ###d and also double checking the removal of 0 nuclei 
    for iii in range(len(filesList)):
        fname = respath + 'result_' + filesList[iii].replace('_','') + '.npz'
        ffiles = np.load(fname)
        DataExp = ffiles['DataExp']
        DataPred = ffiles['DataPred']
        PosPred = ffiles['PosPred']
        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        
        n2 = DataExp.shape
    
        fname = xlsfile + filesList[iii] + extension
        
        ### loading time into mitosis    
        rawData = pd.read_excel(fname)
        rawData.columns = [x.lower() for x in rawData.columns]
        time_into_mitosis = rawData['time'].dropna(how='all').to_numpy()
        
        
        if len(time_into_mitosis) ==0:
            true_FrameLen[iii] =  3.86 
            print('!!!!! framelen pre-assigned to 3.86')
            tstart[iii] = 0
            tend[iii] = n2[0]*3.86
        else:
            true_FrameLen[iii] = np.unique(np.round(np.diff(time_into_mitosis),3))
            tstart[iii] = time_into_mitosis[0]
            tend[iii] = time_into_mitosis[-1]
            
        nfiles_l[iii] = n2[1]
    
    
    ##############
    
    
    nnuclei = round(sum(nfiles_l)) # total number of nucleis
    cumsumnuclei = np.insert(np.cumsum(nfiles_l),0,0)
    cumsumnuclei = cumsumnuclei.astype(int)

    
    frame_num_combined = round(min(tend/true_FrameLen))
    
    
    DureeSimu = frame_num_combined*FrameLen ### film duration in s
    DureeAnalysee = DureeSimu + DureeSignal;
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed))
    
    DataExp_combined = np.zeros((frame_num_combined, nnuclei))
    DataPred_combined= np.zeros((frame_num_combined, nnuclei))
    PosPred_combined= np.zeros((num_possible_poly, nnuclei))
    nn_combined = PosPred_combined.shape
    tmax_combined = np.ones((nnuclei,))
 
    
    #### combining movies according to their time into mitosis
    for iii in range(len(filesList)):

        fname = respath + 'result_' + filesList[iii].replace('_','') + '.npz'
        ffiles = np.load(fname)
        DataExp = ffiles['DataExp']
        DataPred = ffiles['DataPred']
        PosPred = ffiles['PosPred']
        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        DataPred = np.delete(DataPred,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        DataPred = np.delete(DataPred, check_nan_columns, 1)
        PosPred = np.delete(PosPred, check_nan_columns, 1)
        
        n2 = DataExp.shape
        n3 = PosPred.shape   
    ############## concatenating nuclei according to time into mitosis
        if tstart[iii]==0:

            DataExp_combined[:n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:frame_num_combined,:]
            DataPred_combined[:n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataPred[:frame_num_combined,:]
            t0_posPred = 0;
            PosPred_combined[t0_posPred: t0_posPred+n3[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = PosPred[:num_possible_poly,:]
    
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]

        else:

            DataExp_combined[round(tstart[iii]/FrameLen): round(tstart[iii]/FrameLen)+n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:frame_num_combined-round(tstart[iii]/FrameLen), :]
            DataPred_combined[round(tstart[iii]/FrameLen): round(tstart[iii]/FrameLen)+n2[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataPred[:frame_num_combined-round(tstart[iii]/FrameLen), :]
            t0_posPred = round(tstart[iii]*FreqEchSimu)
            
            if (t0_posPred+n3[0]-num_possible_poly) <=0:
                PosPred_combined[t0_posPred: t0_posPred+n3[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = PosPred[:num_possible_poly-t0_posPred, :]
            else:
                PosPred_combined[t0_posPred: t0_posPred+n3[0], cumsumnuclei[iii]:cumsumnuclei[iii+1]] = PosPred[:-(t0_posPred+n3[0]-num_possible_poly),:]
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]
        
    
        if round(tstart[iii]/FrameLen)+n2[0] != frame_num_combined:
            DataExp_combined[round(tstart[iii]/FrameLen)+n2[0]:,cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan
            DataPred_combined[round(tstart[iii]/FrameLen)+n2[0]:, cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan
            PosPred_combined[t0_posPred+n3[0]-1:, cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan
    
    
    T0_combined=np.zeros((nn_combined[1],))  ##### will contain the start of the analyzed region
    
    
    ##### correct PosPred: eliminate positions not contributing to intensity
    ##### compute T0_combined
    for i in range(nn_combined[1]): ### for all cells
        pospol =  np.where(PosPred_combined[:,i] == 1)[0] 
        times = pospol / FreqEchSimu  ; ### starting times of polymerases in seconds 
        ################################## are at times -  (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed)
        ################################## so seq start at times -  (TailleSeqMarq+TaillePostMarq)/Polym_speed)
        ind = np.where( times -  (TailleSeqMarq+TaillePostMarq)/Polym_speed > tmax_combined[i] )[0] #### positions that have no influence
        PosPred_combined[pospol[ind],i] = 0 
        max_intensity = max(DataPred_combined[:,i])
        ihit = np.where(DataPred_combined[:,i] > max_intensity/5 )[0]
        if len(ihit)!=0:
            ihi_min =  min(ihit)
            T0_combined[i] = (ihi_min-1)/FreqEchImg; #### T0_combined    
        else:
            T0_combined[i]= tmax_combined[i]

    return [DataExp_combined, DataPred_combined, PosPred_combined, T0_combined, tmax_combined, FrameLen]












def movies_combining_rawData(FilePath, xlsPath, calibrtion, extension):

        
    """ extracting files names """
    file_name_list = np.array(os.listdir(xlsPath )) # list of the data
    filesList = list(map(lambda x: x.replace(extension,'') ,file_name_list))
        
    
    ###############################################
    
    tend = np.zeros((len(file_name_list)))
    tstart = np.zeros((len(file_name_list)))
    true_FrameLen = np.zeros((len(file_name_list))) # some movies may have a different time resolution
    nfiles_l = np.ones((len(file_name_list)))
    
    
    ### this iteration is to fix the total number of nuclei and movie length
    ###d and also double checking the removal of 0 nuclei 
    for iii in range(len(file_name_list)):
        fname = FilePath + 'data_' + filesList[iii].replace('_','') + '.npz'
        ffiles = np.load(fname)
        DataExp = ffiles['DataExp']*calibrtion
        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        
        n2 = DataExp.shape
    
        fname = xlsPath +  filesList[iii] + extension
        
        ### loading time into mitosis    
        rawData = pd.read_excel(fname)
        rawData.columns = [x.lower() for x in rawData.columns]
        time_into_mitosis = rawData['time'].dropna(how='all').to_numpy()
        
        if len(time_into_mitosis) ==0:
            true_FrameLen[iii] =  3.86 
            print('!!!!! framelen pre-assigned to 3.86')
            tstart[iii] = 0
            tend[iii] = n2[0]*3.86
        else:
            true_FrameLen[iii] = np.unique(np.round(np.diff(time_into_mitosis),3))
            tstart[iii] = time_into_mitosis[0]
            tend[iii] = time_into_mitosis[-1]
        
        nfiles_l[iii] = n2[1]
    
    
    ##############
    
    
    nnuclei = round(sum(nfiles_l)) # total number of nucleis
    cumsumnuclei = np.insert(np.cumsum(nfiles_l),0,0)
    cumsumnuclei = cumsumnuclei.astype(int)
    
    
    frame_num_combined = round(min(tend/true_FrameLen))-1
    
    dataExp = np.zeros((frame_num_combined+1, nnuclei))
    tmax_combined = np.ones((nnuclei,))
    
    FrameLen = np.unique(true_FrameLen)[0]
    
    
    #### combining movies according to their time into mitosis
    for iii in range(len(filesList)):
        fname = FilePath + 'data_' + filesList[iii].replace('_','') + '.npz'
        ffiles = np.load(fname)
        DataExp = ffiles['DataExp']*calibrtion
        
        
        check_nan_rows = np.where(np.sum(np.isnan(DataExp), axis=1)== DataExp.shape[1])[0]
        DataExp = np.delete(DataExp,check_nan_rows,0)
        
        check_nan_columns = np.where(np.sum(DataExp==0, axis =0)+np.sum(np.isnan(DataExp), axis = 0)==DataExp.shape[0])[0] # this checks nuclei that are full with 0 and nan
        DataExp = np.delete(DataExp, check_nan_columns, 1)
        
        n2 = DataExp.shape 
        
    ############## concatenating nuclei according to time into mitosis
        if tstart[iii]==0:
            
            dataExp[:frame_num_combined, cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:min(frame_num_combined, DataExp.shape[0]),:]
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]
        else:
#            shape
            dataExp[int(np.round(tstart[iii]/FrameLen)): round(tstart[iii]/FrameLen)+frame_num_combined, 
                    cumsumnuclei[iii]:cumsumnuclei[iii+1]] = DataExp[:min((frame_num_combined-round(tstart[iii]/FrameLen))+1, DataExp.shape[0]),:]
            tmax_combined[cumsumnuclei[iii]:cumsumnuclei[iii+1]] = tend[iii]
    
        if np.round(tstart[iii]/FrameLen)+n2[0] != frame_num_combined:
            dataExp[int(np.round(tstart[iii]/FrameLen))+n2[0]:,cumsumnuclei[iii]:cumsumnuclei[iii+1]] = np.nan

            
    return [FrameLen, dataExp, tmax_combined, tstart]
