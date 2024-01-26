#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:13:06 2022

@author: mdouaihy
"""

import matplotlib.pyplot as plt
import numpy as np
from   matplotlib.colors import LogNorm
import sys
from copy import copy
import copy


def retriving_time(R,p):
    time_where_change_point_occured=[]
    probability_at_each_run_length=[]
    for i in range(1,len(R)):
        indx_data_i=np.where(R[:,i]>=p)
        time_where_change_point_occured.append(indx_data_i[0])    
        probability_at_each_run_length.append(i+R[indx_data_i,i])
        
    
    test_time = copy.deepcopy(time_where_change_point_occured)
    lines = []
    index_starting_point = []
    
    for i in range(len(test_time)):
        curent_initial_points = test_time[i]
        nbr_points_in_curent_position = len(curent_initial_points)
        testing = test_time[i+1:]
    
        
        for point_i in range(nbr_points_in_curent_position):
            curent_line=[]
            if curent_initial_points[point_i]!=-1:
                
                curent_initial_pt = curent_initial_points[point_i]
                curent_line.append(curent_initial_pt)
                test_time[i][point_i]=-1
                for j in range(len(testing)):
                    next_possibilities = testing[j] - curent_initial_pt
                    indx_in_same_line = np.where(next_possibilities==j+1)[0]
                    if len(indx_in_same_line )>0:
                        # print(np.float(test_time[i+1+j][indx_in_same_line]))
                        curent_line.append(int(test_time[i+1+j][indx_in_same_line]))
                        test_time[i+1+j][indx_in_same_line]=-1
              
                curent_line=np.array(curent_line)-1
                lines.append(curent_line.tolist())
    
                index_starting_point.append(probability_at_each_run_length[i][0,point_i])
    return lines, index_starting_point




def retriving_time_100723(R,p):

    time_where_change_point_occured=[]
    probability_at_each_run_length=[]
    for i in range(len(R)):
        indx_data_i=np.where((R[i,:])>=p)[0]
    
        time_where_change_point_occured.append(indx_data_i)    
        probability_at_each_run_length.append(R[i, indx_data_i])
    
    
    test_time = copy.deepcopy(time_where_change_point_occured)
    lines = []
    proba_line = []
    x_axis_lines = []
    
    
    for i in range(len(test_time)):
        curent_initial_points = test_time[i]
        nbr_points_in_curent_position = len(curent_initial_points)
        testing = test_time[i+1:]
    
        for point_i in range(nbr_points_in_curent_position):
            curent_line=[]
            current_prob = []
            x_current_line = []
            if curent_initial_points[point_i]!=-1:
                
                curent_initial_pt = curent_initial_points[point_i]
                curent_line.append(curent_initial_pt)
                x_current_line.append(i)
                current_prob.append(probability_at_each_run_length[point_i][0])
                test_time[i][point_i]=-1
                for j in range(len(testing)):
                    next_possibilities = testing[j] - curent_initial_pt
                    indx_in_same_line = np.where(next_possibilities==j+1)[0]
                    if len(indx_in_same_line )>0:
                        # print(np.float(test_time[i+1+j][indx_in_same_line]))
                        curent_line.append(int(test_time[i+1+j][indx_in_same_line]))
                        current_prob.append(probability_at_each_run_length[i+1+j][indx_in_same_line][0])
                        
                        test_time[i+1+j][indx_in_same_line]=-1
              
                curent_line=np.array(curent_line)-1
                
                lines.append(curent_line.tolist())
                proba_line.append(current_prob)
                x_axis_lines.append(np.unique(np.array(x_current_line)))
    return x_axis_lines, proba_line

# plots to verify
# plt.figure(1)
# for i in range(len(R)-1):
#     plt.scatter(86-time_where_change_point_occured[i],87-np.rot90(probability_at_each_run_length[i]),c='b', s=10)
# plt.imshow(R, aspect='auto', cmap='gray_r')#,norm=LogNorm(vmin=0.00001, vmax=1))

# plt.figure(30)
# for i in range(len(R)-1):
#     plt.scatter(time_where_change_point_occured[i],probability_at_each_run_length[i],c='b', s=10)

# # sys.exit()
# plt.figure(2)
# plt.imshow(np.rot90(R), aspect='auto', cmap='gray_r', 
#                 norm=LogNorm(vmin=0.00001, vmax=1))

# h = plt.figure(10)
# for i in range(len(lines)):
#     plt.plot(lines[i],index_starting_point[i]+np.arange(len(lines[i])))
# for i in range(len(R)-1):
#     plt.scatter(time_where_change_point_occured[i],probability_at_each_run_length[i],c='b', s=2)
