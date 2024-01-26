#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:56:37 2022

@author: mdouaihy
"""

def MLE_gamma(x):
    mean = x.mean()
    var  = x.var()

    alpha = (mean**2)/var
    beta  = alpha / mean
    return alpha,beta
