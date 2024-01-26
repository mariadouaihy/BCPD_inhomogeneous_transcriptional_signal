#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 14:04:35 2021

@author: mdouaihy
"""
import numpy as np
from scipy.stats import t as Tdistribution
class GaussianUnknownMeanUnknownVariance:
    
    def __init__(self, mu0, alpha0, k0, beta0):

        self.mu0=mu0
        self.alpha0=alpha0
        self.beta0=beta0
        self.k0=k0
        
        self.mu_params=np.array([mu0])
        self.alpha_params=np.array([alpha0])
        self.beta_params=np.array([beta0])
        self.k_params=np.array([k0])
        
    def log_pred_prob(self, t, x):
=
        post_mean=self.mu_params[:t]
        post_precision=(self.alpha_params[:t]*self.k_params[:t])/(self.beta_params[:t]*(self.k_params[:t]+1))
        post_degreeOfFreedom=self.alpha_params[:t]*2
        
        return Tdistribution(df=post_degreeOfFreedom, loc=post_mean,scale=1/post_precision).logpdf(x)
    
    def update_params(self, t, x):
        """Upon observing a new datum x at time t, update all run length 
        hypotheses.
        """=
        new_alpha_params=self.alpha_params+1/2
        new_k_params=self.k_params+1
        new_beta_params=self.beta_params+((self.k_params*(x-self.mu_params)**2)/(2*(self.k_params+1)))
        

        new_mu_params= self.mu_params*(self.k_params/new_k_params)+1/new_k_params
        self.alpha_params = np.append([self.alpha0],new_alpha_params)
        self.k_params = np.append([self.k0],new_k_params)
        self.beta_params = np.append([self.beta0],new_beta_params)
        self.mu_params = np.append([self.mu0],new_mu_params)
