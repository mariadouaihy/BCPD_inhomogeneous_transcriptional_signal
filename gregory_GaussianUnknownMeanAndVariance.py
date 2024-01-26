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
        """Initialize model.
        
        meanx is unknown; 
        varx is unknown;
        
        for appearance we're gonna use the same notation in Murphy 2007
        
        the parameters used to updeate meanx and varx are the following
        mu_n
        k_n
        alpha_n
        beta_n
        
        p(meanx) is computed in equation 90
        p(varx) is computed in equation 91
        p(x) is computed in equation 110
        """
        self.mu0=mu0
        self.alpha0=alpha0
        self.beta0=beta0
        self.k0=k0
        
        self.mu_params=np.array([mu0])
        self.alpha_params=np.array([alpha0])
        self.beta_params=np.array([beta0])
        self.k_params=np.array([k0])
        
        # self.mean_params = np.array([mean0])
        # self.prec_params = np.array([1/var0])
    
    def log_pred_prob(self, t, x):
        """Compute predictive probabilities \pi, i.e. the posterior predictive
        for each run length hypothesis.
        """
        """
        the implementation of the t distribution based on eq 110 in (Murphy 2007)
        is based on the link 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html
        
        df is the degree of freedom or 2*alpha
        loc is the mean 
        scale is 1/precision 
        
        in terms of the normal distribution for a degree of freedom high enouogh (at least 10),
        the t distribution with mean mu and variance sigma is almost the same as a normal distribution
        with mean mu and variance sigma
        """
        # Posterior predictive: see eq. 110 in (Murphy 2007).
        
        post_mean=self.mu_params[:t]
        post_precision=(self.alpha_params[:t]*self.k_params[:t])/(self.beta_params[:t]*(self.k_params[:t]+1))
        post_degreeOfFreedom=self.alpha_params[:t]*2
        
        # post_means = self.mean_params[:t]
        # post_stds  = np.sqrt(self.var_params[:t])
        return Tdistribution(df=post_degreeOfFreedom, loc=post_mean,scale=1/post_precision).logpdf(x)
    
    def update_params(self, t, x):
        """Upon observing a new datum x at time t, update all run length 
        hypotheses.
        """
        # See eq. 101 in (Murphy 2007).
        new_alpha_params=self.alpha_params+1/2
        # See eq. 102 in (Murphy 2007).
        new_k_params=self.k_params+1
        # See eq. 104 in (Murphy 2007).
        new_beta_params=self.beta_params+((self.k_params*(x-self.mu_params)**2)/(2*(self.k_params+1)))
        
        # See eq. 86 in (Murphy 2007). Setting the recurssion relation based on the following proof
        # \begin{equation}
        #     \mu_n=\frac{k_0\mu_0+n\bar{x}}{k_n}
        # \end{equation}
        # knowing that $\bar{x} = 1$ we have
        # \begin{equation}
        #     \begin{split}
        #         \mu_{n+1} &=\frac{k_0\mu_0+(n+1)}{k_{n+1}} \\
        #         &=\frac{k_0\mu_0+n}{k_{n+1}}+\frac{1}{k_{n+1}} \\
        #         &=\frac{k_0\mu_0+n}{k_n}\frac{k_n}{k_{n+1}}+\frac{1}{k_{n+1}} \\
        #         &=\mu_n\frac{k_n}{k_{n+1}}+\frac{1}{k_{n+1}}
        #     \end{split}
        # \end{equation} 
#        print(self.k0)
        new_mu_params= self.mu_params*(self.k_params/new_k_params)+1/new_k_params
        self.alpha_params = np.append([self.alpha0],new_alpha_params)
        self.k_params = np.append([self.k0],new_k_params)
        self.beta_params = np.append([self.beta0],new_beta_params)
        self.mu_params = np.append([self.mu0],new_mu_params)