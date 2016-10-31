##
 #  Project: HAMPtOn
 #
 #  File: mlr_learn.py
 #  Created: Oct 10, 2016
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

## multiple linear regression learning (two variables: ntotcells (x1) and nverttotcells (x2))
## model: y = X b
## where, y = (y1, y2, ..., yn)
##        X = (X1, X2, ..., Xn) = (x11, x21, ..., xn1)
##                                (x12, x22, ..., xn2)
##        b = (b1, b2)
## find b using ordinary least-squares (olr)

## output set:
##
## Version: hs
## Residual: [ 83934818.90832797]
## Solution: [ 7.05857794  0.11344668]
## Version: new7
## Residual: [ 81898385.86956327]
## Solution: [ 6.98452379  0.11773018]
## Version: new8
## Residual: [ 87998168.81527595]
## Solution: [ 7.35814467  0.10984191]
## Version: new9
## Residual: [  1.02529304e+08]
## Solution: [ 7.2506767   0.11230974]
##
## Final Solution: [ 7.16298077  0.11333213]
##
## b_ntotcells      = 7.16298077
## b_nverttotcells  = 0.11333213

import sys
import argparse
import numpy as np
import pandas as pd
import itertools as it

from numpy import linalg as la
from sets import Set


class MLRModel:
  
  def __init__(self, nvars):
    if nvars != 2: print 'error: only 2 variables are supported'
    self.num_vars = nvars
    self.X = [[]]
    self.y = []
    self.b = []
    
  def fit(self, X1, X2, y):
    if len(X1) != len(X2): print 'error: the two variable dimensions should be equal'
    self.X = np.matrix([ X1, X2 ]).T
    self.y = y
    # print 'X:', self.X
    # print 'y:', self.y
    b, residual, rank, s = la.lstsq(self.X, self.y)
    print 'Residual:', residual
    return b


def load_raw_data(filename, vars):
  columns = [ '%', 'excl', 'incl', 'ncall', 'nsubr', 'inclpcall' ]  
  ff = open(filename)
  head1 = ff.readline()
  ff.close()
  head1 = head1.strip().split()[3:]
  l = lambda x: [ (x, c) for c in columns ]
  head1 = list(it.chain(*map(l, head1)))
  df = pd.read_table(filename, sep='\s+', skiprows=1, index_col=[0, 1, 2], na_values='-')
  df.columns = pd.MultiIndex.from_tuples(head1)
  cols = list(set.intersection(set(vars), set(list(df.columns.levels[0]))))
  rvals = df[cols].iloc[df[cols].index.get_level_values('node') >= 0].xs('excl', level=1, axis=1)
  rvals = rvals.sum(axis=1).dropna(axis=0, how='any')
  return rvals

def load_data(infile, v1, v2=None):
  data = pd.read_table(infile, sep='\s+')
  v1_data = data[v1]
  if v2 is not None:
    v2_data = data[v2]
    return v1_data, v2_data
  else: return v1_data

def read_function_list(f):
  funclist = []
  ff = open(f)
  if not ff: return funclist
  while True:
    line = ff.readline()
    if not line: break
    line = line.strip().replace(' ', '$')
    if not line.startswith('#'): funclist.append(line)
  return funclist

def notcontains_filter(slist, pattern):
  flist = []
  for s in slist:
    if pattern not in s: flist.append(s)
  return flist

##
## the main part
##

dataset = 'oEC60to30'
stat_prefix = dataset + '/stats.'
time_prefix = dataset + '/tau.titan.'
stat_file = 'graph.info.part.64.stats'
time_file = 'extracted/N2d1.32.static.GET_TIME_OF_DAY.dat'
# versions = [ 'base', 'hs', 'new7', 'new8', 'new9' ]
versions = [ 'hs', 'new7', 'new8', 'new9' ]

if __name__ == '__main__':
  bsum = None
  for version in versions:
    print 'Version:', version
    infile1 = stat_prefix+version+'/'+stat_file
    infile2 = 'oEC60to30/tau.titan.'+version+'/'+time_file
    X1, X2 = load_data(infile1, v1='ntotcells', v2='nverttotcells')
    vars = read_function_list('compute.functions')
    vars = notcontains_filter(vars, '=>')
    y = load_raw_data(infile2, vars=vars)
    mlr = MLRModel(2)
    b = mlr.fit(X1=X1, X2=X2, y=y)
    if bsum is None: bsum = b
    else: bsum += b
    print 'Solution:', b
  print 'Final Solution:', bsum / len(versions)
