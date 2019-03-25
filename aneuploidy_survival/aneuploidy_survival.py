#!/usr/bin/env python
# encoding: utf-8
'''
aneuploidy_survival.py

Created by Joan Smith on 2018-12-9.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob

import biomarker_survival as surv


def get_options():
  parser = argparse.ArgumentParser(description='Aneuploidy survival')
  parser.add_argument('-i', action='store', dest='input_dir')
  parser.add_argument('-o', action='store', dest='outdir', default=None)

  ns = parser.parse_args()

  if ns.outdir == None:
    ns.outdir = ns.input_dir

  return (ns.input_dir, ns.outdir)

def maybe_mkdir(path):
  if not os.path.exists(path):
    os.mkdir(path)

def endpoint_data(clinical_data, outdir):
  events = clinical_data.groupby('type')['OS', 'DSS', 'DFI', 'PFI'].sum()
  counts = clinical_data.groupby('type')['OS.time', 'DSS.time', 'DFI.time', 'PFI.time'].count()
  endpoint_data = counts.join(events)
  endpoint_data = endpoint_data.reindex_axis(sorted(endpoint_data.columns), axis=1)
  endpoint_data.to_csv(os.path.join(outdir, 'endpoint_data.csv'))

  return clinical_data[['bcr_patient_barcode', 'type'] + list(endpoint_data.columns)]

def pancan(name, data_dir, pancan_dir):
  endpoints = ['OS', 'DSS', 'DFI', 'PFI']
  cancer_type = lambda f: os.path.basename(os.path.dirname(f))
  for end in endpoints:
    files = glob.glob(os.path.join(data_dir, '*',  '*' + end + '*.csv'))
    pancan_data = {cancer_type(f): pd.read_csv(f, usecols=['var', 'z'], index_col=0)['z'] for f in files}
    pancan_out = pd.DataFrame(pancan_data)
    pancan_out['stouffer'] = surv.stouffer_unweighted(pancan_out)
    pancan_out.to_csv(os.path.join(pancan_dir, name + '_pancan_' + end + '.csv'), index_label='var')

def clean_aneuploidy_data(aneuploidy_data, ctype):
  cancer_type_aneuploidy = aneuploidy_data[aneuploidy_data['Type'] == ctype].copy()
  cancer_type_aneuploidy = surv.maybe_clear_non_01s(cancer_type_aneuploidy, 'Sample', ctype)
  cancer_type_aneuploidy = surv.add_identifier_column(cancer_type_aneuploidy, 'Sample')
  cancer_type_aneuploidy.set_index('identifier', inplace=True)
  cancer_type_aneuploidy = cancer_type_aneuploidy.drop(['Sample', 'Type'], axis=1)
  return cancer_type_aneuploidy


def zscores(name, aneuploidy_data, patient_clinical, prep_data_fn, column_subset, indir, outdir):
  for ctype in aneuploidy_data['Type'].unique():
    print 'Cancer Type:', ctype
    maybe_mkdir(os.path.join(outdir, ctype))

    # Clean up data, remove non-01s where appropriate
    cancer_type_aneuploidy = clean_aneuploidy_data(cancer_type_aneuploidy, ctype)
    cancer_type_aneuploidy = cancer_type_aneuploidy[column_subset(cancer_type_aneuploidy)]

    cancer_type_clinical = patient_clinical[patient_clinical['type'] == ctype].set_index('bcr_patient_barcode')
    cancer_type_clinical = cancer_type_clinical.drop('type', axis=1)
    cancer_type_aneuploidy_clin = cancer_type_aneuploidy.join(cancer_type_clinical, how='outer')

    endpoint_cols = [('OS', 'OS.time'), ('DSS', 'DSS.time'), ('DFI', 'DFI.time'), ('PFI', 'PFI.time')]

    for c, t in endpoint_cols:
      print '  Endpoint:', c
      output = {}
      tmp = pd.DataFrame()
      for var in cancer_type_aneuploidy.columns:
        data = cancer_type_aneuploidy_clin[[t, c, var]].dropna(how='any')
        data = prep_data_fn(data, var)
        if len(data) > 10:
          output[str(var)] = surv.do_cox(data[t], data[c], data[var])
          output[str(var)]['censor count'] = data[c].sum()
      if len(output) > 0:
        pd.DataFrame(output).T.to_csv(os.path.join(outdir, ctype, name + '_' + c + '_and_' + t  + '.csv'), index_label='var')

def prep_data(indir, outdir):
  clinical_data = pd.read_excel(os.path.join(indir, 'TCGA_clinical_data.xlsx'), sheet_name=0)
  patient_clinical = endpoint_data(clinical_data, outdir)

  aneuploidy_data = pd.read_excel(os.path.join(indir, 'TCGA_aneuploidy_data.xlsx'), sheet_name=0,
                                  header=1)
  aneuploidy_data['Type'] = aneuploidy_data['Type'].str.strip()
  return aneuploidy_data, patient_clinical


def prep_analysis_dirs(name, outdir):
  print name
  data_dir = os.path.join(outdir, name)
  pancan_dir = os.path.join(outdir, name + '_pancan')
  maybe_mkdir(data_dir)
  maybe_mkdir(pancan_dir)
  return data_dir, pancan_dir


def main():
  indir, outdir = get_options()

  aneuploidy_data, clinical_data = prep_data(indir, outdir)

  # Analysis:  everything univariate, straightforward.
  name = 'zscores'
  data_dir, pancan_dir = prep_analysis_dirs(name, outdir)
  prep_data_fn = lambda d, var: d
  all_columns = lambda d: d.columns
  zscores(name, aneuploidy_data, clinical_data, prep_data_fn, all_columns, indir, data_dir)
  pancan(name, data_dir, data_dir)

  only_chromosome_scores_columns = lambda d: [x for x in d.columns if (type(x) == int or x[0].isdigit())]
  # Analysis: Bin loss and neutral together
  name = '-1s_are_0'
  data_dir, pancan_dir = prep_analysis_dirs(name, outdir)
  prep_data_fn = lambda d, var: d.replace({var: -1}, 0)
  zscores(name, aneuploidy_data, clinical_data, prep_data_fn, only_chromosome_scores_columns, indir, data_dir)
  pancan(name, data_dir, pancan_dir)

  # Analysis: Bin gain and neutral together
  name = '1s_are_0'
  data_dir, pancan_dir = prep_analysis_dirs(name, outdir)
  prep_data_fn = lambda d, var: d.replace({var: 1}, 0)
  zscores(name, aneuploidy_data, clinical_data, prep_data_fn, only_chromosome_scores_columns, indir, data_dir)
  pancan(name, data_dir, pancan_dir)

  aneuploidy_score = lambda d: ['AneuploidyScore(AS)']
  # Analysis: bulk aneuploidy
  name = 'as_20p_v_80p'
  data_dir, pancan_dir = prep_analysis_dirs(name, outdir)
  as_20p_v_80p = lambda data, var: as_percentile(20, 80, data, var)
  zscores(name, aneuploidy_data, clinical_data, as_20p_v_80p, aneuploidy_score, indir, data_dir)
  pancan(name, data_dir, pancan_dir)

def as_percentile(lower, upper, data, var):
  if len(data) == 0:
    return data
  lower = np.percentile(data[var], lower)
  upper = np.percentile(data[var], upper)
  data = data.drop(data[data[var].between(lower, upper, inclusive=False)].index).copy()
  data.loc[data[var] <= lower, var] = 0
  data.loc[data[var] >= upper, var] = 1
  return data

if __name__ == "__main__":
  main()
