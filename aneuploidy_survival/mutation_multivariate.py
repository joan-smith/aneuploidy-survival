# coding: utf-8
import sys
import pandas as pd
import sys

import biomarker_survival as surv

import aneuploidy_survival

def multivariate_mutation():
  df = surv.prep_mutation_data_alone('data/mutation.mc3.v0.2.8.tsv')
  types = clinical['type'].unique()

  aneuploidy_data, clinical_data = aneuploidy_survival.prep_data('data', '.')
  print aneuploidy_data

  output = pd.DataFrame()
  for t in types:
    print t
    clinical_for_type = clinical[clinical['type'] == t]
    clinical_for_type = clinical_for_type.set_index('bcr_patient_barcode')

    muts = surv.prep_mutation_data(df, clinical_for_type, t) #note mutation data doesn't include type, so process all
    clinical_and_muts = muts.join(clinical_for_type)
    endpoints = ['OS', 'PFI', 'DFI', 'DSS']
    for e in endpoints:
      print e
      if '\'TP53' not in clinical_and_muts.columns:
        continue
      d = clinical_and_muts[[e + '.time', e, '\'TP53', ]]
      d = d.dropna(how='any')
      if len(d) > 10:
        out = surv.do_cox(d[e + '.time'], d[e], d['\'TP53'])
        out['endpoint'] = e
        out['type'] = t
        output = output.append(out, ignore_index=True)
  output.set_index(['type', 'endpoint'], inplace=True)
  output.to_csv('data/p53_mutation_zscores.csv', index_label=['Type', 'Endpoint'])
  print output


multivariate_mutation()
