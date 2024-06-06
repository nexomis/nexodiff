#!/usr/bin/env python3

import h5py
import numpy as np
import os
import sys
from scipy.stats import gamma, norm
import yaml

# what is read from kallisto files :

# est_counts 64-bits floating-points
# aux/ids strings
# aux/eff_lengths 64-bits floating-points

def generate_file(ids,est_counts,eff_lengths,filename):
  """
  Generate hdf5 file that mimicks kallisto outputs
  This function simply take ids, effective lenths and estimated counts as
  inputs and create an hdf5 file that mimicks kallisto output

  Parameters
  ----------
  ids : str[]
      transcript ids
  est_counts : float[]
      estimated counts
  eff_lengths : float[]
      effective length
  filename : str
      filename to output

  Returns
  -------
  None

  """
  est_counts = np.asarray(est_counts, dtype="float64")
  eff_lengths = np.asarray(eff_lengths, dtype="float64")
  ids = np.string_(ids)

  if eff_lengths.shape != est_counts.shape or ids.shape != est_counts.shape:
      raise ValueError("Passed arrays are not of the same shapes")

  h5file = h5py.File(filename, "w")

  d_est_counts = h5file.create_dataset("est_counts", data=est_counts)
  aux = h5file.create_group("aux")
  d_eff_lengths = aux.create_dataset("eff_lengths", data=eff_lengths)
  d_ids = aux.create_dataset("ids",data=ids)

  h5file.close()

def generate_annotation(annotation_file, id_mapping_file):

  """
  Generate a transcriptome annotation file

  Parameters
  ----------
  annotation_file : str
    filename to write the transcriptome annotation

  id_mapping_file : str
    filename to write the transcriptome annotation
  Returns
  -------
  str[]: The list of transcript ids

  """

  with open(annotation_file, "w") as ann_file:
    with open(id_mapping_file, "w") as id_file:

      # generate 5000 genes with 1 transcripts for human of size 1000
      taxid="9606"
      taxname="human"
      txids=[]
      id_file.write("\t".join(["Entry","Status","Protein names","Gene names","Cross-reference (GeneID)","Annotation"]) + os.linesep)
      for i in range(1,501):
        geneid=str(i)
        genename="gene"+geneid
        uid="UN"+geneid
        protname="long name of protein for gene "+geneid
        genenames=geneid+"_1 "+geneid+"_2"
        
        txid="tx_g"+geneid
        txids.append(txid)
        if i<5:
          txtype="rRNA"
        elif i<400:
          txtype="mRNA"
        elif i<450:
          txtype="ncRNA"
        else:
          txtype="miscRNA"
        ann_file.write( " ".join([txid,geneid,txtype,genename,taxid,taxname]) + os.linesep)
        id_file.write("\t".join([uid,"reviewed",protname,genenames,geneid+";","5 out of 5"]) + os.linesep)

      taxid="2697049"
      taxname="SARS2"
      txtype="mRNA"
      for i in range(501,516):
        geneid=str(i)
        genename="gene"+geneid
        txid="tx_g"+geneid
        txids.append(txid)
        ann_file.write( " ".join([txid,geneid,txtype,genename,taxid,taxname]) + os.linesep)
  return(txids)

def generate_data(txids,folder):
  """
  Generate data 

  Parameters
  ----------
  folder : str
    folder to write in the h5 files per sample
    
  Returns
  -------
  dict: simulated estimated count per sample

  """
  n_genes = len(txids)
  np.random.seed(seed=1)
  samplesheet=dict()
  samplesheet["variables"]=["donnorID","TissueType"]
  design=dict()
  design["import"]= { 
    "batch1": {
      "A": {
        "b1_A1": ["D1","lung"],
        "b1_A2": ["D2","lung"],
        "b1_A3": ["D3","lung"]
      },
      "B":{
        "b1_B1": ["D1","lung"],
        "b1_B2": ["D2","lung"],
        "b1_B3": ["D3","lung"]
      },
      "C":{
        "b1_C1": ["D1","lung"],
        "b1_C2": ["D2","lung"],
        "b1_C3": ["D3","lung"]
      }
    },
    "batch2":{
      "A":{
        "b2_A1": ["D3","liver"],
        "b2_A2": ["D4","liver"],
        "b2_A3": ["D5","liver"]
      },
      "B":{
        "b2_B1": ["D3","liver"],
        "b2_B2": ["D4","liver"],
        "b2_B3": ["D5","liver"]
      },
      "C":{
        "b2_C1": ["D3","liver"],
        "b2_C2": ["D4","liver"],
        "b2_C3": ["D5","liver"]
      }
    }
  }

  samplesheet["design"]=[design]

  samplesheet["g_labels"]={
    "A": "condition A",
    "B": "condition B",
    "C": "control"
  }

  samplesheet["b_labels"]={
    "batch1": "Experiment batch 1",
    "batch2": "Experiment batch 2"
  }

  samplesheet["ctrlGroups"]={
    "batch1": "C",
    "batch2": "C"
  }

  samplesheet["paired_py"]={
    "batch1": "donnorID"
  }

  fc_sd = 1

  # base expr with gamma 2 2
  base = np.array(gamma.rvs(a=2, loc=2, size=n_genes))
  fc_conds=dict()
  fc_conds["C"] = norm.rvs(loc=0,scale=fc_sd,size=n_genes)
  fc_conds["A"] = fc_conds["C"].copy()
  fc_conds["B"] = fc_conds["C"].copy()
  fc_conds["A"][:10] = norm.rvs(loc=2,scale=fc_sd,size=10)
  fc_conds["A"][10:20] = norm.rvs(loc=-2,scale=fc_sd,size=10)
  fc_conds["B"][20:30] = norm.rvs(loc=2,scale=fc_sd,size=10)
  fc_conds["B"][30:40] = norm.rvs(loc=-2,scale=fc_sd,size=10)
  fc_tissue=dict()
  fc_tissue["lung"] = norm.rvs(loc=0,scale=fc_sd,size=n_genes)
  fc_tissue["liver"]  = fc_tissue["lung"].copy()
  fc_tissue["lung"][40:50] = norm.rvs(loc=2,scale=fc_sd,size=10)
  fc_tissue["lung"][50:60] = norm.rvs(loc=-2,scale=fc_sd,size=10)
  fc_tissue["liver"][60:70] = norm.rvs(loc=2,scale=fc_sd,size=10)
  fc_tissue["liver"][70:80] = norm.rvs(loc=-2,scale=fc_sd,size=10)
  fc_donnor=dict()
  fc_donnor["D1"]=norm.rvs(loc=0,scale=fc_sd,size=n_genes)
  for d in ["D2","D3","D4","D5"]:
    fc_donnor[d]=fc_donnor["D1"].copy()
  start=80
  for d in ["D1","D2","D3","D4","D5"]:
    fc_donnor[d][start:start+10]= norm.rvs(loc=2,scale=fc_sd,size=10)
    start=start+10
    fc_donnor[d][start:start+10]= norm.rvs(loc=-2,scale=fc_sd,size=10)
    start=start+10
  
  samples_expr=dict()

  # batch 1 lung D1 D2 D3
  donnors=["D1","D2","D3"]
  for i in range(1,4):
    for cond in ["A","B","C"]:
      sample_name = "b1_"+cond+str(i)
      samples_expr[sample_name]=base.copy()
      samples_expr[sample_name]=np.multiply(samples_expr[sample_name],
        np.power([2]*n_genes,fc_conds[cond]))
      samples_expr[sample_name]=np.multiply(samples_expr[sample_name],
        np.power([2]*n_genes,fc_tissue["lung"]))
      samples_expr[sample_name]=np.multiply(samples_expr[sample_name],
        np.power([2]*n_genes,fc_donnor[donnors[i-1]]))
  # batch 1 lung D1 D2 D3
  donnors=["D3","D4","D5"]
  for i in range(1,4):
    for cond in ["A","B","C"]:
      sample_name = "b2_"+cond+str(i)
      samples_expr[sample_name]=base.copy()
      samples_expr[sample_name]=np.multiply(samples_expr[sample_name],
        np.power([2]*n_genes,fc_conds[cond]))
      samples_expr[sample_name]=np.multiply(samples_expr[sample_name],
        np.power([2]*n_genes,fc_tissue["liver"]))
      samples_expr[sample_name]=np.multiply(samples_expr[sample_name],
        np.power([2]*n_genes,fc_donnor[donnors[i-1]]))

  est_count=dict()

  for sample in samples_expr.keys():
    N=norm.rvs(loc=10000000,scale=100000,size=1)
    est_count[sample]=N*samples_expr[sample]/samples_expr[sample].sum()
    generate_file(txids,est_count[sample],[1000]*n_genes,folder+"/"+sample+".h5")

  with open(folder+"/samples.yml","w") as file:
    yaml.dump(samplesheet, file)

  samplesheet["design"][0]["path"] = "/fake/path/to/test"

  with open(folder+"/samples_with_fake_path.yml","w") as file:
    yaml.dump(samplesheet, file)

  return(est_count)
  
def main():
  os.makedirs("inst/extdata/sim_data", mode = 0o777, exist_ok = True)
  txids=generate_annotation("inst/extdata/sim_data/annotation.txt","inst/extdata/sim_data/id_mapping.tab")
  generate_data(txids,"inst/extdata/sim_data")
  sys.exit()

if __name__ == "__main__":
    main()
