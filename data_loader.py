import numpy as np
import pandas as pd
import sys
import argparse
import os

def get_data():
	cur_path = os.getcwd()
	parser = argparse.ArgumentParser(description='Process display arguments')
	parser.add_argument("-f", "--jupyter-json")
	parser.add_argument("-adipose-file", "--adipose-file", default=cur_path+"/data/Adipose_Subcutaneous_Analysis.v6p.normalized.expression.bed")
	parser.add_argument("-muscle-file", "--muscle-file", default=cur_path+"/data/Muscle_Skeletal_Analysis.v6p.normalized.expression.bed")
	parser.add_argument("-thyroid-file", "--thyroid-file", default=cur_path+"/data/Thyroid_Analysis.v6p.normalized.expression.bed")
	parser.add_argument("-whole-blood-file", "--whole-blood-file", default=cur_path+"/data/Whole_Blood_Analysis.v6p.normalized.expression.bed")
	parser.add_argument("-label-file", "--label-file", default=cur_path+"/data/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt")
	args = parser.parse_args()

	# Matrix generation
	def gen_matrix(file_path):
	    df = pd.read_csv(file_path, header=0, sep='\t', dtype=str)
	    df.drop(["#chr", "start", "end"], axis=1, inplace=True)
	    df.set_index("gene_id", inplace=True)
	    df = df.transpose()
	    df.columns.name = None
	    return df.apply(pd.to_numeric)
	    
	adipose = gen_matrix(args.adipose_file)
	muscle = gen_matrix(args.muscle_file)
	thyroid = gen_matrix(args.thyroid_file)
	whole_blood = gen_matrix(args.whole_blood_file)

	labels = pd.read_csv(args.label_file, header=0, sep='\t', dtype=str)
	labels = labels.set_index("SUBJID").drop(["GENDER", "DTHHRDY"], axis=1)
	labels.index.name = None
	labels["AGE"] = labels["AGE"].apply(lambda x: float((float(x.split("-")[1]) +  float(x.split("-")[0])) / 2))

	adipose_labels = labels.loc[adipose.index].copy()
	muscle_labels = labels.loc[muscle.index].copy()
	thyroid_labels = labels.loc[thyroid.index].copy()
	whole_blood_labels = labels.loc[whole_blood.index].copy()

	return {"adipose_x": adipose, "adipose_y": adipose_labels, "muscle_x": muscle, "muscle_y": muscle_labels, "thyroid_x": thyroid, "thyroid_y": thyroid_labels, "whole_blood_x": whole_blood, "whole_blood_y": whole_blood_labels}

