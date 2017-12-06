import numpy as np
import pandas as pd
import sys
import argparse
from sklearn import linear_model
from sklearn.model_selection import train_test_split
import os
from sklearn.feature_selection import SelectKBest, f_regression, VarianceThreshold
from sklearn.metrics import mean_squared_error
import math

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


features = list(set(adipose.columns).intersection(set(muscle.columns)).intersection(set(thyroid.columns)).intersection(set(whole_blood.columns)))

adipose = adipose[features]
muscle = muscle[features]
thyroid = thyroid[features]
whole_blood = whole_blood[features]

# Feature selection--solely based on training data and labels
def select_features(data_tup):
    train_x, test_x, train_y, test_y = data_tup
    selector = SelectKBest(f_regression, k=5000)
    selector.fit(train_x, train_y)
    col_indices = selector.get_support(indices=True)
    return train_x.iloc[:, col_indices], test_x.iloc[:, col_indices], train_y, test_y

# Training-testing splits for each tissue
adipose_train, adipose_test, adipose_labels_train, adipose_labels_test = select_features(train_test_split(adipose, adipose_labels, test_size=0.30, random_state=1))
muscle_train, muscle_test, muscle_labels_train, muscle_labels_test = select_features(train_test_split(muscle, muscle_labels, test_size=0.30, random_state=1))
thyroid_train, thyroid_test, thyroid_labels_train, thyroid_labels_test = select_features(train_test_split(thyroid, thyroid_labels, test_size=0.30, random_state=1))
whole_blood_train, whole_blood_test, whole_blood_labels_train, whole_blood_labels_test = select_features(train_test_split(whole_blood, whole_blood_labels, test_size=0.30, random_state=1))

features = (list(set(adipose_train.columns).intersection(set(muscle_train.columns)).intersection(set(thyroid_train.columns)).intersection(set(whole_blood_train.columns))))

adipose_train, adipose_test = adipose_train[features], adipose_test[features]
muscle_train, muscle_test = muscle_train[features], muscle_test[features]
thyroid_train, thyroid_test = thyroid_train[features], thyroid_test[features]
whole_blood_train, whole_blood_test = whole_blood_train[features], whole_blood_test[features]

adipose_train_out = pd.concat([adipose_train, adipose_labels_train], axis=1)
adipose_test_out = pd.concat([adipose_test, adipose_labels_test], axis=1)
muscle_train_out = pd.concat([muscle_train, muscle_labels_train], axis=1)
muscle_test_out = pd.concat([muscle_test, muscle_labels_test], axis=1)
thyroid_train_out = pd.concat([thyroid_train, thyroid_labels_train], axis=1)
thyroid_test_out = pd.concat([thyroid_test, thyroid_labels_test], axis=1)
whole_blood_train_out = pd.concat([whole_blood_train, whole_blood_labels_train], axis=1)
whole_blood_test_out = pd.concat([whole_blood_test, whole_blood_labels_test], axis=1)


adipose_train_out.to_csv("processed_data/adipose_train.csv")
adipose_test_out.to_csv("processed_data/adipose_test.csv")
muscle_train_out.to_csv("processed_data/muscle_train.csv")
muscle_test_out.to_csv("processed_data/muscle_test.csv")
thyroid_train_out.to_csv("processed_data/thyroid_train.csv")
thyroid_test_out.to_csv("processed_data/thyroid_test.csv")
whole_blood_train_out.to_csv("processed_data/whole_blood_train.csv")
whole_blood_test_out.to_csv("processed_data/whole_blood_test.csv")

