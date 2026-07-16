## Overview
This repository contains the source code for the paper **"Microbiome-Based Classification of Soil Conditions Using Machine Learning and Explainable AI"**. This project implements a Machine Learning pipeline using Linear SVM and Random Forest to classify high-dimensional microbiome data. It leverages Explainable AI (XAI) through SHAP and ShapG to quantify feature contributions and identify the most influential genera across multiple sequencing platforms.

## Dataset Overview
We frist used the raw (unnormalized) version of 3 datasets (Aviti_raw_data.txt, Miseq_raw_data.txt, Novaseq_raw_data.txt). Second, we used the percentage-based version of the datasets (Aviti_perc_data.txt, Miseq_perc_data.txt, Novaseq_perc_data.txt). However, the final version of the datasets were:
- Aviti_DESeq2_normalized_data.txt
- Miseq_DESeq2_normalized_data.txt
- Novaseq_DESeq2_normalized_data.txt
These datasets contain 36 samples each and are analyzed independently rather than being mixed. Each time one dataset is considered, split into train and test, and the rest of analysis is done.

- **Train Dataset:** we have 3 datasets. 27 samples (3/4) of each was used as train set. 
- **Test Datasets:**  we have 3 datasets. 9 samples (1/4) of each was used as train set. 
  

## Implementation
- **Notebook:** `Soil_Classification-shapg(supercomputer).ipynb`

## Environment
**Python version:** 3.11.5
**and Core libraries include:**
- NumPy  
- Pandas  
- Scikit‑learn (Version: 1.5.2)  
- SHAP & shapG
- Matplotlib  
- Seaborn

## Contact information
- For data preprocessing or R/DESeq2 questions, contact Samuele.
- For Python code, machine learning models, or SHAP/ShapG questions, contact Fatemeh.

## Citation
If you use this repository, please cite:
**"Microbiome-Based Classification of Soil Conditions using Machine Learning and Explainable AI"**
(BibTeX entry will be added once available).
