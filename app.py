# imports
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# helper functions
def outlier_remover(data_df):
    """
    This function can remove outliers from a pandas dataframe of numerical values for all columns
    data_df: pandas.dataframe
    output: clean dataframe
    """
    template_df = data_df.copy()
    columns = template_df.columns
    
    for column in columns:
        
        # compute the first and third quartile
        thirth_quartile, first_quartile = np.percentile(template_df.loc[:,column],[75,25])
        interquatile_range = thirth_quartile - first_quartile

        max = thirth_quartile+(1.5*interquatile_range)
        min = first_quartile-(1.5*interquatile_range)

        template_df.loc[template_df[column] < min,column] = np.nan
        template_df.loc[template_df[column] > max,column] = np.nan
        
    template_df.dropna(inplace = True)
    try:
        n = int(template_df.index[0])
        template_df.reset_index(inplace = True, drop = True)
    except:
        pass
    return template_df

def df_normalize(data_df, option = "columns"):
    """
    This function helps to normalize a dataframe that contains only numerical value.
    data_df: pandas.dataframe
    option: can be 'columns' or 'rows' to determine if normalizationj should be done across rows or columns
    
    output: normalized dataframe
    """
    
    norm_df = data_df.copy()
    
    if option == "rows":
        norm_df = norm_df.transpose()
        norm_df = (norm_df - norm_df.mean())/(norm_df.std())
        norm_df = norm_df.transpose()
        return norm_df
    
    norm_df = (norm_df - norm_df.mean())/(norm_df.std())
    return norm_df

st.set_page_config(page_title="Gene Expression Analysis", page_icon=None, layout="centered", initial_sidebar_state="auto", menu_items=None)

st.title("Gene Expression Analysis")

# st.sidebar.title("Navigation")

st.markdown("""
The purpose of this notebook is to analyze a gene expression dataset using python. We will be answering questions such as:
* Can gene expression data help us group patients or indiviauls for disease diagnosis?
* Can such data help identify relationships between patients? (differences and similarities genome wide).

### Dataset Description

The dataset for this analysis is obtained from OmicsLogic Bioinformatics course and can be downloaded [here](https://raw.githubusercontent.com/pine-bio-support/DataScience/main/Final_cell_lines_RNA-expression_FPKM_values_1000genes_with_NA.txt). The data is in tabular format and contains 1000 entries of gene expression data for 11 samples. We assume that each sample corresponds to a unique individual's data.

In this embedded table you will see the gene expression dataset, where rows are gene Ensembl ids and columns are sample IDs. Expression values are in the form of FPKM (Fragments Per Kilobase of transcript per Million mapped reads) values, which the relative expression of a transcript (proportional to the number of cDNA fragments that originate from it). 

By looking at the gene expression table, it is difficult to derive any conclusion about the data. Thus, visualization plots help us to understand the differences and derive conclusion about the data. Before, performing data visualization, it is important to process this data. Since, this type of gene expression data contains a lot of noise that we first need to address. For this task, we will clean the data by filtering (for example if you have missing data) and by transforming it’s scale into logarithm for a more “normal” distribution.

View data and transform for visualization: This type of gene expression data contains a lot of noise that we first need to address. For this task, we will clean the data by filtering (for example if you have missing data) and by transforming it’s scale into logarithm for a more “normal” distribution.


### Data Analysis
The code for this analysis can be viewed [here](https://github.com/CyrilleMesue/gene_expression_analysis/Gene_Expression_Analysis.ipynb)

Steps: 
* Display Table using pandas
* Remove Missing Values
* Transform on log scale
* Box plots
* Histograms
* Cluster Maps
* Conclusion

""")

# load data
data = pd.read_table("datasets/gene_expression analysis.txt")

st.write("Let’s look at the data and see what variables it contains. ")

st.dataframe(data)
st.write(data.shape)

#Remove ID column to kep only numeric data
data.index = data.id
data=data.drop(['id'], axis = 1) 

st.markdown("""
""")

# remove missing entries
data.dropna(inplace = True)
st.write("##### Remove missing values")
st.dataframe(data)
st.write(data.shape)

#Convert integers to floats 
datafinal = data.astype(float) 
#Perform log transformation using numpy package and show data description
log = np.log(datafinal+1) 
log.describe() 

st.markdown("##### Boxplots")
x_labels = [x.split("-")[0][0].upper() + x.split("-")[1][0].upper() for x in log.columns]
no_outlierlog = outlier_remover(log)
fig, (ax1, ax2) = plt.subplots(2, figsize=(20,12));
ax1.set_title("Box plots with outliers")
ax2.set_title("Box plots with outliers removed")
ax1.boxplot(log, labels = x_labels);
ax2.boxplot(no_outlierlog, labels = x_labels);
st.pyplot(fig)

st.markdown("""
###### Observation
* It can be clearly seen that the average gene expression values for non-malignant samples is similar amongst all non-malignant samples and different from other samples.
* The same can be said of Claudin-low samples. They all have similar gene expression mean values. 
* The difference in gene expression profile amongst similar samples is not very noticeable with box plots. However, it can be observed that the removal of outliers has helped to remove unnecessary variations in the dataset.
""")

st.markdown("##### Histogram plots")
fig, (ax1, ax2) = plt.subplots(2, figsize=(20,12));
ax1.set_title("Histogram plot with outliers")
ax2.set_title("Historgram plot with outliers removed")
ax1.hist(log, label = x_labels);
ax1.legend();
ax2.hist(no_outlierlog, label = x_labels);
ax2.legend();
st.pyplot(fig)

st.markdown("""###### Observation
* The histograms above show the number of genes with different gene expression values across samples. On the x axis is the gene expression value scaled to between 0 and 8.
* You can see that very few genes have a very high expression profile across all samples. 
* Fewer genes have a high expression level across all samples.
* Not much can be concluded from this histogram.

##### Clustermap
""")
gene_norm_df = df_normalize(no_outlierlog, option = "rows");
sns.set_context("paper", font_scale = 2);
fig = sns.clustermap(gene_norm_df, xticklabels = gene_norm_df.columns,
               figsize = (20,16)
              );
st.write("###### cluster map1")
st.pyplot(fig)

st.markdown("""
###### Observation
* The gene expression values for each gene are first standardized across samples
* This map can be used to compare samples
* As can see, the gene expression profile for each disease group is clearly different. For example, genes that are highly expressed in non-malignant samples are lowly expressed in Claudin-low samples. Also, there are more highly expressed genes in non-malignat samples than in Claudin-low samples. 
* Also, similar diseases are clustered together.
""")

sample_norm_df = df_normalize(no_outlierlog, option = "columns");
sns.set_context("paper", font_scale = 2);
fig = sns.clustermap(sample_norm_df, xticklabels = sample_norm_df.columns,
               figsize = (20,16), 
               );

st.write("###### cluster map2")
st.pyplot(fig)

st.markdown("""
#### Observation
* The gene expression values for each sample are first standardized across genes
* This map can be used to compare the expression of genes across samples
* Again, the difference between smaples is easily seen but the difference between genes is not very apparantly as there are many genes with a similar profile.


### 3. Conclusion
From the first map, it can be seen that clustering algorithm accurately clusters samples with a similar gene expression profile. This is very important for molecular diagnostics recomendations different patients with different genomic profile respond differently to diseases. Therefore, a person can be recommended a treatment that was successful on patients who had a similar gene expression profile. 

Also, for patients or individuals who have identical gene expression profiles, information from one such patient may be helpful in making prognostic decisions for another such patient.
""")

