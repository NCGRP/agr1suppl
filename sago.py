#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
from sklearn.pipeline import Pipeline
from sklearn.metrics import pairwise_distances

#parse the command line input
corpusf = sys.argv[1]
queryf = sys.argv[2]
outf = sys.argv[3]


#set up the method to convert word counts to a table
vectorizer = TfidfVectorizer(stop_words='english',strip_accents=unicode,lowercase=True,ngram_range=(1,2),use_idf=True,smooth_idf=True,sublinear_tf=True)

#set up the SVD model
#svd_model = TruncatedSVD(n_components=300,algorithm='randomized',n_iter=10,random_state=42)
svd_model = TruncatedSVD(n_components=2500,algorithm='randomized',n_iter=10)

#build a pipeline to vectorize then apply SVD
svd_transformer = Pipeline([('tfidf', vectorizer),('svd', svd_model)])

#set up the corpus
c1=pd.read_csv(corpusf,sep='@',header=None,names=['id','name','namespace','ddef']) #read in as pandas dataframe
document_corpus = c1.name.str.cat(c1.ddef, sep='. ') #extract column 2 ('name') and column 4 ('ddef') which contains deflines, as the corpus
 
#set up the query, read from a text file:
with open(queryf, "r") as f:
    query = [f.read()]
    
#apply pipeline to vectorize and SVD corpus
svd_matrix = svd_transformer.fit_transform(document_corpus)

#map query onto SVD coordinate system of corpus, this puts the words from query into the vector system previously
#generated for the document_corpus.
query_vector = svd_transformer.transform(query)

#calculate cosine difference between query and all documents of corpus
distance_matrix = pairwise_distances(query_vector, 
                                     svd_matrix, 
                                     metric='cosine', 
                                     n_jobs=-1)

df1=pd.DataFrame(distance_matrix).T #transpose distance_matrix to set up for adding to data frame
df=pd.concat([c1, df1], axis=1) #assemble the output data frame
df.columns = ['id', 'name', 'namespace', 'ddef', 'cosine'] #add/rename columns
d=df.sort_values('cosine') #sort by column 'cosine', containing the pairwise distances
d.to_csv(outf, sep='\t')