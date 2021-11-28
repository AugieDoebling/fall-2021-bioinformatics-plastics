import pandas as pd
import numpy as np
# UNCOMMENT LINE BELOW IF SERVER DOESN'T HAVE BIO 
# !pip install biopython
from Bio import SeqIO
import re
import itertools
import cProfile
from urllib.request import urlretrieve
import time
import pickle
import sys

# get the file from the args
fasta_file = sys.argv[1]
print('using fasta file', fasta_file)

# !head ./bio_sampledata.txt

# loads fasta file into pandas dataframe
def read_fasta(filename):
  descriptor_line_pattern = re.compile(">([\w|\.]*) ([\w|\s|\:|\(|\)|\-]*) \[(.*)\]")

  with open(filename, 'r') as file:
    records = []

    id = None
    enzyme = None
    species = None
    genetic_data = ''

    for line in file.readlines():
      if line[0] == '>':
        if id is not None:
          records.append([id, enzyme, species, genetic_data])
        # either use the regex matching or just use the line as the id

        # match = descriptor_line_pattern.match(line)
        # if match is None:
        #   print('HELP', line)
        #   return
        # id = match.group(1)
        # enzyme = match.group(2)
        # species = match.group(3)

        id = line
        genetic_data = ''
      else:
        genetic_data += line.strip()
    records.append([id, enzyme, species, genetic_data, []])

  return pd.DataFrame(records, columns=['ID', 'Enzyme', 'Species', 'Genes', 'Heatmap'])

# load the data into the dataframe
data = read_fasta(fasta_file)

# load the blosum scoring matrix
url = 'https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt'
blosum62 = pd.read_csv(url, skiprows=6, delim_whitespace=True, index_col=0)
blosum62['-'] = blosum62['*']

def matrix(a, b, gap_cost=4):
    H = np.zeros((len(a) + 1, len(b) + 1), int)
    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = H[i - 1, j - 1] + blosum62.loc[a[i-1], b[j-1]] #Line causing an index error, split above to see which part is caising error
        delete = H[i - 1, j] - gap_cost
        insert = H[i, j - 1] - gap_cost
        H[i, j] = max(match, delete, insert, 0)

    return H

def traceback(H, b, b_='', old_i=0):
    # flip H to get index of **last** occurrence of H.max() with np.argmax()
    H_flip = np.flip(np.flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(H[0:i, 0:j], b, b_, i)

def smith_waterman(a, b, gap_cost=4):
    # a, b = a.upper(), b.upper()
    H = matrix(a, b, gap_cost)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_)

# cache is outside of function to allow for errors
cache = {}

def generate_heatmap(df):
  # cache each comparison - so we only do each comparison once
  # cache = {}
  for i, enzyme in df[pd.isnull(df['Heatmap'])].iterrows():
    print(f'{enzyme["ID"]} - {i} / {len(df)}')
    # init empty heatmap
    length = len(enzyme['Genes'])
    heatmap = np.zeros(length, int)
    # compare with every other enzyme
    for j, compare in df.iterrows():
      # skip same comparison
      if i == j:
        continue

      local_alignment = None
      # check if we've already done this alignment
      if (j,i) in cache:
        # print('FOUND IN CACHE')
        local_alignment = cache[(j,i)]
      elif (i,j) in cache:
        local_alignment = cache[(i,j)]
      else:
        # align comparison
        comparison_start_index = 20
        local_alignment = smith_waterman(enzyme['Genes'][comparison_start_index:], compare['Genes'][comparison_start_index:])
        # add comparison to cache
        cache[(i, j)] = local_alignment
      # update heatmap
      for hit in range(local_alignment[0], min(length-1,local_alignment[1])):
        heatmap[hit] += 1 

    # add heatmap to dataframe
    df.iloc[i]['Heatmap'] = heatmap

# generate the heatmap
generate_heatmap(data)

finish_time = time.strftime("%Y%m%d-%H%M%S")

# save the cache to a pickle - this can then be reloaded without rerunning analysis
with open(f'./cache_generated_{fasta_file}_{finish_time}.pickle', 'wb') as cache_file:
  pickle.dump(cache, cache_file, protocol=4)

# save the heatmap data to pickle
data.to_pickle(f'./heatmap_generated_{fasta_file}_{finish_time}.pickle', protocol=4)
