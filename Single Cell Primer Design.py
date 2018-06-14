
# coding: utf-8

# In[24]:


import os
import sys
import primer3
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq


# #### Collect consensus regions longer then 19 bp

# In[25]:



align = AlignIO.read("DV4_all_cleanup.fasta","fasta") # change the filename and the file format to match you alignment result.


summary_align = AlignInfo.SummaryInfo(align)
consensus = summary_align.dumb_consensus(threshold=0.95,ambiguous='X') # 0.95  as the threshold of considered as consensus.

a = consensus.split('X')

pr_list = []
for i in a:
	if len(i) > 17:
		pr_list.append(i)


# #### Tools

# In[26]:


# OH for overhang of the primer
OH = 'AAGCAGTGGTATCAACGCAGAGT'
# TSO for template switching oligo
TSO = 'AAGCAGTGGTATCAACGCAGAGTACATGGG'
# Reverse complement tool
def rc (my_sequence):
	my_dictionary = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return "".join([my_dictionary[base] for base in reversed(my_sequence)])


# #### Chop up the consensus regions and collect fragments between 18 and 24 bp as primers.

# In[27]:


full_list = []
for i in pr_list:
	for j in range(0,len(i)):
		for k in range(0,len(i)):
			if 17 < k-j < 25:
				full_list.append(i[j:k])


# #### Reverse complement the primers and test for criterias for TSO compatibility by Primer3. Select only the primers with Tm > 50. Remove any primers with "CC" or "TTT".

# In[28]:


rc_pr_list = []
for i in full_list:
	
	if ("CC" not in rc(i) and "TTT" not in rc(i) and primer3.calcTm(rc(i)) > 50):
		l = OH + rc(i)
		rc_pr_list.append(l)


# #### Select the primers with least tendency to form heterodimers with TSO. dG > -3000 was chosen acccording to Fabio's DENV2 primer.

# In[64]:


dg_3000 = []
for i in rc_pr_list:
    result = primer3.calcHeterodimer(i,TSO)
    pin = primer3.calcHairpin(i)
    if result.dg > -3000:
        #print(i, result.tm,result.dg, primer3.calcTm(i[23:]))
        #Check for the formation of hairpins.
        #print(i, primer3.calcHairpin(i))
        tttt = [i, primer3.calcTm(i[23:]), result.tm, result.dg, pin.tm, pin.dg]
        dg_3000.append(tttt)
dg_3000 = pd.DataFrame(np.array(dg_3000),columns=["Primer","Annealing Tm","HeteroDimer Tm","HeteroDimer dG","Hairpin Tm","Hairpin dG"])
dg_3000.iloc[:,1:] = dg_3000.iloc[:,1:].astype(float).round(2)
print(dg_3000)

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('dv4_primers.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
dg_3000.to_excel(writer, sheet_name='Sheet1')

# Close the Pandas Excel writer and output the Excel file.
writer.save()


# #### DENV primers used in eLife paper for comparison.

# In[11]:


DENV2 = 'AAGCAGTGGTATCAACGCAGAGTACGAACCTGTTGATTCAACAGC'

print(DENV2, 'Heterodimer', primer3.calcHeterodimer(DENV2,TSO))
print(DENV2, 'Homodimer', primer3.calcHomodimer(DENV2))
print(DENV2, 'Hairpin', primer3.calcHairpin(DENV2))
print(DENV2, 'Tm', primer3.calcTm(DENV2[23:]))

