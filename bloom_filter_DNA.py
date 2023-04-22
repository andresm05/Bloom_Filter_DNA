#!/usr/bin/env python
# coding: utf-8

# # DNA matching using BloomFilter

# In[8]:


pip install biopython


# In[16]:


pip install bitmap


# In[14]:


get_ipython().system('pip install primesieve==2.3.0')


# In[44]:


import pandas as pd
from bitmap import BitMap
import numpy as np
from primesieve import n_primes
import math
import hashlib
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt


# In[45]:


class DNA_bloom_filter:
    
    def __init__(self, size, file_name, p):
        self._m = size
        self._capacity = math.ceil((-self._m*math.log(p))/math.log(2)**2)
        self._buckets = n_primes(1,self._capacity)[0]
        self._file_name = file_name
        
        #find the value of k functions
        self._k = math.ceil((self._capacity/self._m)*math.log(2))
        
        self._table = BitMap(self._buckets)
        self._functions = [None]*self._k
        
        #Create k functions
        
        for i in range(self._k):
            a = np.random.randint(1, self._buckets-1)
            b = np.random.randint(0,self._buckets-1)
            self._functions[i] = (a, b)
    
    #Get hash functions
    def getFunctions_(self):
        return self._functions, self._buckets
    
    #Get the quantity of bits set in BloomFilter
    def getCount(self):
        return self._table.count()
    
    def map_(self, key):
        # Map using hash SHA-256
        valor_hash = hashlib.sha256(key.encode()).hexdigest()
        
        # From Hex to Decimal
        mapped_key = int(valor_hash, 16)
        
        return mapped_key
    
    #Set the sequence inside BloomFilter
    def setItem_(self, key):
        
        #Map the key
        mapped_key = self.map_(key)
        
        #Place the key inside BloomFilter
        for f in self._functions:
            a, b = f
            position = (a * mapped_key + b) % self._buckets
            self._table.set(position)
    
    #Check if sequence is in BloomFilter
    def findItem_(self, mapped_key):
        for f in self._functions:
            a, b = f
            position = (a * mapped_key + b) % self._buckets
            if(not self._table.test(position)):
                return False
        return True
    
    def checkSample_(self, records, k=31):
        
        #hits[0] counts polluted samples
        #hits[1] counts not polluted samples
        hits = [0,0]
        
        for read in records:
            DNA_seq = str(read.seq)
            
            #Get rid of N characteres in the sequence
            if 'N' in DNA_seq:
                 DNA_seq = ''.join([word for word in DNA_seq if word != "N"])
            
            #Get the K-mers and check if they are inside the BloomFilter
            for i in range(len(DNA_seq)-k+1):
                kmer = DNA_seq[i:i+k]
                mapped_kmer = self.map_(kmer)
                if(self.findItem_(mapped_kmer)):
                    hits[1]+=1
                else:
                    hits[0]+=1
        return hits
    
    #Process the original DNA Sequence
    def process_(self, mapped = False):     
        with open (self._file_name, 'r') as file:
            for line in file:
                self.setItem_(line)
            


# ## Data Base inspection

# In[42]:


df = pd.read_csv('kMer-BloomFilter.txt', header=None, names = ['Cadena'])


# In[26]:


df.head()


# In[43]:


df.tail()


# In[5]:


df.shape[0]


# ## Create BloomFilter using false positive rate of 0.1%

# In[46]:


pure_DNA = DNA_bloom_filter(10000000, 'kMer-BloomFilter.txt', 0.001)


# In[5]:


pure_DNA._k


# In[6]:


pure_DNA._capacity


# ## Complete BloomFilter hash table

# In[47]:


pure_DNA.process_()


# In[12]:


pure_DNA.getCount()


# ## Testing

# ### Physarum Polycephalum

# In[13]:


#Get DNA sequences
records_physarum = list(SeqIO.parse("physarum.fastq", "fastq"))


# In[17]:


print(records_physarum[6].seq)


# In[19]:


#Compare sequence with BloomFilter
hits_physarum = pure_DNA.checkSample_(records_physarum)


# In[20]:


hits_physarum[1]


# In[21]:


colors = sns.color_palette('pastel')
plt.pie(hits_physarum, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Physarum polycephalum', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# ### Salmonella enterica subsp. enterica serovar

# In[22]:


#Get DNA sequences
records_salmonella = list(SeqIO.parse("salmonella.fastq", "fastq"))


# In[23]:


#Compare sequence with BloomFilter
hits_salmonella = pure_DNA.checkSample_(records_salmonella)


# In[35]:


hits_salmonella[0]


# In[24]:


colors = sns.color_palette('pastel')
plt.pie(hits_salmonella, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Salmonella enterica subsp. enterica serovar', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# ### Homo sapiens sapiens (Blood sample)

# In[25]:


records_homo2 = list(SeqIO.parse("homo_sapiens2.fastq", "fastq"))


# In[26]:


hits_homo2 = pure_DNA.checkSample_(records_homo2)


# In[27]:


colors = sns.color_palette('pastel')
plt.pie(hits_homo2, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Blood sample)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# ### Homo Sapiens Sapiens (Dental Calculus)

# In[8]:


records_oral9974 = list(SeqIO.parse("oral9974.fastq", "fastq"))


# In[9]:


str(records_oral9974[1].seq)


# In[10]:


hits_oral9974 = pure_DNA.checkSample_(records_oral9974)


# In[11]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral9974, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Ancient dental calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[43]:


records_oral7742 = list(SeqIO.parse("oral7742.fastq", "fastq"))


# In[44]:


hits_oral7742 = pure_DNA.checkSample_(records_oral7742)


# In[61]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7742, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens \n(Human metagenome isolated from archaeological dental calculus of adult female with periodontal disease.)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[47]:


records_oral7743 = list(SeqIO.parse("oral7743.fastq", "fastq"))


# In[48]:


hits_oral7743 = pure_DNA.checkSample_(records_oral7743)


# In[59]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7743, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens \n(Human metagenome isolated from archaeological dental calculus of adult male with periodontal disease)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[50]:


records_oral7131 = list(SeqIO.parse("oral7131.fastq", "fastq"))


# In[54]:


hits_oral7131 = pure_DNA.checkSample_(records_oral7131)


# In[58]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7131, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Pre-historic California Dental Calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[12]:


records_oral6712 = list(SeqIO.parse("oral6712.fastq", "fastq"))


# In[13]:


hits_oral6712 = pure_DNA.checkSample_(records_oral6712)


# In[19]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral6712, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Late 19th and early 20th century dental calculus metagenomes)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[15]:


records_oral7008 = list(SeqIO.parse("oral7008.fastq", "fastq"))


# In[16]:


hits_oral7008 = pure_DNA.checkSample_(records_oral7008)


# In[20]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7008, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Ancient human oral microbiome)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[39]:


records_oral9976 = list(SeqIO.parse("oral9976.fastq", "fastq"))
print(records_oral9976[0].seq)


# In[48]:


hits_oral9976 = pure_DNA.checkSample_(records_oral9976)
hits_oral9976[0]


# In[49]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral9976, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Ancient Dental Calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[24]:


records_oral7305 = list(SeqIO.parse("oral7305.fastq", "fastq"))


# In[25]:


hits_oral7305 = pure_DNA.checkSample_(records_oral7305)


# In[34]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7305, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (archaeological dental calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[27]:


records_oral7132 = list(SeqIO.parse("oral7132.fastq", "fastq"))


# In[28]:


hits_oral7132 = pure_DNA.checkSample_(records_oral7132)


# In[33]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7132, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (Pre-historic California Dental Calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[35]:


records_oral7303 = list(SeqIO.parse("oral7303.fastq", "fastq"))


# In[36]:


hits_oral7303 = pure_DNA.checkSample_(records_oral7303)


# In[37]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral7303, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (archaeological dental calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[39]:


records_oral6633 = list(SeqIO.parse("oral6633.fastq", "fastq"))


# In[40]:


hits_oral6633 = pure_DNA.checkSample_(records_oral6633)


# In[41]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral6633, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (archaeological dental calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# In[44]:


records_oral2944 = list(SeqIO.parse("oral2944.fastq", "fastq"))


# In[45]:


hits_oral2944 = pure_DNA.checkSample_(records_oral2944)


# In[46]:


colors = sns.color_palette('pastel')
plt.pie(hits_oral2944, labels = ['Contaminado','No Contaminado'], colors = colors, autopct='%0.3f%%')
plt.title('Homo Sapiens (archaeological dental calculus)', bbox={'facecolor':'0.8', 'pad':5})
plt.show()


# # Testing Kmer Function

# In[28]:


def get_kmers_from_file(file, k=31):
    records= list(SeqIO.parse(file, "fastq"))
    kmers = []
    
    for r in records:
        record = str(r.seq)
        for i in range(len(record)-k+1):
            kmer = record[i:i+k]
            kmers.append(kmer)
    
    return kmers


# In[39]:


kmers_oral = get_kmers_from_file("oral1.fastq")


# In[32]:


kmers_oral = pd.DataFrame(kmers_oral)
kmers_oral.head()


# In[34]:


kmers_oral.shape


# In[ ]:


kmers_oral.to_csv('base_datos_prueba.cvs', index = False)

