
# Applied data science in biology

Artificial intelligence কি ও কেন তা নিয়ে খেজুরে গল্পটা আজ আর করবনা । অর্থাৎ গুগল গাড়ি চালাচ্ছে, রোবট সিজারিয়ান করছে, মোবাইল ফোন কড়াই মেজে দিচ্ছে ইত্যাদি প্রভৃতি বিবিধ । এসবের পিছনে যে অঙ্কগুলো সেগুলোও খুব কঠিন । সেসব কিছুই না করে আজকের প্রকল্পঃ এই ঘরে বসে আগামী এক ঘন্টার মধ্যে এই পাতি কম্পিউটারটাকে কিছু '


## Plan

1. Mining well organised data
2. Genome assembly and alignment by hand
2. Genomic data query and manipulation with basic Linux commands
3. Genome assembly software toolchain - FastQC, TrimGalore, SPAdes, QUAST
4. Finding the origin (Ori) site of replication in a bacterial genome with python
5. Finding the regulatory DNA in an operon with python
3. Alignment and homology lookup - Basic local alignment and search tool (BLAST)


```python
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
```

## Warm up exercises

### Calculate

\begin{align}
{5 \choose 2}
\end{align}

### কোন অঙ্ক না কষে, শুধু ভেবে প্রমাণ করো

\begin{align}
{n \choose r} = {n \choose n-r}
\end{align}

### What does this code output?


```python
for x in np.arange(-pi/2,pi/2,pi/4): #x এর মান -pi/2 to +pi/2, pi/4 ধাপ ছেড়ে ছেড়ে
    y = x/x
    s = np.sin(x)
    c = np.cos(x)
    #print ("x: %s, y: %s, s: %s, c: %s" % (x,y,s,c))
    if (y == s**2 + c**2): # '**'' মানে 'to the power of'
        print("Success")
    else:
        print("Error")
```

### এই লেখচিত্রটা (graph) কেমন দেখতে হবে?


```python
x = np.linspace(-1,1,1000)
y = x**2
plt.plot(x,y)
```

# Organised data

প্রথমে খুব সহজ সরল একটা table । Table এর সুবিধা হল সবকিছু বেশ গোছানো, এতগুলো row আর এতগুলো column, এখানে শুরু আর ওখানে শেষ । কোন noise নেই । এই বিশেষ table টা কিছু লোকের শারীরিক কিছু mineral এর মাপ । 

[R markdown workbook](file:///home/tintin/core/code/R/multivariate/multivar_pca.html)

এই table থেকে পাওয়া গেল, বিভিন্ন metabolite এর সাথে ডায়াবেটিস এর একটা সম্পর্ক (association) । কিন্তু মনে রাখতে হবে, যদি এসব electrolyte না মেপে লোকটার বাড়ির নম্বর, মোবাইল নম্বরের শেষ অঙ্ক, কটা দাঁত নেই, মাধ্যমিকের ইতিহাসের নম্বর এইসব দিয়েও একটা table বানাতাম, তাহলেও একটা সম্পর্ক পাওয়া যেত, এবং লোকটার ডায়াবেটিস আছে কি নেই তা এসব থেকেও মোটামুটি বলা যেত । Common sense বাদ দিয়ে জটিল কুটিল অঙ্ক করে 'একটা কিছু' বার করা, এ বিজ্ঞান নয়, বিড়ম্বনা মাত্র । 

# Noisy data

## Genome assembly by hand

আজকাল সব ভারী ভারী যন্ত্র হয়েছে, তারা একখানা বিশাল বড় DNA দিলে নিমেষের মধ্যে এক বিলিয়ন (একের পরে ছটা শূন্য) টুকরো করে ফেলে । তারপর সেগুলোকে 'পড়ে' বলে দেয় 'CACGATCGATATAT...' ইত্যাদি প্রভৃতি । কিভাবে করে জানিনা । উজালা ওয়াশিং লিকুইডের মত । ভেজাও ধোও আর হয়ে গেল । 'কিভাবে' প্রশ্ন করা যাবেনা । সব কপিরাইট । এবার এই টুকরোগুলো থেকে জুড়ে 'অরিজিনাল' টা কে বানাতে হবে । সমস্যাটা জটিল । এই সমস্যটার কত ভাল সমাধান হতে পারে, তার ওপর নির্ভর করছে handheld DNA sequencer তৈরি হতে পারবে কিনা । 

Next generation sequencers (NGS) use proprietary technologies to slice up a DNA strand and then sequence. However, the reads obtained from NGS machines are usually very small and piecemeal. The following are the fragments (reads) returned by an NGS machine from a DNA strand. 

```
AATCC
TCCAA
CAATC
CCAAC
ATCCA
```

**সমস্যা** Reconstruct the original DNA from these fragments.

**NOTE** : _All NGS machines read every nucleotide several times (coverage). Thus, some overlap is expected_

### The first approach

Supposing you have a strand 

```
TAATGCCATGGGATGTT
```

If the sequencer produces 3-mers in random order

```
ATG TGC TAA GCC CCA ATG TGG AAT GGG GGA CAT GAT TGT GTT ATG
```

The mind automatically makes an overlap graph of all the 3-mers.

![](overlap_graph_alternate.png)

The challenge is to find a path that visits every node *exactly once*. (**Hamiltonian path**). Although with 3-mers it looks easy, with increasing read length the problem scales beyond comprehension, and is **one of the NP-hard problems**, i.e. there is no general solution in finite time.

### An alternative approach

If we place the 3-mers in the  what you have done is construct a **De Bruijn Graph** of the 3-mers, where each 3-mer is one *edge* of the graph, rather than a *node*. The nodes are occupied by two 'halves' of each 3-mer, and the 3-mers 'connect' related nodes.

![](composition_graph.png)
![](composition_gluing.png)

This problem can be solved by the machine in finite time.

## Querying the genome at the linux command line

[Click here](./plain_text_life.pdf)

## Genome assembly toolchain

এবার পরপর কতগুলো রেডিমেড প্রোগ্রাম চালাব, যারা _Salmonella enterica_ র জিনোমের টুকরো থেকে পুরোটা তৈরি করে দেবে। শুধু তাই নয়, একটা reference genome এর সাথে তুলনা করে কোথায় কোন জিন কোথায় সেটাও বলে দেবে । আগের অঙ্কটা হাতে করে করার পরে, আশাকরি এই গোটা ব্যাপারটা এতটা রহস্যময় লাগবেনা ।

একটা কথা খুব জরুরি । জীবনবিজ্ঞান কিছুমাত্র না জেনে শুধুমাত্র python দিয়ে bioinformatics সম্ভব নয় । নিদেনপক্ষে এক দুবার দুলালচন্দ্র সাঁতরার বই খুলতেই হবে । তার কারণ, কম্পিউটার ভুল করবে প্রচুর । হাতে গরম উদাহরণ: এই ব্যাক্টিরিয়ার জিনোমেই একজায়গায় একটা জিন শুরু হয়েছে TAT দিয়ে (অর্থাৎ কম্পিউটার গুণে বার করেছে যে ওই জিন TAT দিয়ে শুরু) । অতএব তার mRNA শুরু হবে TAT এর বিপ্রতীপ (compliment) দিয়ে, অর্থাৎ AUA. যা বাস্তবে হয়না । সাধারণতঃ ব্যাক্টিরিয়ার সমস্ত mRNA শুরু হয় AUG (মেথিওনিন) দিয়ে । এজাতীয় ভুল নিজের হাতেই সংশোধন করে নিতে হবে । 

## Finding the origin of replication

https://stepik.org/lesson/2/step/1?unit=8231


```python
# Replication
from collections import Counter
flatten = lambda l: [item for sublist in l for item in sublist]

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def SkewArray(Genome):
    skew = []
    skew.append(0)
    l = len(Genome)
    for i in range(l):
        if Genome[i] == 'G':
            skew.append(skew[i] + 1)
        elif Genome[i] == 'C':
            skew.append(skew[i] - 1)
        else:
            skew.append(skew[i])
    return skew




def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def Reverse(Pattern):
    n = len(Pattern)
    rev = ''
    for i in range(n):
        rev += Pattern[n-1-i]
    return rev
    
def Complement(Pattern):
    comp = {'A':'T','G':'C','T':'A','C':'G'}
    comp_str = ''
    for i in range(len(Pattern)):
        comp_str += comp[Pattern[i]]
    return comp_str
    

def ReverseComplement(Pattern):
    return Complement(Reverse(Pattern))
    # your code here

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    k = len(Pattern)
    l = len(Genome)
    for i in range(l):
        if Genome[i:i+k] == Pattern:
            positions.append(i)
    return positions

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key] == m:
        	words.append(key)
    return words

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for j in range(n-k+1):
            if Text[j:j+k] == Pattern:
                freq[Pattern] += 1
    return freq

# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
        # add each key to words whose corresponding frequency value is equal to m
    return words

def MinimumSkew(Genome):
	positions = []
	skew_array = SkewArray(Genome)
	lowest_vals = []
	m = min(skew_array)
	positions = [k for k,v in enumerate(skew_array) if v == m]
	return positions

def MaximumSkew(Genome):
	positions = []
	skew_array = SkewArray(Genome)
	highest_vals = []
	m = max(skew_array)
	positions = [k for k,v in enumerate(skew_array) if v == m]
	return positions

#print(MinimumSkew('CATTCCAGTACTTCGATGATGGCGTGAAGA'))
# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    l = len(p)
    hd = 0
    for i in range(l):
        if p[i] != q[i]:
            hd +=1
    return hd
    # your code here

def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    l = len(Text)
    p = len(Pattern)
    for i in range(l):
        if len(Text[i:i+p]) == p:
            if HammingDistance(Text[i:i+p],Pattern) <= d:
                positions.append(i)
    # your code here
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    positions = [] # initializing list of positions
    l = len(Text)
    p = len(Pattern)
    for i in range(l):
        if len(Text[i:i+p]) == p:
            if HammingDistance(Text[i:i+p],Pattern) <= d:
                positions.append(i)
    # your code here
    return len(positions)

dict = {'A':0,'C':1,'G':2,'T':3}
def PatternToNumber(Pattern):
	k = len(Pattern)
	base_4_number = ''
	for i in range(k):
		base_4_number = base_4_number + str(dict[Pattern[i]])
	#print(base_4_number)
	return int(base_4_number,4)

def base_convert(i, b):
    result = []
    while i > 0:
            result.insert(0, i % b)
            i = i // b
    return result

def NumberToPattern(Number, k):
	base_4_num = []
	for j in base_convert(Number,4):
		for key,v in dict.items():
			if v == j:
				base_4_num.append(key)
	diff = k - len(base_4_num)
	if diff > 0:
		base_4_num.insert(0,'A'*diff)
	return ''.join(base_4_num)


def ComputingFrequencies(Text, k):
	FrequencyArray = []
	for i in range(4**k):
		FrequencyArray.append(0)
	#print(len(FrequencyArray))
	for i in range(len(Text) - k + 1):
		Pattern = Text[i:i+k]
		j = PatternToNumber(Pattern)
		#print(Pattern,j)
		FrequencyArray[j] = FrequencyArray[j] + 1
	return FrequencyArray


def FasterFrequentWords(Text, k):
    FrequentPatterns = []
    FrequencyArray = ComputingFrequencies(Text, k)
    maxCount = max(FrequencyArray)
    for i in range(4**k):
    	if FrequencyArray[i] == maxCount:
    		Pattern = NumberToPattern(i, k)
    		FrequentPatterns.append(Pattern)
    return FrequentPatterns

def FindingFrequentWordsBySorting(Text , k):
	FrequentPatterns = []
	Index = []
	Count = []
	for i in range(len(Text) -k + 1):
		Pattern = Text(i, k)
		Index.append(PatternToNumber(Pattern))
		Count[i] = 1
		SortedIndex = sorted =(Index)
	for i in range(1, len(Text) - k + 1):
		if SortedIndex[i] == SortedIndex[i - 1]:
			Count[i] = Count[i - 1] + 1
	maxCount = max(Count)
	for i in range(len(Text) - k +1):
		if Count[i] == maxCount:
			Pattern = NumberToPattern(SortedIndex[i], k)
			FrequentPatterns.append(Pattern)
	return FrequentPatterns

def BetterClumpFinding(Genome, k, t, L):
    FrequentPatterns = []
    Clump = []
    for i in range (4**k):
        Clump.append(0)
    Text = Genome[0:L]
    FrequencyArray = ComputingFrequencies(Text, k)
    for i in range(4**k):
        if FrequencyArray[i] >= t:
            Clump[i] = 1
    for i in range(1, len(Genome) - L +1):
        FirstPattern = Genome[i - 1: i - 1 + k]
        index = PatternToNumber(FirstPattern)
        FrequencyArray[index] = FrequencyArray[index] - 1
        LastPattern = Genome[i + L - k: i+L]
        index = PatternToNumber(LastPattern)
        FrequencyArray[index] = FrequencyArray[index] + 1
        if FrequencyArray[index] >= t:
            Clump[index] = 1
    for i in range(4**k):
        if Clump[i] == 1:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

#print(BetterClumpFinding('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA',5,4,50))

def ImmediateNeighbors(Pattern):
	nucleotides = "ACGT"
	Neighborhood = [Pattern]
	for i in range (1, len(Pattern)):
		symbol = Pattern[i]
		for n in nucleotides:
			if n != symbol:
				Neighbor = Pattern[:i-1] + n + Pattern[i:]
				Neighborhood.append(Neighbor)
	return Neighborhood

def FirstSymbol(Pattern):
	return Pattern[0]

def Suffix(Pattern):
	return Pattern[1:]

def Neighbors(Pattern, d):
	nucleotides = "ACGT"
	if d == 0:
		return Pattern
	if len(Pattern) == 1: 
		return ['A', 'C', 'G', 'T']
	Neighborhood = []
	SuffixNeighbors = Neighbors(Suffix(Pattern), d)
	for Text in SuffixNeighbors:
		if HammingDistance(Suffix(Pattern), Text) < d:
			for n in nucleotides:
				Neighborhood.append(n + Text)
		else:
			Neighborhood.append(FirstSymbol(Pattern) + Text)
	return Neighborhood

def FrequentWordsWithMismatches(Text, k, d):
	Neighborhoods  = []
	for i in range(len(Text) - k +1):
		Neighborhoods.append(Neighbors(Text[i:i+ k], d))
	NeighborhoodArray = flatten(Neighborhoods)
	c = Counter(NeighborhoodArray)
	highest_count = c.most_common(1)[0][1]
	l=list(filter(lambda t: t[1]>=highest_count, c.most_common()))
	return [x[0] for x in l ]

def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
	FrequentPatterns = []
	Neighborhoods  = []
	Index = []
	Count = []
	for i in range(len(Text) - k +1):
		Neighborhoods.append(Neighbors(Text[i:i+ k], d))
	rc = ReverseComplement(Text)
	for i in range(len(rc) - k + 1):
		Neighborhoods.append(Neighbors(rc[i:i+k],d))
	NeighborhoodArray = flatten(Neighborhoods)
	c = Counter(NeighborhoodArray)
	highest_count = c.most_common(1)[0][1]
	l=list(filter(lambda t: t[1]>=highest_count, c.most_common()))
	return [x[0] for x in l ]
```


```python
FrequentWordsWithMismatchesAndReverseComplements('ATGCCGGAATCGCATGTAGCTATTGACATGCAGATGCATGGATGTACCATGTGCACTGAT',5,1)
```




    ['ATGCA', 'TGCAT']




```python
PatternCount('TATACGACTGTACTCGATCACAGTTGGATTAGGCCACCATTTATATCA','ATC')
```




    2



## Finding the regulatory gene in an Operon

https://stepik.org/lesson/23063/step/1?unit=6796


```python
#Regulatory motifs

from random import randint, uniform
from math import log

def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    # your code here
    return count


motifs = [['A','T','G','C','C','G','A'],['A','T','T','C', 'G','T','C']]

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    # insert your code here
    count = Count(Motifs)
    for k in count:
        v = count[k]
        v2 = [x / t for x in v]
        profile[k] = v2
    return profile

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs):
    consensus = Consensus(Motifs)
    l = len(Motifs)
    k = len(Motifs[0])
    scores = []
    for i in range(k):
        scores.append(0)
        for j in range(l):
            if Motifs[j][i] != consensus[i]:
                scores[i] += 1
    return sum(scores)

prof = [[0.4,0.3,0.2,0.1],[0.3,0.2,0.5,0.4],[0.3,0,0.3,0],[0.4,0.7,0.2,0.8]]

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(text, k, profile):
    # Inputs a text string, takes k-mers from it and tries it on profile, finds out which has the most probablity of arising
    l = len(text)
    pr = {}
    for i in range(l-k+1):
        kmer = text[i:i+k]
        p = Pr(kmer,profile)
        pr[kmer] = p
    sorted_by_value = sorted(pr.items(), key=lambda kv: kv[1],reverse=True)
    maxval = sorted_by_value[0][0]
    return maxval

def GreedyMotifSearch(Dna, k, t):
    # takes Dna, which is a set of strings, k is length of k-mer, and t is number of k Dna strings
    # 
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def CountWithPseudocounts(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    # your code here
    return count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    # insert your code here
    count = CountWithPseudocounts(Motifs)
    sum_cols = []
    count_vals = []
    for key in count:
        count_vals.append(count[key])
    for i in range(k):
        sum_col = 0
        for j in range(4):
            sum_col += count_vals[j][i]
        sum_cols.append(sum_col)
    str = "ACGT"
    for key in str:
        profile[key] = []
    for i in range(k):
        for j in range(4):
            key = str[j]
            profile[key].append(count_vals[j][i]/sum_cols[i])
    return profile

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    # takes Dna, which is a set of strings, k is length of k-mer, and t is number of k Dna strings
    # 
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def Motifs(Profile, Dna):
    new_motif = []
    t = len(Dna)
    k = len(Profile['A'])
    for i in range(t):
        m = ProfileMostProbableKmer(Dna[i],k,Profile)
        new_motif.append(m)
    return new_motif

def RandomMotifs(Dna, k, t):
    arr = []
    for i in range(t):
        j = randint(1, len(Dna[i]) - k)
        v = Dna[i][j-1:j-1+k]
        arr.append(v)
    return arr

def RandomizedMotifSearch(Dna, k, t):
	M = RandomMotifs(Dna, k, t)
	BestMotifs = M
	while True:
		Profile = ProfileWithPseudocounts(M)
		M = Motifs(Profile, Dna)
		if Score(M) < Score(BestMotifs):
			BestMotifs = M
		else:
			return BestMotifs

def Normalize(Probabilities):
    normalized = {}
    sum_prob = 0
    for k,v in Probabilities.items():
        sum_prob += v
    for k,v in Probabilities.items():
        normalized[k] = v/sum_prob
    return normalized


def WeightedDie(Probabilities):
    kmer = '' # output variable
    # your code here
    prob_ranges = {}
    running_sum = 0
    for k,v in Probabilities.items():
        running_sum_old = running_sum
        running_sum += v 
        prob_ranges[k] = (running_sum_old, running_sum)
    r = uniform(0,running_sum)
    for k,v in prob_ranges.items():
        if r > v[0] and r <= v[1]:
            kmer = k
    return kmer

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    Motifs = RandomizedMotifSearch(Dna, k, t)
    BestMotifs = Motifs
    # your code here
    for j in range(N):
        i = randint(0,t-1)
        motifs_except_ith = Motifs[:i] + Motifs[i:]
        profile = ProfileWithPseudocounts(motifs_except_ith)
        Motifs[i] = ProfileGeneratedString(Dna[i], profile, k)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def MotifEnumeration(Dna, k, d):
	Patterns = []
	for i in range(len(Dna[0])-k + 1):
		pattern = Dna[0][i:i+k]
		neighbors = Neighbors(pattern,d)
		for p in neighbors:
			present_in_str = []
			for dna in Dna:
				j = ApproximatePatternMatching(dna,p,d)
				if len(j) == 0:
					present_in_str.append(False)
				else:
					present_in_str.append(True)
			if all(s==True for s in present_in_str):
				Patterns.append(p)
	return set(Patterns)


def MinDist(Pattern, Text):
	k = len(Pattern)
	dist = []
	for i in range(len(Text) -k + 1):
		h = HammingDistance(Pattern, Text[i:i+k])
		dist.append(h)
	return min(dist)

def DistanceBetweenPatternsAndStrings(Pattern, Dna):
	k = len(Pattern)
	dist = 0
	for dna in Dna:
		h = len(Pattern)
		for i in range(len(dna) -k + 1):
			kmer = dna[i:i+k]
			if h > HammingDistance(Pattern, kmer):
				h = HammingDistance(Pattern, kmer)
		dist += h
	return dist

def MedianString(Dna, k):
	dist = k * len(Dna) # Max possible distance
	for i in range(4**k):
		Pattern = NumberToPattern(i,k)
		currentDist = DistanceBetweenPatternsAndStrings(Pattern, Dna)
		if dist > currentDist:
			dist = currentDist
			Median = Pattern
	return Median


def newlines_to_list(str):
	return str.split('\n')

def list_to_newlines(l):
	return '\n'.join([i for i in l])

def Entropy(prob):
	sum = 0
	for p in prob:
		if (p == 0):
			x = 0
		else:
			x = p * log(p,2)
		sum += x
	return -sum

str='''GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG
CTCGCTGGTCATAATGCATCCTTGCGATTGGGGTGTAACCTCCACAGACAAATCCGACCACCGCGCCGGAGCCCTCCGATAGGGCCTGAAGTCCCGGCCACTGGGGTCATTCGTCCGCCAGGGGCGAAAGGTCTGCCACATTGGGGACTTCCGGCCCTAACTGAGCCGGCGCCATCGAAGCCAGGGTTAGCGCATGGCTGTCCGGGCGGGGCGCGGAGTTCGCCGGAGTCCGAGACCCCGGATCGTGTCG
GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC
GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG
ACGGGACTGGCGACCGAGGTCAGCGATCCCGAGCAGGTTGCCCGCTACCAGCGGCTGCTACACCCGTGGGTGAACATGGCGATGGACACCGTGGTCGCGATCGAACCCGAGATCGTCACCGGCATCCGCATCGTTGCTGACTCGCGTACGCCGTAGCCGATTGGCCGCGGGCGGCCCGCACGCATCCGCACTATCTGATAAATTCTTCAACTCGTCAACCGATGTAACGCTGAAGCTCTCAGGAGACGCG
GCATCCGCGTCCGGGCAACCAGGGACGTTTCGGCCACGTCGCAAGCGACAGCCCAGGGTGCCGGATGCCCAGCCGCTCGCGGTGGCGGGCCGTCCAGAGCCGTCGGCAGCCAGTGACCAGCTGGGGTATCCGCATGGGGTCGCCCAGCGGGTCCCGAGGGGACTTTTGGCCACCGGCGCTGGTGGCCTACTGCCCTCCCGCCGTTGCGCCGGGTGCGTGCACGATTGAAGTCCCCAAGGAAGGGACGCTC
ACGATAGATCAACCGGACCACCGAGGGTTGATTATTGAGGTGCGCTCATCCGATGGTTCGCCGCCGTATGTGGTGCGCTGGCTCGAGACCGACCATGTGGCGACGGTGATTCCGGGTCCGGATGCGGTCGTGGTCACTGCGGAGGAGCAGAATGCGGCCGACGAGCGGGCGCAGCATCGGTTCGGCGCGGTTCAGTCGGCGATCCTCCATGCCAGGGGAACGTAGGCGATTCGCTCAAGCGACGAAGTCG
CGATCACCTCAATCCGATGACCTATGAAGCGCGCGGGCCCGGCCGCCATCGGCCCGTCGATCCGAGTGCGCACGGCCGAAGTGAGCCACCACCGTAGCGCCGCCGAGTTCGCTTCCGCGGACGCAAGCCCGGGATTTGCGGAGTAGCGTACAAGGCGGGATGCGCAATCCGGTCGCTGAGGGACTGTCGTCAGTCGGCCCAGGGACCTAATTCCATATTTGGAGCCATGAGATGGCGACACACTCGGCAG
GCAAACCGTGAGGCTAGGGAAGCGAGGAGCACATGGCCGCCGACCCGCAATGTACACGCTGCAAGCAAACCATCGAACCCGGATGGCTATACATCACCGCCCATCGCCGCGGTCAAGCCGGGATCGTCGATGACGGCGCAGTACTGATTCACGTGCCCGGTGAATGCCCGCACCCCGGGGAGCACGTTCCGCGAAGCTAGCCGGTTGGCATGCTAGGCAACGCAGCGTCAAGTGAGAGGGCAAAACCAAA
CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA
TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC
ATGGCGATCCGCTCGATACTCCGCTGGCAGCCCGAGCGGGGGGTGTCCCCTACCCTGGCTGCCCCGGCTGACGGTATTCACCGGTATCCCATGATGCACTCAATCCCGATGCCTGCACCTTCTGTAAGCTCGCGACGGCTAACGGGGTTATCGGATGGCCCGGGCTCGAAGGAGGTTGGCACGTGTGGATGTTCCTGCCCGGTCGCCACGCGGCGGCGGCGATACCCCAGAGTCGTCAGGTGTTTAGCCG
GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC
TGTCCGCTGTCGAGTTTCGCGACATGGCCGAGACCTTTGAAGACGAGGAGCACCGGCGCTTTGGCGAGGCCGGTTTTCAATCGGTGGTCGACAAGGTCGCCGATATCGAAAAAAGCCTTGGCATCTACGACCTGAGCCAGTTCACCCCCAGCTAAAGACACTAATGCCCTTGGGTTAGGGACCATCGCCTCCTGACGCGATCGCGACAGCTGGCTAACGTCGGTAGTACACCCATGCAGAGGGGACGCCA
TCCGGTTATCGTCGCCCGAATCCCCCAAGATCCGGCAGTGCCGGCCTGAGGGCCTGTGCGATCTGCTCGGGTGGTGCCCACCCGCGCGGAAAGCCCCGTCCGAACCGTGATTGGGCAACGTCGGGCCGGGCCAGCAGCGCTGGACCGTAGGTCCCTGCAGTGGATGACTTACGGCCCTGATCCACACCGGCGACCGTTAGGCAGGGTTGAGCCAACCGTCGGTTGAGCGTCTGGCTGCGAGGTGAGGTGA
TTCGGCGACGGTTGACGTGATGGCCGACAGAATCGTTGCTAGAGCCGGGGGCAACTCCGACGCCACCGCCGAGATCGCGGCGGCCTTGGCCGCCCGGCAAGCGGACTGGGACACCGGGCACCGGATCGACACGGCGGGGCCGCGTGAGCGCTCCGTGGGACAGGCCTACCACATCTGGCGCAGCGCGATCTGAGAACGCCGAAAGGAAAACCGATGCCAACCATCACTGTCAGCAGCACATCGTCGCTGT
CGTCCATGCCGCGGCGTGCACTCGCGTGGCCTTGTCGACCACGTTGTCGAGGCCGACGATGACAGCGTCGTGGTAGCGCCGGTCGATGGTGGGCGGCGGTCCGGTCACATTGACGACGATGCCACAGCGCACGCTTGGATTCGGGCGTCGAATCCAATCATGGTCGGTTGGCCGTCCGATTGGGGACTAAAGCCTCATGACCGGTGACTGTCCCGGTTCATAGCGTTGCTCGAGGCGAAGGGAGAATGAG
CACCACGTGGACCACGGTCAGCGGAATGTTCCTCATCGCCGCATCGGTGGCACCCCAACAGGCGGCGGCATCCGATTCGAGCGAACCATCTACCCCGACGACAACTCCGTGCTGCTTGCGGGGTTTAGACATCTCATTCTCCCTTCGCCTCGAGCAACGCTATGAACCGGGACAGTCACCGGTCATGAGGCTTTAGTCCCCAATCGGACGGCCAACCGACCATGATTGGATTCGACGCCCGAATCCAAGC
TTCCCGCGGATCAGATCTTGACCACCGGGAGTGTCGATGAACTTCTCGCGCTCTTGAAATGACGGGCTATCGTAAGTTTATGGCCTGGGGGAGCGTGAATCCCGCTGGCGGTCGGGTGAACCGCCCCGGTTTTCTTGCACCCCGCGTCGACGTGCCAGTGACGAACTTGACGAATAAGGCCTTTGGTCCTTTCCGGTAGGGGTCTTTGGATAGGCGCGATCCTCGGCATCGGGCCGGTAGCTTGCCGTTT
GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA
AAACTCGGGGAAGAGGGACCGCGGGTGGCGCTGAACGGGAAGGGTGGTGGCCATTTGATGCCTCCTAATCGATGGAAACGGATGCCTTTGATCCGACCAGCCCATCGTGGCCAGGGCTAGGGACAGAAGTCCCCGAAGCGCGGGCCATTTGTCCGCGCCCGTCGGTGATCCACTTGGGGACCATTGACCCTGTTGTCTGCCAACCGCCGTTCAGAAAGATCGGGGTGATATCGAACAGCGGAGGTTGATC
CATCCATCGTGCCGCGCGTCGGCGAGTCCTGTTGATGCGCCACACATTCCGCAGGCATCGTGAACGCTTGACAGCCGTCGCCATTGTCGCGCACAAACCGCACGTCGATTCGTGATCCATTGAGGACCTAAGCCCGTTGGGCTAGTGACAAACGCCTCCTGAGCAAAACCCTCCTCCCCCGTTACCGTCGTGCGGTAGGGACAAGCCACATCGGCCGAGCGGGCGATCAGCCAACGACAGGAGGACCGCG
GTCGATGCGAACGCAAGCATCCAGGAGATGCTCAACGTCATGGAAGAACATCAGGTCCGCCGTGTTCCGGTCATCTCAGAGCACCGCTTGGTCGGAATCGTCACCGAAGCCGACATCGCCCGACACCTGCCCGAGCACGCCATTGTGCAGTTCGTCAAGGCAATCTGCTCGCCCATGGCCCTCGCCAGCTAGCGACCTCGGATCATCCGCCGGCACGAAAGGCTTTATCCGCCAGCATCGGAGGTACCCA
CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG
ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC
CACCTGTGTGGTGGTGGCCACTGTCGAGACCGGCGGCCCGTTGTGGTGGCCCAAGTGCCCTAAGGTGATCAGGTGCCGCAGCCCGGCCAGCACGCCGTCAGAGTTTCACGGGGCTTGGTCGCGGCCGATGGCGTCCTCATCGTGGGGTCGATGACCGAGGTGGACGCGGCGCGACCGGGCACATCGACGGTCCCCGGGGCTTTGTGGGCCAGTGAAGTGACGAAAGACCCCAGTGGACACGGACTTCGGC
GTGCCACGACGCATAGAAGACATGGCCATGCCACACCCTGATAGCATTGCAGCAAGCTACATGTACTGCTCTACCAGGATCCTTATGGGCAACAGTGGGTTTGAGTTATGAAACCCGTGGGCACATACCCTTCCGCGTCGTACTGGTCAGTCTCGACAGCGAAGAGATCACCGGTTGATCCACCAAGCATGCATTGGCGGGCATCTGCATAAACGGTGACGTATCAGCACAAAACAGCGGAGAGAACAAC
GCGAGTACGCTAGTTCAGGTGGGTGCCGTGCCGAAGGCGGTGTCACTCAACGAACTTCGGTTCTCGCAGGGTCGCCACGGCTGGCGATGTGCGGTAACGCTCGATGTGTGAATTGAGACCTGATTCATGAAAATCGTCGAGGAGACCCCATACCGGTTCCGGATCGAACAAGAGGGCGCGATGCGGGTGCCCGGGATCGTGTTCGCGTCCAGGTCGTTGCTGCCTCGTGACGAAGGCGACATGGCCCTTG
CTTGAGGTCCAGTCCGGCGCAGAACACCGGATCGGCGCCGGTGAGGATGACGACGTCGATGTCGTCGTCGGCCTCGGCGTCGGCCAACGCCGCGAAAAACCGATCCCGTAGCGCCGCCGAGAGCGCGTTGCGGGACTGCGGCCGGTTGAGGGTGAGGGTTCGCACCCGTTCGTCGGTGTCGATCAGCAGGATGTCGTCGGTCATTCGATCACCGTAACAGGACCGTTAGACTGTCCTAATGACCAGAAAA
CGCGCGAGATCTCATCGACAACCCACTTCCCATGCCTCACGACGGTCACCATGTCGCGGGCATATTTACGTGAGGCACCGAGGGTGTTTCGCGGGCATTCTTGGTGAGTCAAGTCGAACGGTTGAGCCATGATCGACGATTCCGTTACCGTGCTGTCAGAAGACGAAAGTTGGCACCGGCTGGGCAGCGTTGCACTCGGTCGGCTAGTTACCACCTTTGCTGATGAGCCTGGGATCTTCCAGTCAATTTC
GGGGGCCCGGACGGCCAGGGTGAGAACCGTTCGCACGGTTTCGGCGTCCGGGAAATGGGTGTTCATGGCTGCCGGGCCTTTCCCATCAACGGCTGCGGTCATCTATTAAGGATCGCGCGTCGACAGCCGCGGTGGCAGAGCACGAAGGCTCGCCAGCGGAGGACCTTTGGCCCTGCCTCAACGGCTCCTCGCAGCGGAGAGTGGTACTGACGGTTGCAGATAATCCGGATCACCGGGGAAGGCGCTGACC
ACTCACGTGCCGATCCACGTCTTCTGCCTTGAGAAACCCGGCGTCAAGTGTCGTTAGGTGATTCATGGTCAGCGCCTTCCCCGGTGATCCGGATTATCTGCAACCGTCAGTACCACTCTCCGCTGCGAGGAGCCGTTGAGGCAGGGCCAAAGGTCCTCCGCTGGCGAGCCTTCGTGCTCTGCCACCGCGGCTGTCGACGCGCGATCCTTAATAGATGACCGCAGCCGTTGATGGGAAAGGCCCGGCAGCC
GTTGGCGCATGTACACCTGAGCCGTCGGCTCGCCCACTGGACCCGGCTCTACCCCGAGGTGCGGGTGGATCGGGCCATCGCCGGCGGCAGTGCGTGCCGTCATCTGGCCGCCAACGCAAAGCCGGGTCAGCTGTTCGTCGCGGACTCACACTCCGCGCACGAATTGTGCGGTGCATACCAGCCCGGATGCGCCGTACTTACGGTACGCAGTGCCAACTTGTAGGGAGCGGATCTTGGGAGTGGTGCCCTG
GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG
GTTGTCGATGAGATCTCGCGCGGAGGTCACCAGCAGGTACGCCAAGGCGTATGTGCAGGCTTTGAAGAAGAGCCGGGGCCGGATTTTCGACCAGGTGGTTGACCTGACGGGCTAGTCACGTGATAACGCGCGGCGCCGGCTTGTCGCAGCGGCCAAGCTATCGCCGGGGCTGGGCCGCAGTGTTGCCAAGCGGCGGCGCAAACCGCGTTCGCTGAAGTACTCCTATGACGCGCTGAAGGTGTTGCAGAGG
TAGAGGGCTCCGACGTGCCGGTGCCAGCCGCCGCGTTCGAAACACAGCCCTAACGACACGCTGCCGAATATGACCCGTGTCGGAAATTAGGGCGACAAGAGTAATGCGGCTCAACATAGCCTTGCTTTACTTAGGCAAACCTGCCTTCAACCAGGAGGTTATTATCATCCTGTGGTAACTAGGAAAGCCTTTCCTGAGTAAGTATTGCCTTCGTTGCATACCGCCCTTTACCTGCGTTAATCTGCATTTT'''
print(GibbsSampler(newlines_to_list(str),50,len(newlines_to_list(str)),100))
```

    ['GTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACC', 'GCTGTCCGGGCGGGGCGCGGAGTTCGCCGGAGTCCGAGACCCCGGATCGT', 'CCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAG', 'GTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATT', 'CTGACTCGCGTACGCCGTAGCCGATTGGCCGCGGGCGGCCCGCACGCATC', 'GCCGGATGCCCAGCCGCTCGCGGTGGCGGGCCGTCCAGAGCCGTCGGCAG', 'ATGCGGCCGACGAGCGGGCGCAGCATCGGTTCGGCGCGGTTCAGTCGGCG', 'GCGCGCGGGCCCGGCCGCCATCGGCCCGTCGATCCGAGTGCGCACGGCCG', 'AATGCCCGCACCCCGGGGAGCACGTTCCGCGAAGCTAGCCGGTTGGCATG', 'CCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGG', 'GTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTC', 'GTTCCTGCCCGGTCGCCACGCGGCGGCGGCGATACCCCAGAGTCGTCAGG', 'CCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTG', 'GCCTCCTGACGCGATCGCGACAGCTGGCTAACGTCGGTAGTACACCCATG', 'GCCCTGATCCACACCGGCGACCGTTAGGCAGGGTTGAGCCAACCGTCGGT', 'TCGCGGCGGCCTTGGCCGCCCGGCAAGCGGACTGGGACACCGGGCACCGG', 'GGTGGGCGGCGGTCCGGTCACATTGACGACGATGCCACAGCGCACGCTTG', 'GGCACCCCAACAGGCGGCGGCATCCGATTCGAGCGAACCATCTACCCCGA', 'GGCGGTCGGGTGAACCGCCCCGGTTTTCTTGCACCCCGCGTCGACGTGCC', 'CGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTG', 'GTGGCCAGGGCTAGGGACAGAAGTCCCCGAAGCGCGGGCCATTTGTCCGC', 'GCCCGTTGGGCTAGTGACAAACGCCTCCTGAGCAAAACCCTCCTCCCCCG', 'CTTGGTCGGAATCGTCACCGAAGCCGACATCGCCCGACACCTGCCCGAGC', 'CCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACG', 'GGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGG', 'GCCCGGCCAGCACGCCGTCAGAGTTTCACGGGGCTTGGTCGCGGCCGATG', 'GATCACCGGTTGATCCACCAAGCATGCATTGGCGGGCATCTGCATAAACG', 'CGGTTCTCGCAGGGTCGCCACGGCTGGCGATGTGCGGTAACGCTCGATGT', 'TGTCGTCGGTCATTCGATCACCGTAACAGGACCGTTAGACTGTCCTAATG', 'CCACTTCCCATGCCTCACGACGGTCACCATGTCGCGGGCATATTTACGTG', 'GGATCGCGCGTCGACAGCCGCGGTGGCAGAGCACGAAGGCTCGCCAGCGG', 'GGTCCTCCGCTGGCGAGCCTTCGTGCTCTGCCACCGCGGCTGTCGACGCG', 'CCTGAGCCGTCGGCTCGCCCACTGGACCCGGCTCTACCCCGAGGTGCGGG', 'CGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAG', 'CGCGCGGCGCCGGCTTGTCGCAGCGGCCAAGCTATCGCCGGGGCTGGGCC', 'GGGCTCCGACGTGCCGGTGCCAGCCGCCGCGTTCGAAACACAGCCCTAAC']

