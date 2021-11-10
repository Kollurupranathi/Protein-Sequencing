"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename)
    text=''
    temp = f.read()
    x=temp.split('\n')
    for i in range(len(x)): 
         text+= x[i] 
    return text
   


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    stop=['TAA', 'TAG', 'TGA']
    rna=[]
    codon=''
    for i in dna[startIndex:]:
        codon+=i
        if len(codon)==3:
            rna.append(codon.replace('T','U'))
            if codon in stop:
                break
            codon=''
    return rna


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    t=open(filename)
    f = json.load(t)
    sample={}
    ac=list(f.keys())
    c=list(f.values())
    for i in range(len(c)):
        for j in range(len(c[i])):
            if c[i][j] not in sample:
                sample[c[i][j].replace('T','U')]= ac[i]
    return sample
    


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein=[]
    for i in codons:
        if i == 'AUG' and len(protein)==0:
            temp='Start'
        else:
            temp = codonD[i]
        protein.append(temp)
        if codonD[i]=='Stop':
            break
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    file=readFile(dnaFilename)
    codondict=makeCodonDictionary(codonFilename)
    synpro=[]
    i=0
    while i<len(file):
        if file[i:i+3]=='ATG':
            u=dnaToRna(file,i)
            synpro.append(generateProtein(u,codondict))
            i+=(3*len(u))
        else:
            i+=1
    return synpro    

def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    protein=[]
    for i in proteinList1:
        for j in proteinList2:
            if i==j and i not in protein:
                     protein.append(i)             
    return protein


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    protein_comb=[]
    for proteins in proteinList:
        for amino_acid in proteins:
            protein_comb.append(amino_acid)
    return protein_comb


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aa_dict={}
    for i in aaList:
        if i not in aa_dict:
            aa_dict[i]=0
        aa_dict[i]+=1
    return aa_dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    k1=combineProteins(proteinList1)
    k2=combineProteins(proteinList2)
    list1=aminoAcidDictionary(k1)
    list2=aminoAcidDictionary(k2)
    key,temp,final=[],[],[]
    for i in list(list1.keys()):
        for j in list(list2.keys()):
            rem=['Start','Stop']
            if i not in rem and j not in rem:
                if i==j and i not in key:
                            key.append(i) 
                if i!=j and i not in key:
                            key.append(i)
                if i!=j and j not in key:
                            key.append(j)   
    for p in key:
        if p in list1 and p in list2:
            if abs((list1[p]/len(k1))-(list2[p]/len(k2)))>cutoff:
                temp.append(p)
                temp.append(list1[p]/len(k1))
                temp.append(list2[p]/len(k2))
                final.append(temp)
        if p in list1 and p not in list2:
            if (list1[p]/len(k1))>cutoff:
                temp.append(p)
                temp.append(list1[p]/len(k1))
                temp.append(0)
                final.append(temp)
        if p not in list1 and p in list2:
            if (list2[p]/len(k2))>cutoff:
                temp.append(p)
                temp.append(0)
                temp.append(list2[p]/len(k2))
                final.append(temp)
        temp=[]
    return final


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print('The following proteins occurred in both DNA Sequences:')
    for o in commonalities:
        for w in o:
            rem=['Start','Stop']
            if w not in rem:
                print(w)
    print('The following amino acids occurred at very different rates in the two DNA sequences:')
    for row in differences:
        tempo=row[0] +':' +str(round(row[1]*100,2)) +"% in Seq1,"+ str(round((row[2]*100),2)) +'% in Seq2'
        print(tempo)
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    list1=aminoAcidDictionary(combineProteins(proteinList1))
    list2=aminoAcidDictionary(combineProteins(proteinList2))
    labels=[]
    for i in list(list1.keys()):
        for j in list(list2.keys()):
                if i==j and i not in labels:
                            labels.append(i) 
                if i!=j and i not in labels:
                            labels.append(i)
                if i!=j and j not in labels:
                            labels.append(j) 
    labels.sort()
    return labels


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    f=combineProteins(proteinList)
    list_dict=aminoAcidDictionary(f)
    chart_data=[]
    for i in labels:
        if i in list_dict:
            chart_data.append(list_dict[i]/len(f))
        else:
            chart_data.append(0)
    return chart_data


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    test.testSetupChartData()

    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    
    test.runWeek2()
    
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
