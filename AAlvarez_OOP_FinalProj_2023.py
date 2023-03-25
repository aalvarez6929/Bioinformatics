import re

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

## RENAME this file YourLastName_OOP_FinalProject_2023.py
# after done, import this into a jupitor notebook - exeecutive running. 
##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties. NOTE: You are not allowed
##            to import any specialized libraries for this project (e.g., no Biopython)
##            The idea is for you to write these methods from scratch.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 
### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and strip to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3. 
#
#     -----------------------notes-----------------------------------------------------------
#     kmers = k length of monomers --- can be DNA, RNA, protein sequences of various k lengths
#         i = wherever the monomer starts on the seq AKA kmer_start
#       i:i+k <------this means from i to i+k (0 to 0+3, 1 to 1+3, 2 to 2+3, ...), we need      							#			to stop some distance b/f the sequence since the monomer is k in 						#                       length. it would go over if this is not included
#
#     self.kmers.extend(kmers) <------ this line extends the self.kmers [list] by appending all 								#					the kmers madee in the current method (make_kmers). This 			
#					way self.kmers stores all the kmers generated so far.
#     ----------------------notes end--------------------------------------------------------------                
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT

class Seq:
    def __init__(self, sequence, gene, species):
        self.sequence = sequence.upper().strip()
        self.gene = gene
        self.species = species
        self.kmers = []

    def __str__(self):
        return f"{self.species}, {self.gene}: {self.sequence}"

    def print_record(self):
        return self.sequence

    def make_kmers(self, k=3):
        kmers = [self.sequence[i:i + k] for i in range(len(self.sequence) - k + 1)]
        self.kmers.extend(kmers)
        return self.kmers

    def fasta(self):
        fasta_string = ">" + self.species + " " + self.gene + "\n" + self.sequence
        return fasta_string

### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#    ------->  re.sub('[^ATGC]','N',sequence) <---------will change any character that is not a
#     							 capital A, T, G, C into an N. (Seq already uppercases and strips.)

#  Methods:
#  (1) Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.
#  (2) Add a method called reverse_complement that returns the reverse complement of
#      self.sequence (watson-crick base pairing) 
#  (3) Add a method called six_frames that returns all 6 frames of self.sequence
#      This include the 3 forward frames, and the 3 reverse complement frames

class DNA(Seq):
    def __init__(self, sequence, gene, species, gene_id):
        super().__init__(sequence, gene, species)
        self.gene_id = gene_id
        self.sequence = re.sub('[^ATGC]', 'N', sequence.upper().strip())

    def analysis(self):
        gc = len(re.findall('G', self.sequence) + re.findall('C', self.sequence))
        return gc

    def print_info(self):
        return self.gene_id + " " + self.sequence

    def reverse_complement(self):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        # Reverse the original sequence
        reversed_sequence = self.sequence[::-1]
        reverse_complement_sequence = ""
        # for loop lets us go one at a time through each base and use the complement dictionary.
        for base in reversed_sequence:
            complement_base = complement[base]
            reverse_complement_sequence += complement_base
        return reverse_complement_sequence

    def six_frames(self):
        frames = []
        # adds to frame list <-Slice the sequence starting from the current index (0,1,2...) 
        # andthen for reverse_comp
        for frame_start in range(3):
            frames.append(self.sequence[frame_start:])
        reverse_comp = self.reverse_complement()
        for frame_start in range(3):
            frames.append(reverse_comp[frame_start:])
        return frames
        
        
        
### RNA Class:  INHERITS DNA class
"""---->this means that all of the objects created from the DNA will be brought into RNA class<------"""

#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
"""------>  import re lets us use .replace("","") function to change charachte<-----"""

#  (2) Add self.codons equals to an empty list

#  Methods: make codons, and translate 
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
"""------>self.coodon = []  then use for loop to iterate by 3 through the whole sequence -----------> rnage(0, len(), 4) <---------""" # adds codons to an empty list, must come first. 
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.prtoein sequence will #      be a string of AA also **kwargs will let us add this to init objects

class RNA(DNA):

    def __init__(self, sequence, gene, species, gene_id, **kwargs):
    	super().__init__(sequence, gene, species, gene_id, **kwargs)
    	self.sequence = self.sequence.replace("T", "U")
    	self.codons = []
        
    def make_codons(self):
        self.codons = []
        for codon_start in range(0, len(self.sequence), 3):
            codon = self.sequence[codon_start:codon_start+3]
            # if codon is equal== 3 nt in length add(append) to self.codon [] list 
            if len(codon) == 3:
                self.codons.append(codon)
        return self.codons

    def translate(self):
        standard_code = {
         "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
         "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
         "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
         "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
         "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
         "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
         "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
         "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
         "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
         "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
         "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
         "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
         "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
        protein_sequence = ""
        for codon in self.make_codons():
            if codon in standard_code:
                protein_sequence += standard_code[codon]
        return protein_sequence  


### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  Use re.sub to change any non LETTER characters in self.sequence into an 'X'.

#  Methods:
#  The next 2 methods use a kyte_doolittle and the aa_mol_weights dictionaries.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a self.sequence
#  (3) Add mol_weight, which returns the total molecular weight of the protein
#      sequence assigned to the protein object. 
class Protein(Seq):
    def __init__(self, sequence, gene, species, gene_id):
        super().__init__(sequence, gene, species)
        sequence_upper = sequence.upper().strip()
        self.gene_id = gene_id
        cleaned_sequence = ''
        for aa in sequence_upper:
            if aa in 'ACDEFGHIKLMNPQRSTVWY':
                cleaned_sequence += aa
            else:
                cleaned_sequence += 'X'
        
        self.sequence = cleaned_sequence
    def total_hydro(self):
        kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}
        total_hydrophobicity = 0
        for aa in self.sequence:
            hydrophobicity = kyte_doolittle[aa]
            total_hydrophobicity += hydrophobicity
        return total_hydrophobicity

    def mol_weight(self):
        aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}
        total_mol_weight = sum([aa_mol_weights[aa] for aa in self.sequence])
        return total_mol_weight
    
"""#############------------Lets test it----------------#########################################  
#create the DNA object
dna = DNA(sequence="ATGAGCACCTTGTTCTGACTTATGCTTACCTTATTAGGACCTGACCTGGTTATTTGAACTCGATGTTTTGGTCTAA", gene="tmp_gene", species="temp_species", gene_id=3456)

# create the RNA object with  DNA object
rna = RNA(sequence=dna.sequence, gene=dna.gene, species=dna.species, gene_id=dna.gene_id)

# create the Protein object using  DNA object
protein = Protein(sequence=dna.sequence, gene=dna.gene, species=dna.species, gene_id=dna.gene_id)

#use these generated objects to call the methods to preform special taskss
reverse_complement = dna.reverse_complement()
fasta = dna.fasta()
print("FASTA format:", fasta)
print("\n kmers :", dna.make_kmers())

print("\n Reverse complement:", reverse_complement)

print("\n RNA sequence:", rna.sequence)

make_codons = rna.make_codons()
print("\n codons: ", make_codons)

translated_seq = rna.translate()
print("\n Translated:", translated_seq)

print("\n six frames : ",dna.six_frames())

total_hydro_value = protein.total_hydro()
print("\n Hydrophobicity value: ", total_hydro_value)
    
print("\n mole weight(g/mole): ", protein.mol_weight()) #<--- ez way to print
"""


