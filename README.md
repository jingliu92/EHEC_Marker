# EHEC_Marker

## Genes interested
`eae` (intimin; include a few alleles if you care about subtyping, e.g., β, γ, ε). \
`stx1` and `stx2` (include major subtypes: stx1a/1c/1d; stx2a–2g to avoid false negatives). \
`espK, espV and espN` (T3SS effectors; often on OI-122 / related islands). \
`aggR`: Encodes the AggR protein, which is involved in the expression of several virulence determinants in both plasmids and chromosomes. The AggR protein is crucial for the adherence of EAEC to host cells and is associated with the pathogenesis of diarrhea in both children and adults. \
`subA`: Subtilase cytotoxin. \
`Z2098` (OI-122 marker; locus tag from the O157 Sakai island). \
`ureD` (urease accessory subunit—present in some clades; include a couple of E. coli hits to avoid over-stringency). \

<img width="585" height="361" alt="image" src="https://github.com/user-attachments/assets/6ebcdb81-8a39-4c00-b477-650429ab189d" />


## Tasks:
1. Identify all eae+ isolates
2. Find their linkage to stx, esp, ureD and z2098. Then create the linkage table.
3. Find aggR that are positive for stx or stx/subA
4. Find eae negative but stx+ isolates (STEC). STEC is not important for food safety. 

## Bioinformatics Analysis

`ehec_markers.fna`
```
>eae
CTTCGTCACAGTTGCAGGCCTGGTTACAACATTATGGAACGGCAGAGGTTAATCTGCAGAGTGGTAATAACTTTGACGGTAGTTCACTGGA
>stx1
CGGTTACATTGTCTGGTGACAGTAGCTATACCACGTTACAGCGTGTTGCAGGGATCAGTCGTACGGGGATGCAGATAAATCGCCATTCGTTGACTACTTCTTATCTGGATTT
>stx2
GCAACGTGTCGCAGCGCTGGAACGTTCCGGAATGCAAATCAGTCGTCACTCACTGGTTTCATCATATCTGGCGTTAATGGAGT
>espK
GGACTATCTTCACCCAGAGAAAAACTTGTCAACAGCATTAAACTAGTTAACGACTTTTCCAAAAAAGATATTATAACTCAGAATACGTTAACAAACGCAGTATGGGACCCCAGAACACCTCGCAAATAT-AAGCAAGATCCACTTATCAAAAGAGCATTAAATGAACACGAAAGAGG
>espV
AAAAAGATTGGGGAGGGCTAGAGCCTTCTGCACGAAAATAGTGATGCGACAGATCTGAGATGGCTTTGTATCCTTCTGCCTGGTTCAGGCCTGGAGCAGTCCATGTTTTACAAAAATACTCGGCAGCGTAACGTGTTGCAAGAGCGGCATCAGACGAGGAACCTGAGGTTCCGCTCAT
>espN
ATCCAAATGTAAAAGCACCAGATGAATATATGGCACAACAGCAGGTACGTATTGGAAGAATATTGCCTGAATACACATCGGATCTTAGCGTTGCTGATCGGTTTAGTCGTGAGCATTATCTAATAGTTGTTAAAGTAAAGGCAAAATATATCACACGAGGAAGTGTTACAGAGAGTGGTTG
>subA
CCCTTCTGGCCTCCCGATATATTAATGATTTTCTCTTCTCCTGGTGCAAGAAAAACATTTGACATAATGCTTTCAGTTGCTCTAATCCATGAGTCCTGTACACCGTCAGGAATAACCCTTCTGCTGGATATCAGTGCATGAGGATAAA
>aggR
TATTACTTTTAACCAGATCCTTATGCAATCAAGAATGAGCAAAGCAGCATTGCTGTTGCTTGATAACTCATATCAGATATCACAGATATCTAATATGATAGGATTTT
>Z2098
CTGAGCGAATGCCTGAAAAGAGCCAGAACGTGCTTATTTCGATGAATATCGATAGCGAGGCTGGGCCATTAATATATTCCGCACGCTATCTCGGAGGCACGTTCCGGCGCGGAGGTATAGCAGTTAGTCCGGGTAATGATCTTAGGCAAGCAACCCACTGGATGTCGCTACCAGAACCGCCGCAGGAGGTGAATCAATGA
>ureD
AGAAACCTTCAGCCACGGAACACTCGCTAACCGTCTGGAGATATGGGTCGACGGTTCTCCGCTGTTGATTGAGCGTCTGCACGTGGCCAATGGAGAACTTACCTGCGTTGCCCGGAAACCATGGGTGGGAACGATGCTGTTCTACCTGGGGAATGAAACGCAACTGGAAGACATACGCGAAAAACTCA
```
## Download assemblies from NCBI
```
# Get ALL complete genomes for Escherichia coli 
./datasets download genome taxon "Escherichia coli" \
  --filename E.coli_complete.zip
```
<img width="1326" height="66" alt="image" src="https://github.com/user-attachments/assets/54e1ffe6-672e-4e96-bb2d-cc054760e37e" />
In total, 11,321 genome assemblies were retrieved and downloaded. The assemblies are organized into folders, each named according to the respective NCBI accession number: 
<img width="839" height="108" alt="image" src="https://github.com/user-attachments/assets/57a48e37-c570-4dd6-a83d-6e3bbc03f457" />

## Search target genes in the downloaded assemblies
### Option 1: Use Blast-----large dataset takes forever, very slow, requires high computing compacity
### Option 2: Use package ABRicate.
#### What is ABRicate? (https://github.com/tseemann/abricate)
ABRicate is a bioinformatics tool specifically designed to screen bacterial genome assemblies (or contigs) for the presence of known genes of interest, such as:
- Antibiotic resistance genes (ARGs)
- Virulence factors
- Plasmid-associated genes

It does this by comparing the genome sequences against curated databases. Some commonly used databases include:

- ResFinder – for antimicrobial resistance genes

- VFDB – for virulence factors

- PlasmidFinder – for plasmid replicons

- CARD – for comprehensive antibiotic resistance \
- In our case, we will customize a databased by using the genes that we are interested.

Key features of ABRicate:
- Works on genome assemblies (FASTA files).
- Reports hits with identity and coverage thresholds, so you know how well a gene matches.
- Can scan multiple genomes in batch mode.
- Focused on functional interpretation—e.g., it tells you which resistance genes are present, rather than showing raw sequence similarity.
## Install `ABRicate` in Our Server
```
# Create a clean environment for abricate
conda create -n abricate -c bioconda -c conda-forge abricate

# Activate it
conda activate abricate

# Check that it works
abricate --version
```
## Make Our Custom Marker DB
Once ABRicate is installed, make a custom database for our interested markers:
```
mkdir -p /home/jing/miniconda3/envs/abricate/db/ehec_markers
cp ehec_markers.fna /home/jing/miniconda3/envs/abricate/db/ehec_markers/sequences
abricate --setupdb
abricate --list   # should now show ehec_markers
```
## ABRicate with a custom DB:
```
mkdir -p abricate_out
find all_sample -type f -name "*.fna" -print0 \
| while IFS= read -r -d '' f; do
  base=$(basename "$f" .fna)
  echo "[ABRicate] $base"
  abricate --db ehec_markers --minid 80 --mincov 90 "$f" > "abricate_out/${base}.tsv"
done

abricate --summary abricate_out/*.tsv > abricate_summary.tsv
echo "Wrote abricate_summary.tsv"
```

## Make presence/absence matrix + classification
```
# save as make_matrix_and_classify.py; run: python make_matrix_and_classify.py
import csv, pandas as pd

genes = ['stx1','stx2','eae','espK','espV','Z2098','ureD']
rows=[]
with open("abricate_summary.tsv") as fh:
    r=csv.DictReader(fh, delimiter="\t")
    for row in r:
        rec={'sample':row['FILE']}
        for g in genes:
            rec[g]=int(row.get(g,'0') not in ('','0','0.00','NA'))
        rows.append(rec)

df=pd.DataFrame(rows).sort_values('sample')

def classify(r):
    stx = r['stx1'] or r['stx2']
    support = r['espK']+r['espV']+r['Z2098']+r['ureD']
    if not stx: return "Non-EHEC"
    if r['eae'] or support>=2: return "EHEC"
    if support==1: return "Possible EHEC"
    return "Non-EHEC"

df['CLASS']=df.apply(classify,axis=1)
df.to_csv("marker_matrix.tsv", sep="\t", index=False)
print("Wrote marker_matrix.tsv")
```
