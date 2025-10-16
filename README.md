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
🧬 Step 1: Prepare your custom marker database
We already have all your marker gene FASTA sequences (e.g., eae, espK, espV, Z2098, stx1, stx2, ureD, etc.), put them together into one FASTA file:
```
cat eae.fasta espK.fasta espV.fasta Z2098.fasta stx1.fasta stx2.fasta ureD.fasta > ehec_markers.fna
```

Then create a BLAST database:
```
makeblastdb -in ehec_markers.fna -dbtype nucl -out ehec_markers
```

✅ This will create files:
```
ehec_markers.nsq
ehec_markers.nin
ehec_markers.nhr
```
🧫 Step 2: Run BLAST for one genome assembly for testing purpose:
```
blastn -query ehec_markers.fna \
       -subject /home/jing/E.coli/EHEC/EHEC_O111_H8/EHEC_O111_H8_GCA_049489415.1_ASM4948941v1_genomic.fna \
       -out EHEC_O111_H8_GCA_049489415.1_markers.tsv \
       -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
```

🔍 Explanation of output columns:

Column	Meaning
qseqid	Query name (marker gene ID)
sseqid	Subject (genome) sequence ID
pident	Percent identity
length	Alignment length
qlen / slen	Query and subject sequence lengths
qstart/qend	Start and end of query alignment
sstart/send	Start and end on subject
evalue	BLAST significance score
bitscore	Strength of alignment

🧮 Step 3: Add filtering (optional)
If you only want matches ≥80% identity and ≥80% coverage, you can filter like this:
```
blastn -query ehec_markers.fna \
       -subject /home/jing/E.coli/EHEC/EHEC_O111_H8/EHEC_O111_H8_GCA_049489415.1_ASM4948941v1_genomic.fna \
       -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
       | awk '{cov=($4/$5)*100; if($3>=80 && cov>=80) print $0}' > EHEC_O111_H8_GCA_049489415.1_filtered_hits.tsv
```
🌀 Step 4: Run on all assemblies

If you have many genome folders, you can loop through all .fna files automatically:
```
mkdir -p blast_out
mkdir -p blast_hits_only

for f in $(find /home/jing/E.coli/ecoli_all -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")
  base=$(basename "$f" .fna)
  out_file="blast_out/${folder}_${base}_hits.tsv"

  echo "Running BLAST on $folder ..."

  blastn -query ehec_markers.fna \
         -subject "$f" \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         | awk '{cov=($4/$5)*100; if($3>=80 && cov>=80) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only/
      echo "✅ Hits found in $folder — copied to blast_hits_only/"
  else
      echo "❌ No hits for $folder"
  fi
done
```
<img width="1047" height="103" alt="image" src="https://github.com/user-attachments/assets/54048c16-85f2-4a3a-b6e2-f865b4e398f3" />

### Check how many assemblies have hits? (1984)
```
cd blast_hits_only
ls -l . | wc -l
# 1984
```
<img width="879" height="49" alt="image" src="https://github.com/user-attachments/assets/d7c68fe3-3219-4edc-a5f5-ac21c1ddc0b7" />

### Add headers to each tsv file
```
# Define both possible headers
hdr11=$'qseqid\tsseqid\tpident\tlength\tqlen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'
hdr10=$'qseqid\tsseqid\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'

# Walk all *_hits.tsv (adjust the root path if needed)
find ./blast_hits_only -type f -name "*_hits.tsv" -print0 | while IFS= read -r -d '' f; do
  # Normalize line endings (avoids weird header checks)
  perl -pi -e 's/\r$//' "$f"

  # Skip if header is already present
  if head -n1 "$f" | grep -q $'^qseqid\t'; then
    echo "Header already present: $f"
    continue
  fi

  # Detect column count from the first non-empty line
  nfields=$(awk 'NF{print NF; exit}' "$f")
  if [ -z "$nfields" ]; then
    # empty file: choose 11-col header by default
    header="$hdr11"
  elif [ "$nfields" -eq 11 ]; then
    header="$hdr11"
  elif [ "$nfields" -eq 10 ]; then
    header="$hdr10"
  else
    echo "⚠️  Unknown column count ($nfields) in $f — not modifying."
    continue
  fi

  # Prepend header safely
  tmp=$(mktemp)
  { printf "%s\n" "$header"; cat "$f"; } > "$tmp" && mv "$tmp" "$f"
  echo "Added header to: $f"
done
```

🧾 Step 5: Summarize presence/absence (Python)

- Reads all .tsv files in a folder

- Filters for hits ≥80% identity and ≥80% coverage

- Outputs a presence/absence matrix

nano summarize_blast_results.py
```
import pandas as pd
import glob, os

# Parameters
identity_cutoff = 80
coverage_cutoff = 0.8

all_data = []
for f in glob.glob("*.tsv"):
    df = pd.read_csv(f, sep="\t", header=None, names=[
        "qseqid","sseqid","pident","length","qlen","qstart","qend",
        "sstart","send","evalue","bitscore"
    ])
    df["coverage"] = df["length"] / df["qlen"]
    df = df[(df["pident"] >= identity_cutoff) & (df["coverage"] >= coverage_cutoff)]
    genes = df["qseqid"].unique().tolist()
    all_data.append({"Sample": os.path.splitext(os.path.basename(f))[0], **{g:1 for g in genes}})

# Combine into a single DataFrame
summary = pd.DataFrame(all_data).fillna(0).set_index("Sample")
summary = summary.astype(int)

# Save to file
summary.to_csv("blast_presence_absence.tsv", sep="\t")
print(summary)
```
# Run 
```
python summarize_blast_results.py --glob "/home/jing/E.coli/blast_hits_only/*.tsv"
```
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
# Create output directory
mkdir -p abricate_out

# Loop through all subfolders containing .fna files
for f in $(find . -type f -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")   # get folder name (e.g., GCA_037094545.1)
  base=$(basename "$f" .fna)             # get file name without extension
  outname="${folder}_${base}.tsv"        # combine both for unique output name
  echo "Processing: $folder / $base"
  
  abricate --db ehec_markers --minid 80 --mincov 80 "$f" > "abricate_out/$outname"
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
