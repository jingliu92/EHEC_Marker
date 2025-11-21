# EHEC_Marker

## 1. Genes interested
`eae` (intimin; include a few alleles if you care about subtyping, e.g., β, γ, ε). \
`stx1` and `stx2` (include major subtypes: stx1a/1c/1d; stx2a–2g to avoid false negatives). \
`espK, espV and espN` (T3SS effectors; often on OI-122 / related islands). \
`aggR`: Encodes the AggR protein, which is involved in the expression of several virulence determinants in both plasmids and chromosomes. The AggR protein is crucial for the adherence of EAEC to host cells and is associated with the pathogenesis of diarrhea in both children and adults. \
`subA`: Subtilase cytotoxin. \
`Z2098` (OI-122 marker; locus tag from the O157 Sakai island). \
`ureD` (urease accessory subunit—present in some clades; include a couple of E. coli hits to avoid over-stringency). \

<img width="585" height="361" alt="image" src="https://github.com/user-attachments/assets/6ebcdb81-8a39-4c00-b477-650429ab189d" />


## 2. Tasks
1. Identify all eae+ isolates
2. Find their linkage to stx, esp, ureD and z2098. Then create the linkage table.
3. Find aggR that are positive for stx or stx/subA
4. Find eae negative but stx+ isolates (STEC). STEC is not important for food safety. 

## 3. Bioinformatics Analysis

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
### 3.1 Download assemblies from NCBI
```
# Get ALL complete genomes for Escherichia coli 
./datasets download genome taxon "Escherichia coli" \
  --filename E.coli_complete.zip
```
<img width="1326" height="66" alt="image" src="https://github.com/user-attachments/assets/54e1ffe6-672e-4e96-bb2d-cc054760e37e" />
In total, 11,321 genome assemblies were retrieved and downloaded. The assemblies are organized into folders, each named according to the respective NCBI accession number: 
<img width="839" height="108" alt="image" src="https://github.com/user-attachments/assets/57a48e37-c570-4dd6-a83d-6e3bbc03f457" />

### 3.2 Search target genes in the downloaded assemblies
### 3.2.1 Option 1: Use Blast-----large dataset takes forever, very slow, requires high computing compacity
🌀 Step 1: Prepare your custom marker database
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
🌀 Step 2: Run BLAST for one genome assembly for testing purpose:
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

🌀 Step 3: Add filtering (optional)
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

🌀 Step 5 Check how many assemblies have hits? (1984)
```
cd blast_hits_only
ls -l . | wc -l
# 1984
```
<img width="879" height="49" alt="image" src="https://github.com/user-attachments/assets/d7c68fe3-3219-4edc-a5f5-ac21c1ddc0b7" />

🌀 Step 6 Add headers to each tsv file
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
🌀 Step 7: Summarize presence/absence (Python)

- Reads all .tsv files in a folder

- Filters for hits ≥80% identity and ≥80% coverage

- Outputs a presence/absence matrix

nano make_presence_matrix.py
```
#!/usr/bin/env python3
import argparse, glob, os, pandas as pd, re

def canonical_gene(name: str) -> str:
    n = str(name).strip().split()[0]
    l = n.lower()
    if l.startswith("stx1"): return "stx1"
    if l.startswith("stx2"): return "stx2"
    m = {"eae":"eae","espk":"espK","espv":"espV","z2098":"Z2098",
         "ured":"ureD","espn":"espN","suba":"subA","aggr":"aggR"}
    for k,v in m.items():
        if l.startswith(k): return v
    return n  # fallback

def read_blast_tsv(path: str) -> pd.DataFrame:
    # Try 11-column (with qlen)
    cols11 = ["qseqid","sseqid","pident","length","qlen",
              "qstart","qend","sstart","send","evalue","bitscore"]
    cols10 = ["qseqid","sseqid","pident","length",
              "qstart","qend","sstart","send","evalue","bitscore"]
    try:
        df = pd.read_csv(path, sep="\t")
        # If no header was actually present, try to coerce
        if "qseqid" not in df.columns:
            raise ValueError("No header detected")
    except Exception:
        # Headerless or mismatched: sniff by field count
        first = open(path).readline().rstrip("\n")
        nf = len(first.split("\t"))
        if nf == 11:
            df = pd.read_csv(path, sep="\t", header=None, names=cols11)
        elif nf == 10:
            df = pd.read_csv(path, sep="\t", header=None, names=cols10)
        else:
            # empty file or weird shape
            return pd.DataFrame(columns=cols11)
    # Ensure columns exist
    for c in cols11:
        if c not in df.columns: df[c] = pd.NA
    return df

def main():
    ap = argparse.ArgumentParser(description="Build presence/absence matrix from BLAST outfmt 6 TSVs.")
    ap.add_argument("--glob", default="blast_out/*_hits.tsv", help="Glob for BLAST TSVs (default: blast_out/*_hits.tsv)")
    ap.add_argument("--id", type=float, default=80.0, help="Min % identity to keep (default 80)")
    ap.add_argument("--cov", type=float, default=0.80, help="Min query coverage (length/qlen) to keep (default 0.80)")
    ap.add_argument("--outprefix", default="blast", help="Output prefix (default blast)")
    args = ap.parse_args()

    files = sorted(glob.glob(args.glob, recursive=True))
    if not files:
        print(f"No files matched {args.glob}")
        return

    presence_rows = []
    best_rows = []

    for f in files:
        sample = os.path.splitext(os.path.basename(f))[0]
        df = read_blast_tsv(f)

        if df.empty:
            presence_rows.append({"Sample": sample})
            continue

        # Compute coverage when qlen is available; if missing, assume already filtered upstream
        if "qlen" in df.columns and df["qlen"].notna().any():
            df["coverage"] = pd.to_numeric(df["length"], errors="coerce") / pd.to_numeric(df["qlen"], errors="coerce")
            df = df[(pd.to_numeric(df["pident"], errors="coerce") >= args.id) &
                    (pd.to_numeric(df["coverage"], errors="coerce") >= args.cov)]
        else:
            # no qlen -> keep as-is or filter by identity only
            df["coverage"] = 1.0
            df = df[pd.to_numeric(df["pident"], errors="coerce") >= args.id]

        if df.empty:
            presence_rows.append({"Sample": sample})
            continue

        df["Gene"] = df["qseqid"].apply(canonical_gene)

        # Best hit per gene: max identity, then coverage, then length
        best = (df.sort_values(["Gene","pident","coverage","length"], ascending=[True, False, False, False])
                  .groupby("Gene", as_index=False).first())

        # presence row for this sample
        row = {"Sample": sample}
        for g in best["Gene"].unique(): row[g] = 1
        presence_rows.append(row)

        # collect best hits rows
        for _, r in best.iterrows():
            best_rows.append({
                "Sample": sample,
                "Gene": r["Gene"],
                "Identity(%)": round(float(r["pident"]),2) if pd.notna(r["pident"]) else None,
                "Coverage": round(float(r["coverage"]),3) if pd.notna(r["coverage"]) else None,
                "Subject": r.get("sseqid",""),
                "AlignLen": int(r.get("length",0)) if pd.notna(r.get("length",0)) else None
            })

    mat = pd.DataFrame(presence_rows).fillna(0)
    # Order columns if present
    order = ["Sample","stx1","stx2","eae","espK","espV","Z2098","ureD","espN","subA","aggR"]
    cols = [c for c in order if c in mat.columns] + [c for c in mat.columns if c not in order]
    mat = mat[cols]
    for c in mat.columns:
        if c != "Sample": mat[c] = mat[c].astype(int)

    # Save outputs
    pa_path = f"{args.outprefix}_presence_absence.tsv"
    bh_path = f"{args.outprefix}_best_hits.tsv"
    mat.to_csv(pa_path, sep="\t", index=False)
    pd.DataFrame(best_rows).to_csv(bh_path, sep="\t", index=False)

    print(f"Wrote: {pa_path}")
    print(f"Wrote: {bh_path}")
    print("Tip: open the matrix with `column -t blast_presence_absence.tsv | less -S`")

if __name__ == "__main__":
    main()

```
🌀 Step 8: Run the script
```
python make_presence_matrix.py --glob "/home/jing/E.coli/blast_results/blast_hits_only/*.tsv"
```
#### Presence and absence matrix
<img width="819" height="352" alt="image" src="https://github.com/user-attachments/assets/5c40d019-5f9e-4321-a517-aee62f9cfe03" />

### 3.2 Option 2: Use package ABRicate.
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

🌀 Step 1:  Install `ABRicate` in Our Server
```
# Create a clean environment for abricate
conda create -n abricate -c bioconda -c conda-forge abricate

# Activate it
conda activate abricate

# Check that it works
abricate --version
```
🌀 Step 2: Make Our Custom Marker DB
Once ABRicate is installed, make a custom database for our interested markers:
```
mkdir -p /home/jing/miniconda3/envs/abricate/db/ehec_markers
cp ehec_markers.fna /home/jing/miniconda3/envs/abricate/db/ehec_markers/sequences
abricate --setupdb
abricate --list   # should now show ehec_markers
```
🌀 Step 3: ABRicate with a custom DB (threshold: --minid 90 --mincov 90)
```
# Create output directory
mkdir -p abricate_out

# Loop through all subfolders containing .fna files
for f in $(find . -type f -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")   # get folder name (e.g., GCA_037094545.1)
  base=$(basename "$f" .fna)             # get file name without extension
  outname="${folder}_${base}.tsv"        # combine both for unique output name
  echo "Processing: $folder / $base"
  
  abricate --db ehec_markers --minid 90 --mincov 90 "$f" > "abricate_out/$outname"
done

abricate --summary abricate_out/*.tsv > abricate_summary.tsv
echo "Wrote abricate_summary.tsv"
```
<img width="1070" height="242" alt="image" src="https://github.com/user-attachments/assets/e9992ef4-da90-4f1e-8faf-a9e52e31d78c" />

🌀 Step 4:  Make presence/absence matrix + classification

Save as summarize_from_abricate.py in the same folder as abricate_summary.tsv, then run: summarize_from_abricate.py
```
#!/usr/bin/env python3
import re
import pandas as pd

INPUT = "abricate_summary.tsv"
GENES = ["stx1","stx2","eae","espK","espV","Z2098","ureD","espN","subA","aggR"]

# Read the summary
df = pd.read_csv(INPUT, sep="\t", dtype=str)

# Normalize column names (strip spaces)
df.columns = [c.strip() for c in df.columns]

# Accept either '#FILE' or 'FILE'
file_col = None
for cand in ("#FILE", "FILE"):
    if cand in df.columns:
        file_col = cand
        break
if file_col is None:
    raise SystemExit("Error: neither '#FILE' nor 'FILE' column found. Check abricate_summary.tsv header.")

df.rename(columns={file_col: "Sample"}, inplace=True)

# Keep only columns we care about and that actually exist
present_cols = [g for g in GENES if g in df.columns]
if not present_cols:
    raise SystemExit(f"No expected gene columns found. Got columns: {list(df.columns)}")

mat = df[["Sample"] + present_cols].copy()

# Convert ABRicate summary entries to presence/absence:
# - ABRicate puts '.' for none
# - If present, it may list one or more values like '96.00;100.00'
# We'll call presence=1 if any digit appears in the cell.
def to_presence(cell: str) -> int:
    if cell is None:
        return 0
    cell = str(cell).strip()
    return 1 if re.search(r"\d", cell) else 0

for g in present_cols:
    mat[g] = mat[g].apply(to_presence).astype(int)

# Classification
def classify(row) -> str:
    stx = int(row.get("stx1", 0)) or int(row.get("stx2", 0))
    eae = int(row.get("eae", 0))
    support = sum(int(row.get(k, 0)) for k in ("espK","espV","Z2098","ureD") if k in row.index)
    if stx and (eae or support >= 2): return "EHEC"
    if (not stx) and eae and support >= 1: return "Typical EPEC"
    if (not stx) and eae and support == 0: return "Atypical EPEC"
    if stx and not eae: return "STEC (LEE−)"
    return "Commensal/Other"

cls = mat.copy()
cls["CLASS"] = cls.apply(classify, axis=1)

# Write outputs
mat.to_csv("abricate_presence_absence.tsv", sep="\t", index=False)
cls.to_csv("abricate_classification.tsv", sep="\t", index=False)
print("✅ Wrote: abricate_presence_absence.tsv and abricate_classification.tsv")
```
🌀 Step 5:  Remove assemblies that has no genes interested
```
awk 'NR==1{print; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum>0) print}' abricate_presence_absence.tsv > abricate_presence_absence_filtered.tsv
```
🌀 Step 6 Check how many assemblies have hits? 

(2390 if threshold: --minid 80 --mincov 80 )
```
wc -l abricate_presence_absence_filtered.tsv

# 2390
```
<img width="837" height="41" alt="image" src="https://github.com/user-attachments/assets/3b084a46-cd6c-44ad-972a-039052be807c" />

(1966 if threshold: --minid 90 --mincov 90 )

<img width="735" height="43" alt="image" src="https://github.com/user-attachments/assets/6cc6daa9-6cb4-4e34-a9ae-6cbd6220c06f" />

🌀 Step 7 Statistical analysis

nano abricate_statistics_summary.py
```
#!/usr/bin/env python3
import pandas as pd

INPUT  = "abricate_presence_absence_filtered.tsv"
OUTPUT = "abricate_statistics_summary.csv"

GENES = ["stx1","stx2","eae","espK","espV","espN"]

df = pd.read_csv(INPUT, sep="\t")

# Ensure all genes exist and are numeric
for g in GENES:
    if g not in df.columns:
        df[g] = 0
    df[g] = pd.to_numeric(df[g], errors="coerce").fillna(0).astype(int)

n = len(df)
def count(mask): return int(mask.sum())
def pct(mask):   return round(100 * count(mask) / n, 2) if n else 0.0

rows = []
def add_row(stat, genotype, mask):
    rows.append([stat, genotype, count(mask), pct(mask)])

# ---------------- BASIC ----------------
add_row("Total isolates", "", pd.Series([True]*n))
add_row("stx1/stx2(+) (this means Stx1 OR Stx2)", "", (df.stx1==1)|(df.stx2==1))
add_row("stx1(+)", "", df.stx1==1)
add_row("stx2(+)", "", df.stx2==1)
# NEW: double-positive
add_row("stx1(+)&stx2(+) (both stx1 and stx2 present)", "", (df.stx1==1)&(df.stx2==1))

add_row("stx1+stx2(-) (this means both stx1 and stx2 are absent)", "", (df.stx1==0)&(df.stx2==0))
add_row("eae(+)", "", df.eae==1)
add_row("eae(-)", "", df.eae==0)
add_row("espK(+)", "", df.espK==1)
add_row("espK/V(+) (espK AND espV positive)", "", (df.espK==1)&(df.espV==1))
add_row("espV(+)", "", df.espV==1)
add_row("espK AND/OR espV", "", (df.espK==1)|(df.espV==1))
add_row("espN(+)", "", df.espN==1)
add_row("espK/N (+) (espK AND espN positive)", "", (df.espK==1)&(df.espN==1))
add_row("espK AND/OR espN", "", (df.espK==1)|(df.espN==1))

# ---------------- NOT STEC/EPEC/EHEC ----------------
add_row("stx1+stx2(-) AND eae(-)", "NOT STEC/EPEC/EHEC", (df.stx1==0)&(df.stx2==0)&(df.eae==0))
add_row("stx1+stx2(-) AND eae(-) AND espK(+)", "NOT STEC/EPEC/EHEC", (df.stx1==0)&(df.stx2==0)&(df.eae==0)&(df.espK==1))
add_row("stx1+stx2(-) AND eae(-) AND espN(+)", "NOT STEC/EPEC/EHEC", (df.stx1==0)&(df.stx2==0)&(df.eae==0)&(df.espN==1))
add_row("stx1+stx2(-) AND eae(-) AND espV(+)", "NOT STEC/EPEC/EHEC", (df.stx1==0)&(df.stx2==0)&(df.eae==0)&(df.espV==1))

# ---------------- EPEC ----------------
add_row("stx1+stx2(-) AND eae(+)", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1))
add_row("stx1+stx2(-) AND eae(+) AND espK(+)", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&(df.espK==1))
add_row("stx1+stx2(-) AND eae(+) AND espV(+)", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&(df.espV==1))
add_row("stx1+stx2(-) AND eae(+) AND espN(+)", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&(df.espN==1))
add_row("stx1+stx2(-) AND eae(+) AND {espK(+) AND/OR espV(+)}", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&((df.espK==1)|(df.espV==1)))
add_row("stx1+stx2(-) AND eae(+) AND {espK(+) AND espV(+)}", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&(df.espK==1)&(df.espV==1))
add_row("stx1+stx2(-) AND eae(+) AND {espK(+) AND/OR espN(+)}", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&((df.espK==1)|(df.espN==1)))
add_row("stx1+stx2(-) AND eae(+) AND {espK(+) AND espN(+)}", "EPEC cohorts", (df.stx1==0)&(df.stx2==0)&(df.eae==1)&(df.espK==1)&(df.espN==1))

# ---------------- EHEC ----------------
add_row("stx1+stx2(+) AND eae(+)", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1))
add_row("stx1+stx2(+) AND eae(+) AND espK(+)", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&(df.espK==1))
add_row("stx1+stx2(+) AND eae(+) AND espV(+)", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&(df.espV==1))
add_row("stx1+stx2(+) AND eae(+) AND espN(+)", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&(df.espN==1))
add_row("stx1+stx2(+) AND eae(+) AND {espK(+) AND/OR espV(+)}", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&((df.espK==1)|(df.espV==1)))
add_row("stx1+stx2(+) AND eae(+) AND {espK(+) AND espV(+)}", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&(df.espK==1)&(df.espV==1))
add_row("stx1+stx2(+) AND eae(+) AND {espK(+) AND/OR espN(+)}", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&((df.espK==1)|(df.espN==1)))
add_row("stx1+stx2(+) AND eae(+) AND {espK(+) AND espN(+)}", "EHEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==1)&(df.espK==1)&(df.espN==1))

# ---------------- STEC ----------------
add_row("stx1+stx2(+) AND eae(-)", "Pathovars", ((df.stx1==1)|(df.stx2==1))&(df.eae==0))
add_row("stx1+stx2(+) AND eae(-) AND espK(+)", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&(df.espK==1))
add_row("stx1+stx2(+) AND eae(-) AND espV(+)", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&(df.espV==1))
add_row("stx1+stx2(+) AND eae(-) AND espN(+)", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&(df.espN==1))
add_row("stx1+stx2(+) AND eae(-) AND {espK(+) AND/OR espV(+)}", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&((df.espK==1)|(df.espV==1)))
add_row("stx1+stx2(+) AND eae(-) AND {espK(+) AND espV(+)}", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&(df.espK==1)&(df.espV==1))
add_row("stx1+stx2(+) AND eae(-) AND {espK(+) AND/OR espN(+)}", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&((df.espK==1)|(df.espN==1)))
add_row("stx1+stx2(+) AND eae(-) AND {espK(+) AND espN(+)}", "STEC cohorts", ((df.stx1==1)|(df.stx2==1))&(df.eae==0)&(df.espK==1)&(df.espN==1))

# ---------------- eae(+) triple-marker patterns ----------------
E = df.eae==1
K,V,N = (df.espK==1),(df.espV==1),(df.espN==1)
add_row("eae(+) & espK only", "", E &  K & ~V & ~N)
add_row("eae(+) & espV only", "", E & ~K &  V & ~N)
add_row("eae(+) & espN only", "", E & ~K & ~V &  N)
add_row("eae(+) & espK+espV", "", E &  K &  V & ~N)
add_row("eae(+) & espK+espN", "", E &  K & ~V &  N)
add_row("eae(+) & espV+espN", "", E & ~K &  V &  N)
add_row("eae(+) & espK+espV+espN", "", E &  K &  V &  N)
add_row("eae(+) & any(espK/espV/espN)", "", E & (K|V|N))
add_row("eae(+) & none(espK/espV/espN)", "", E & ~(K|V|N))

out = pd.DataFrame(rows, columns=["Statistics for Abricate","Genotype","Count","Percent(%)"])
out.to_csv(OUTPUT, index=False)
print(f"✅ Wrote summary to {OUTPUT}")
```
<img width="645" height="756" alt="image" src="https://github.com/user-attachments/assets/d86474d7-c2a2-47a8-ae6b-5fcf8b920491" />

🌀 Step 8 Venn Diagram
- draws two 3-set Venns (espK/espV/espN) in the two strata **eae+ & (stx1 OR stx2)+ and 2) eae+ & (stx1 AND stx2)−**

- shows counts + percentages inside each Venn region

- also makes UpSet plots for the same two strata (great when labels get crowded)

- writes a tidy TSV with counts and % per region

Script: `venn_eae_strata_with_percents_upset.py`
```
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

# pip install matplotlib-venn upsetplot
from matplotlib_venn import venn3
from upsetplot import UpSet, from_indicators

INPUT = "abricate_presence_absence_filtered.tsv"
OUT_VENN_STXPOS = "venn_eaepos_stxpos.png"
OUT_VENN_STXNEG = "venn_eaepos_stxneg.png"
OUT_UPSET_STXPOS = "upset_eaepos_stxpos.png"
OUT_UPSET_STXNEG = "upset_eaepos_stxneg.png"
OUT_TSV = "venn_eae_strata_counts_percents.tsv"

df = pd.read_csv(INPUT, sep="\t")

# Ensure 0/1 ints for the fields we need
for g in ["stx1","stx2","eae","espK","espV","espN"]:
    df[g] = pd.to_numeric(df.get(g, 0), errors="coerce").fillna(0).astype(int)

# Strata
eae_pos = df["eae"] == 1
stx_pos = (df["stx1"] == 1) | (df["stx2"] == 1)
stx_neg = (df["stx1"] == 0) & (df["stx2"] == 0)

def venn_counts(sub):
    K = sub["espK"] == 1
    V = sub["espV"] == 1
    N = sub["espN"] == 1
    # order for venn3: 100,010,110,001,101,011,111
    a_only = int(( K & ~V & ~N).sum())
    b_only = int((~K &  V & ~N).sum())
    ab     = int(( K &  V & ~N).sum())
    c_only = int((~K & ~V &  N).sum())
    ac     = int(( K & ~V &  N).sum())
    bc     = int((~K &  V &  N).sum())
    abc    = int(( K &  V &  N).sum())
    total  = len(sub)
    return (a_only, b_only, ab, c_only, ac, bc, abc), total

def label_with_counts_percents(v, total, subsets):
    # Map venn region IDs to our tuple order
    id_map = {
        "100": 0, "010": 1, "110": 2,
        "001": 3, "101": 4, "011": 5,
        "111": 6
    }
    for vid, idx in id_map.items():
        lab = v.get_label_by_id(vid)
        if lab is None:
            continue
        cnt = subsets[idx]
        if total > 0:
            pct = 100.0 * cnt / total
            lab.set_text(f"{cnt}\n({pct:.1f}%)")
        else:
            lab.set_text("0\n(0.0%)")

def plot_venn(sub, title, outfile):
    subsets, total = venn_counts(sub)
    plt.figure(figsize=(6.5,6.5))
    v = venn3(subsets=subsets, set_labels=("espK","espV","espN"))
    label_with_counts_percents(v, total, subsets)
    # Title with N
    plt.title(f"{title}\nN={total}")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()
    return subsets, total

def make_upset(sub, title, outfile):
    # Build boolean indicator frame for K,V,N
    ind = sub[["espK","espV","espN"]].astype(bool)
    u = from_indicators(ind)
    plt.figure(figsize=(8,6))
    UpSet(u, sort_by="cardinality").plot()
    plt.suptitle(title)
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close()

def table_rows(name, subsets, total):
    # order: 100,010,110,001,101,011,111
    keys = ["K_only","V_only","K∩V","N_only","K∩N","V∩N","K∩V∩N"]
    rows = []
    for k, c in zip(keys, subsets):
        pct = round(100.0 * c / total, 2) if total else 0.0
        rows.append([name, k, c, pct, total])
    # add any/none
    any_kvn = sum(subsets)  # any of K/V/N
    none_kvn = total - any_kvn
    rows.append([name, "any(K/V/N)", any_kvn, round(100.0*any_kvn/total,2) if total else 0.0, total])
    rows.append([name, "none(K/V/N)", none_kvn, round(100.0*none_kvn/total,2) if total else 0.0, total])
    return rows

# 1) eae+ & stx+
sub_pos = df[eae_pos & stx_pos].copy()
s1, n1 = plot_venn(sub_pos, "eae+ AND (stx1 OR stx2)+", OUT_VENN_STXPOS)
make_upset(sub_pos, "UpSet: eae+ AND (stx1 OR stx2)+", OUT_UPSET_STXPOS)

# 2) eae+ & stx−
sub_neg = df[eae_pos & stx_neg].copy()
s2, n2 = plot_venn(sub_neg, "eae+ AND (stx1 AND stx2)−", OUT_VENN_STXNEG)
make_upset(sub_neg, "UpSet: eae+ AND (stx1 AND stx2)−", OUT_UPSET_STXNEG)

# Write table
rows = []
rows += table_rows("eae+ & stx+", s1, n1)
rows += table_rows("eae+ & stx−", s2, n2)
pd.DataFrame(rows, columns=["Stratum","Region","Count","Percent(%)","N"]).to_csv(OUT_TSV, index=False)

print("✅ Wrote:")
print(f" - {OUT_VENN_STXPOS}")
print(f" - {OUT_VENN_STXNEG}")
print(f" - {OUT_UPSET_STXPOS}")
print(f" - {OUT_UPSET_STXNEG}")
print(f" - {OUT_TSV}")
```

## 4. Comparison of two methods
🧬 1️⃣ Methodological difference
| Tool                          | How it works                                                                                                                                                                                                                                   | Typical behavior                                                                                                                                 |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **ABRicate**                  | Uses **BLASTn under the hood**, but applies its own filters (`--minid`, `--mincov`), and **scans all contigs** for any partial or full matches. It records every region that passes the threshold.                                             | Multiple hits per genome if a gene occurs in several loci, plasmids, or duplicated regions.                                                      |
| **Your manual BLAST command** | You used `blastn -task blastn-short -query markers -subject assembly` with an **`awk` filter** `(pident >= 80 && cov >= 80)` that outputs only alignments satisfying both conditions — **but it stops at the sequence level, not gene level.** | Misses hits if the query covers multiple disjoint regions (e.g., gene split across contigs), or if coverage falls slightly below your threshold. |

🧩 2️⃣ Typical reasons ABRicate reports more hits than raw BLAST
| Cause                                 | Explanation                                                                                                                                                                                           |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Gene fragmentation**                | ABRicate merges multiple high-identity fragments to call a gene “present” if total coverage meets threshold; plain BLAST outputs individual fragments separately (and your awk filter may drop them). |
| **Multiple query matches per genome** | If a marker gene occurs on plasmids or in duplicated loci, ABRicate counts all; you may have been filtering only the best local alignment per query.                                                  |
| **Different BLAST task**              | ABRicate uses **`blastn`**, while you used **`blastn-short`** — that’s optimized for <50 bp sequences. For full-length genes, it can miss long hits.                                                  |
| **Coverage calculation**              | ABRicate calculates coverage relative to *gene length*, not *query alignment length* (and it merges hits before calculating coverage). Your awk filter calculated coverage per fragment.              |
| **Database indexing differences**     | ABRicate indexes marker databases with makeblastdb, sometimes allowing slightly looser word sizes; your raw `blastn` call may be stricter by default.                                                 |
| **Parsing thresholds**                | ABRicate’s internal `--mincov` and `--minid` compare float values differently (≥ vs >); minor, but can change borderline hits.                                                                        |

### Compare ABRicate vs BLAST per genome & gene (show agreements/disagreements, and capture BLAST best-hit details)
1️⃣ Save as compare_abricate_blast.py in the folder with your results.
```
#!/usr/bin/env python3
import argparse, glob, os, re
import pandas as pd

GENES_DEFAULT = ["stx1","stx2","eae","espK","espV","Z2098","ureD","espN","subA","aggR"]

def parse_args():
    ap = argparse.ArgumentParser(description="Compare ABRicate vs BLAST gene calls.")
    ap.add_argument("--abricate", required=True, help="Path to abricate_summary.tsv")
    ap.add_argument("--blast_glob", required=True, help="Glob for BLAST TSVs (e.g., blast_out/*_genomic_hits.tsv)")
    ap.add_argument("--id", type=float, default=80.0, help="Min % identity to count BLAST presence (default 80)")
    ap.add_argument("--cov", type=float, default=0.80, help="Min query coverage (length/qlen) for BLAST presence (default 0.80)")
    ap.add_argument("--genes", nargs="*", default=GENES_DEFAULT, help="Genes to compare (default common EHEC/EPEC markers)")
    ap.add_argument("--outprefix", default="abricate_vs_blast", help="Output prefix (default abricate_vs_blast)")
    ap.add_argument("--all", action="store_true", help="Output all rows (not only discrepancies) in the long table")
    return ap.parse_args()

def extract_accession(s: str) -> str:
    """Extract first GCA_/GCF_ accession from a path/filename."""
    if not isinstance(s, str):
        return str(s)
    m = re.search(r'(G[CF]A_\d+\.\d+)', s)
    return m.group(1) if m else os.path.splitext(os.path.basename(s))[0]

def canonical_gene(name: str) -> str:
    n = str(name).strip().split()[0]
    l = n.lower()
    if l.startswith("stx1"): return "stx1"
    if l.startswith("stx2"): return "stx2"
    m = {"eae":"eae","espk":"espK","espv":"espV","z2098":"Z2098",
         "ured":"ureD","espn":"espN","suba":"subA","aggr":"aggR"}
    for k,v in m.items():
        if l.startswith(k): return v
    return n

def read_abricate_presence(path: str, genes):
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    file_col = "#FILE" if "#FILE" in df.columns else ("FILE" if "FILE" in df.columns else None)
    if file_col is None:
        raise SystemExit("Could not find '#FILE' or 'FILE' column in abricate summary.")
    df["SampleID"] = df[file_col].apply(extract_accession)

    present_cols = [g for g in genes if g in df.columns]
    if not present_cols:
        raise SystemExit(f"No expected gene columns found in abricate summary. Got: {list(df.columns)}")

    # Presence = any digit in the cell (ABRicate uses '.' for none; numbers/semicolon lists for hits)
    def any_digit(x): 
        return 1 if isinstance(x, str) and re.search(r"\d", x) else 0
    pres = df[["SampleID"] + present_cols].copy()
    for g in present_cols:
        pres[g] = pres[g].apply(any_digit).astype(int)

    # If some genes are missing from the summary, add 0 columns
    for g in genes:
        if g not in pres.columns:
            pres[g] = 0
    return pres[["SampleID"] + genes].drop_duplicates()

def read_blast_presence(globpat: str, genes, id_cut: float, cov_cut: float):
    files = sorted(glob.glob(globpat))
    if not files:
        raise SystemExit(f"No files matched {globpat}")
    rows = []
    details = []
    for f in files:
        sid = extract_accession(os.path.basename(f))
        # Try to read with header
        try:
            df = pd.read_csv(f, sep="\t")
            if "qseqid" not in df.columns:
                raise Exception()
        except Exception:
            # Fallback to 11/10 col formats
            cols11 = ["qseqid","sseqid","pident","length","qlen","qstart","qend","sstart","send","evalue","bitscore"]
            cols10 = ["qseqid","sseqid","pident","length","qstart","qend","sstart","send","evalue","bitscore"]
            first = open(f).readline().rstrip("\n")
            nf = len(first.split("\t")) if first else 0
            if nf == 11:
                df = pd.read_csv(f, sep="\t", header=None, names=cols11)
            elif nf == 10:
                df = pd.read_csv(f, sep="\t", header=None, names=cols10)
            else:
                df = pd.DataFrame(columns=cols11)

        if df.empty:
            rows.append({"SampleID": sid})
            continue

        # compute coverage if qlen exists
        if "qlen" in df.columns:
            df["coverage"] = pd.to_numeric(df["length"], errors="coerce") / pd.to_numeric(df["qlen"], errors="coerce")
        else:
            df["coverage"] = 1.0

        df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
        df = df[df["pident"].notna()]
        keep = (df["pident"] >= id_cut) & (df["coverage"] >= cov_cut)
        df = df[keep]
        if df.empty:
            rows.append({"SampleID": sid})
            continue

        df["Gene"] = df["qseqid"].apply(canonical_gene)
        # best per gene
        best = (df.sort_values(["Gene","pident","coverage","length"], ascending=[True,False,False,False])
                  .groupby("Gene", as_index=False).first())
        row = {"SampleID": sid}
        for g in best["Gene"].unique():
            row[g] = 1
        rows.append(row)

        for _, r in best.iterrows():
            details.append({
                "SampleID": sid, "Gene": r["Gene"],
                "BLAST_best_identity(%)": round(float(r["pident"]),2),
                "BLAST_best_coverage": round(float(r["coverage"]),3),
                "Subject": r.get("sseqid",""),
                "AlignLen": int(r.get("length",0)) if pd.notna(r.get("length",0)) else None
            })

    pres = pd.DataFrame(rows).fillna(0)
    for g in genes:
        if g not in pres.columns:
            pres[g] = 0
    # Ensure ints
    for c in pres.columns:
        if c != "SampleID": pres[c] = pres[c].astype(int)
    return pres[["SampleID"] + genes], pd.DataFrame(details)

def main():
    args = parse_args()
    genes = args.genes

    ab = read_abricate_presence(args.abricate, genes)
    bl, bl_details = read_blast_presence(args.blast_glob, genes, args.id, args.cov)

    # Melt to long format for comparison
    ab_long = ab.melt(id_vars="SampleID", var_name="Gene", value_name="ABRicate")
    bl_long = bl.melt(id_vars="SampleID", var_name="Gene", value_name="BLAST")

    merged = pd.merge(ab_long, bl_long, on=["SampleID","Gene"], how="outer").fillna(0)
    merged["ABRicate"] = merged["ABRicate"].astype(int)
    merged["BLAST"] = merged["BLAST"].astype(int)

    # Attach BLAST details (best identity/coverage) if available
    merged = merged.merge(bl_details, on=["SampleID","Gene"], how="left")

    def status(row):
        a, b = row["ABRicate"], row["BLAST"]
        if a==1 and b==1: return "Agree+"
        if a==0 and b==0: return "Agree−"
        if a==1 and b==0: return "ABRicate_only"
        if a==0 and b==1: return "BLAST_only"
        return "NA"
    merged["Status"] = merged.apply(status, axis=1)

    # Discrepancies table
    discrep = merged if args.all else merged[merged["Status"].isin(["ABRicate_only","BLAST_only"])]
    discrep = discrep.sort_values(["Gene","Status","SampleID"])
    discrep.to_csv(f"{args.outprefix}_discrepancies.tsv", sep="\t", index=False)

    # Summary counts by gene
    summary = (merged
               .groupby(["Gene","Status"]).size()
               .reset_index(name="Count")
               .pivot(index="Gene", columns="Status", values="Count")
               .fillna(0).astype(int)
               .reset_index())
    summary.to_csv(f"{args.outprefix}_agreement_summary.tsv", sep="\t", index=False)

    # Wide matrix with both calls side-by-side
    wide_parts = []
    for src in ["ABRicate","BLAST"]:
        w = merged.pivot_table(index="SampleID", columns="Gene", values=src, fill_value=0, aggfunc="max")
        w.columns = [f"{src}_{c}" for c in w.columns]
        wide_parts.append(w)
    wide = pd.concat(wide_parts, axis=1).reset_index()
    wide.to_csv(f"{args.outprefix}_wide.tsv", sep="\t", index=False)

    print("✅ Wrote:")
    print(f" - {args.outprefix}_discrepancies.tsv")
    print(f" - {args.outprefix}_agreement_summary.tsv")
    print(f" - {args.outprefix}_wide.tsv")
    print("\nTip: open the summary to see which genes disagree most between methods.")

if __name__ == "__main__":
    main()
```
2️⃣ Run the script
```
python compare_abricate_blast.py \
  --abricate abricate_summary.tsv \
  --blast_glob "/home/jing/E.coli/blast_out/*_genomic_hits.tsv" \
  --id 80 --cov 0.80 \
  --outprefix abricate_vs_blast
```
#### Outputs
- abricate_vs_blast_discrepancies.tsv → long table of per-sample per-gene comparisons + BLAST best %ID/coverage (only rows where methods differ unless you add --all).

- abricate_vs_blast_agreement_summary.tsv → counts of Agree+/Agree−/ABRicate_only/BLAST_only per gene.

- abricate_vs_blast_wide.tsv → wide matrix: for each sample/gene, ABR (0/1) and BLAST (0/1) side-by-side.

#### Check results

_agreement_summary.tsv → if you see many ABRicate_only for a gene, those are usually partial/fragmented hits that ABRicate merges (or slightly below your BLAST thresholds).

_discrepancies.tsv → check BLAST %identity and coverage columns to see if BLAST is just under cutoff; try rerunning comparison with --id 75 --cov 0.6 to see sensitivity.

_wide.tsv → handy for quick greps: grep -E "GCA_.*\t.*_espK"

#### Any gene confirmed by both
```
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    if($i ~ /^ABRicate_/) { abr[i]=$i }
  }
  # map matching BLAST columns
  for(i in abr){
    gene = substr(abr[i], 10)      # strip "ABRicate_"
    for(j=1;j<=NF;j++) if($j=="BLAST_" gene) bl[j]=i
  }
  next
}
{
  ok=0
  for(b in bl){
    a = bl[b]              # ABRicate col index
    if($(a)==1 && $(b)==1){ ok=1; break }
  }
  if(ok) print $1
}' abricate_vs_blast_wide.tsv | sort -u > both_methods_ANY.txt

```
#### Specific gene confirmed by both (example: espK)
```
awk -F'\t' '
NR==1{for(i=1;i<=NF;i++){if($i=="ABRicate_espK") a=i; if($i=="BLAST_espK") b=i} next}
($a==1 && $b==1){print $1}
' abricate_vs_blast_wide.tsv | sort -u > both_methods_espK.txt
```
#### All target genes confirmed by both
```
genes="eae stx1 stx2 espK espV Z2098 ureD"
awk -F'\t' -v genes="$genes" '
NR==1{
  n=split(genes,G," ");
  for(i=1;i<=n;i++){
    for(j=1;j<=NF;j++){ if($j=="ABRicate_" G[i]) A[i]=j }
    for(j=1;j<=NF;j++){ if($j=="BLAST_"    G[i]) B[i]=j }
  }
  next
}
{
  ok=1
  for(i=1;i<=n;i++){
    if($(A[i])!=1 || $(B[i])!=1){ ok=0; break }
  }
  if(ok) print $1
}' abricate_vs_blast_wide.tsv | sort -u > both_methods_ALL_TARGETS.txt

```
#### Get the presence and absence table that confirmed by both methods
```
{ head -n1 abricate_vs_blast_wide.tsv; awk 'NR==FNR{a[$1]=1; next} ($1 in a)' both_methods_ALL_TARGETS.txt abricate_vs_blast_wide.tsv | tail -n +2; } > both_methods_ALL_TARGETS_table.tsv

```
#### Make a co-occurrence table for every possible combination of markers found in your both_methods_ALL_TARGETS_table.tsv.
nano make_all_marker_combinations.py
```
#!/usr/bin/env python3
import pandas as pd
import itertools

INPUT = "both_methods_ALL_TARGETS_table.tsv"
OUTPUT = "all_marker_combinations.tsv"

# --- Load data ---
df = pd.read_csv(INPUT, sep="\t")

# Identify all genes (based on ABRicate_ and BLAST_ columns)
genes = sorted(set(c.split("_")[1] for c in df.columns if c.startswith("ABRicate_")))

# Merge ABRicate and BLAST results (1 only if both are 1)
for g in genes:
    a, b = f"ABRicate_{g}", f"BLAST_{g}"
    df[g] = ((df.get(a, 0) == 1) & (df.get(b, 0) == 1)).astype(int)

# --- Compute combination frequencies ---
n = len(df)
out = []

# Iterate through all possible combination sizes (1..N)
for r in range(1, len(genes) + 1):
    for combo in itertools.combinations(genes, r):
        # Check rows where all genes in combo = 1
        mask = df[list(combo)].all(axis=1)
        pct = 100 * mask.sum() / n
        out.append({
            "Combination": "/".join(combo),
            "Num_Genes": r,
            "Count": int(mask.sum()),
            "Frequency(%)": round(pct, 2)
        })

# --- Sort for readability ---
tbl = pd.DataFrame(out).sort_values(["Num_Genes", "Frequency(%)"], ascending=[True, False])

# Save output
tbl.to_csv(OUTPUT, sep="\t", index=False)
print(f"✅ Wrote {OUTPUT} with {len(tbl)} combinations across {len(genes)} genes.")
print(tbl.head(15).to_string(index=False))
```
```
python make_all_marker_combinations.py
```
From this step we will get the co-occurrence table for every possible combination of markers: all_marker_combinations.tsv

############################################## NEXT STEPS #########################################################################
# 1. Test the distribution of eae subtypes in EHEC and EPEC
## 1.1 Seacrh NCBI for the eae subtypes sequence
intimin_subtypes.fasta
```
>Gamma
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAGACAACTGCAGTAGCAAATGGTAAGGATGCTATTAAATATACTGTAAAAGTTATGAAAAACGGTCAGCCAGTTAATAATCAATCCGTTACATTCTCAACAAACTTTGGGATGTTCAACGGTAAGTCTCAAACGCAAGCAACCACGGGAAATGATGGTCGTGCGACGATAACACTAACTTCCAGTTCCGCCGGTAAAGCGACTGTTAGTGCGACAGTCAGTGATGGGGCTGAGGTTAAAGCGACTGAGGTCACTTTTTTTGATGAACTGAAAATTGACAACAAGGTTGATATTATTGGTAACAATGTCAGAGGCGAGTTGCCTAATATTTGGCTGCAATATGGTCAGTTTAAACTGAAAGCAAGCGGTGGTGATGGTACATATTCATGGTATTCAGAAAATACCAGTATCGCGACTGTCGATGCATCAGGGAAAGTCACTTTGAATGGTAAAGGCAGTGTCGTAATTAAAGCCACATCTGGTGATAAGCAAACAGTAAGTTACACTATAAAAGCACCGTCGTATATGATAAAAGTGGATAAGCAAGCCTATTATGCTGATGCTATGTCCATTTGCAAAAATTTATTACCATCCACACAGACGGTATTGTCAGATATTTATGACTCATGGGGGGCTGCAAATAAATATAGCCATTATAGTTCTATGAACTCAATAACTGCTTGGATTAAACAGACATCTAGTGAGCAGCGTTCTGGAGTATCAAGCACTTATAACCTAATAACACAAAACCCTCTTCCTGGGGTTAATGTTAATACTCCAAATGTCTATGCGGTTTGTGTAGAATAA
>Sigma
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAGACAACTGCAGTAGCAAATGGTCAGGATGCTGTCACATACACTGTTAAAGTGATGAAAGATGGTAAGCCAGTGCAGGGACACTCTGTTGCATTCTCAACAAACTTTGGGATGTTCAACGGTAAGGCTCAAGCACAAACTGCAATCACGGGAAGTGATGGTCGTGCGACGATTACACTAACTTCCAATTCCGCCGGTAAAGCGACTGTTAGTGCAGCAGTTAGTGGGGGAACTGAGGTTAAAGCGACTGAGGTTACTTTTTTTGATGAACTGAAAATTGACAATAAGGTTGATATTATTGGTAACAATGTCAGTGGCGAGTTGCCTAATATCTGGTTGCAATATGGTCAGTTCAAACTGAAGGCAAGTGGAGGTAACGGTATATATTCATGGTATTCAGAAAATACCAGTATCGCGACTGTTGATGCATCGGGGAAAGTCACCTTGAATGGTAAAGGCAGTGTCGTAATTAAAGCCACATCTGGTGATAAGCAAACAGTAAGTTACACTATTAAAGCACCGTCTTATATGATAAGAGTGGATAATAAAGCCAATTATGCTAGCGCTACGGCCATTTGTAAAAATTCTTTACCATCCTCACAGAAGGTATTAGCAGACATATTTAATTCATGGGGAGCGGCAAACAAATATGGTCATTATGGTTCTATGAACTCAATACCGGCCTGGATTAAACAAACCGAAAGTGATAAAAATTCTGGAGTATCAACCACTTATGATTTAATAAGACAAAATCCTCGTTCTGATGTTAATGTAGACACTCTAAATGTCTATGCCGTTTGTGTAGAATAA
>Mu
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAAACAACGGCAGTAGCAAATGGTCAGGATGCTATTACATACACTGTTAAAGTGATGAAAGATGGTCAGCCAGTGCAGGGACACTCCGTTACATTCTCAACAAACTTTGGGATATTCAACGGTAAGTCTCAGACACAAAATGCAACCACGGGAAATGATGGTCGTGCGACGATAACGCTGACTTCCAGTTCCGCAGGTAAAGCGACTGTTAGTGCGACAGTTAGTGGTGGAACTGATGTTAAAGCGACTGAGGTTACTTTTTTTGATGAACTGAAAATTGACAACAAGGTTGACATTATTGGTAACAATGTCAGTGGCGAATTGCCTGATATCTGGTTGCAATATGGTCAGTTTAAATTGAAGGCAAGTGGTGGTAACGGTACATATTCATGGTATTCAGAAAATACCAGTATCGCGACTGTTGATGCATCGGGGAAAGTCACCTTGAATGGTAAAGGCAGTGTCGTAATTAAAGCCACTTCTGGTGATAAGCAAACAGTAAGCTACACTGTTAAAACTCCAAAATATATGATAAAGGTAGGCCAAAAAGCTTATTATGATGAGTCCATGACTATTTGTAAAGGCTCTCTGCCATCCTCACAGACTGTATTATCTGACATCTTTAACGCATGGGGGGCAGCAAATAAATATAGCCATTATAACACTATGAGCTCAATAACGGCTTGGATTAAACAAACTGAAAAAGATAAAGCATCTGGAGTATCAAGTACTTATAACTTAATAACACAAAACCCTCTTCCCGGGGTTAACGTAAAAACTCCAAATGTCTATGCGGTTTGTGTAGAATAA
>THETA
ACCAAGGCCAGTATTACTGAGATTAAGGCTGATAAGACAACTGCAGTAGCAAATGGTAATGATGCTGTTACATACACTGTTAAAGTGATGAAAGAGGGTCAGCCAGTGCAGGGACACTCCGTTGCATTCACAACAAACTTTGGGATGTTCAACGGTAAGTCTCAGACGCAAAATGCGACCACGGGAAGTGATGGTCGTGCGACGATAACACTGACTTCCAGTTCCGCAGGTAAAGCGACTGTTAGTGCGACTGTTAGTGGTGGGAATGATGTTAAAGCACCTGAGGTTACATTTTTTGATGGACTGAAAATTGACAACAAGGTTGATATTCTTGGTAAGAACGTTACTGGTGACTTACCTAATATCTGGTTGCAATATGGTCAGTTTAAACTGAAGGTAAGCGGTGGTAATGGTACATATTCATGGCATTCAGAGAATACCAATATTGCGACTGTTGATGAATCAGGGAAAGTAACCTTGAAAGGAAAAGGTACTGCAGTAATTAATGTTACATCTGGTGATAAGCAAACAGTAAGCTACACTATTAAAGCTCCGAATTATATGATAAGAGTGGGTAATAAAGCCAGTTATGCAAATGCTATGTCCTTTTGTGGAAATTTATTACCATCCTCACAGACGGTATTATCAAACGTTTATAATTCATGGGGGCCTGCAAACGGATATGACCATTATCGTTCTATGCAGTCAATAACAGCTTGGATTACACAAACTGAAGCTGACAAAATATCAGGAGTATCAACTACTTATGACTTAATAACACAAAACCCTCATAAGGATGTTACGCTAAACGCTCTAAATGTCTATGCAGTTTGTGTAGAATAA
>TAU
AAGTAAGGCCAGCATTACTGAGATTAAGGCTGATAAAACAACTGCGAAGGCAAATGGTTCCGATGCTGTTACATACACTGTTAAAGTGATGAAGGGGGAGCAACCAGTTCAGAATCAATCTGTTACATTCTCAACAAACTTGGGGTTGTTAAACGGTAAGCCTCAAACGCAAACAACAACCACGGGAAGTGATGGTCGTGCGACGATAACGCTGACTTCCAATTCTGCTGGTAAAGCGACTGTGAGTGCGACTGTCAGTGGTGGGGGTGATGTTAAAGCACCTGAGGTTACATTTTTTGATGAACTGAAAATTGACAACAAGGTTGATATTATTGGTAACAACGTTACTGGTGACTTACCTAATATTTGGTTGCAATATGGTCAGTTTAAACTGAAGGTAAGCGGTGGTAACGGTACATATTCATGGCATTCAGAGAATACCAATATTGCGACTGTTGATGAATCAGGGAAAGTAACCTTGAATGGTAAAGGTACGGCAGTAATTAATGTTACATCTGGTGATAAGCAAACAGTAAGCTACACTATTAAAATGCCGACTTATATGATAAGAGTAGATAAGAAAGCTAATTATACTGGTGCCATGACCATTTGTAATAATTCTTTACCATCCTCACAGAAGGTATTATCAGACATTTTTAATTCTTGGGGGCCGGCAAATAAATATGGTAACTATGGTTCTATGAGCTCAATAACTGCCTGGATTAAACAAACTGAAAGTGATAAAAAATCTGGGGTGTCAACTACTTATGACTTAATAAGGCAAAACCCTCATAATAGTGTTACGCTAAACGCTCCAAATGTCTATGCAGTTTGTGTAGAATAA
>alpha_2
AGTAAGGCTAGTATTACCGAGATTAAGGCTGATAAAACAACTGCAGTAGCAAATGGTAAGGATGCTATTACATACACTGTTAAAGTGATGAAAAATGGTCAGCCGGTATCCGGCGAAGAAGTGACCTTTACGAAGACCTTAGGTACATTAAGTAAGCCCACTGAAAAAACGGATGCGAATGGCTATGCTAAAGTAACATTAACATCGGCGACTCAAGGGAAATCCCTCGTTAGTGCCCGTGTTAGCGATGTCGCCGTTGACGTCAAAGCCCCTGAAGTTGAATTTTTTACACCACTTACAGTTGATGACGGTAATATTGAAATTGTTGGAACCGGAGTTAAAGGGAAGTTACCTACTGTATGGTTGCAATATGGACAAGTTAAGCTGAAAGCTAGCGGAGGTAACGGAAAATATACATGGCGCTCAGCAAATACAGCGATTGCTTCGGTGGATGCCTCTTCTGGTCAGGTCACCTTAAAAGATAAGGGAACTACAACTATTACCGTTGTCTCAAGTGATAATCAAACCGCAACTTATACTATTGCAACACCTAATAGCCTGATTGTTCCTAATATGAGCAAGCGTGTGACTTATAATGATGCTGTGAATACATGTAAGATTTTGGGAGGAAGATTACCCTCTTCTCAGGATGAACTGAAAAATGTCTTTAACACATGGGGAGCAGCAAATAAATATGAATATTACAAGTCTCTTCTGACTATAAGCTCTTGGGTACAACAAACAACAGTGGATATGAACAGTGGAGTTGCAAGTACATACGATCTAGTTAGACAAAACCCTCTTAATGGCGTTGAAACTAATACTACTAATGCCTATGCCACTTGCGTAAAATAA
>alpha_1
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAAACAACGGCAGTAGCAAATGGTCAGGATGCTATTACATACACTGTTAAAGTGATGAAGGGGGATAAGCCTGTATCTAATCAGGAAGTGACCTTTACGACGACCTTAGGTAAGTTAAGTAATTCCACTGAAAAAACGGATACGAATGGCTATGCCAAAGTAACATTAACATCGACAACTCCAGGAAAATCACTCGTTAGTGCCCGTGTTAGCGATGTCGCCGTTGATGTCAAAGCACCTGAAGTTGAATTTTTTACAACGCTTACAATTGATGACGGTAATATTGAAATTGTTGGAACCGGAGTTAAAGGGAAGTTACCCACTGTATGGTTGCAATATGGTCAAGTTAATCTGAAAGCCAGCGGAGGTAACGGAAAATATACATGGCGCTCAGCAAATCCAGCAATTGCTTCGGTGGATGCTTCTTCTGGTCAGGTCACCTTAAAAGAGAAGGGAACTACAACTATTTCCGTTATCTCAAGTGATAATCAAACTGCAACTTATACTATTGCAACACCTAATAGTCTGATTGTTCCTAATATGAGCAAGCGTGTGACCTATAATGATGCTGTGAATACATGTAAGAATTTTGGAGGAAAGTTGCCGTCTTCTCAGAATGAACTGGAAAATGTCTTTAAAGCATGGGGGGCTGCAAATAAATATGAATATTATAAGTCTAGTCAGACTATAATTTCATGGGTACAACAAACAGCTCAAGATGTGAAGAGTGGTGTTGCAAGTACATACGATTTAGTTAAACAAAACCCTCTGAATAATATTAAGGCTAGTGAATCTAATGCTTATGCCACTTGTGTAAAATAA
>beta
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAAACAACTGCAGTAGCAAATGGTCAGGATGCTATTACATACACTGTTAAAGTGATGAAAAATGGTCAGCCATTATCCGGTGAAGAAGTGACTTTCTTTACGGATTTTGGGGCATTGGATAAAACTAAAGTAACGACCGATCAGAGTGGTTATGCTACTGTAAAACTGTCATCCAGTACTTCAGGCAAAGCTATTGTTCGCGCAAAAGTGAGTGATGTTGATACAGAAGTTAAGGCTGCCGCTGTTGAGTTTTTTGCCTCCCTGAGTATTGATAGCAATAAAGTTACTGTAATTGGTACTGGTGTCACTGGTTCTCTGCCAAAGAACTGGTTGCAGTATGGTCAGGCTAAGTTACAGGCAACGGGAGGTAATGGGAAATATACATGGAAATCCAGTGATACTAAAATTGCTTCTGTTGATAGTTCGGGAGTGATAACCCTGAATGGGAAAGGGAGTACCATAATTACGGTCGTATCTGGTGATAATCAGAGTGCGACATATACAATTAGTACTCCTGACAAGATTGTAATAGCTGTGGATAAAATTAATCGGATGACATATTCTGCTGCCGAAAACAAGTGCCGAACAATTAGTGCAAATTTAGCACCGTCAAAGAGTATATTGGCAGACACATATTCTAAATGGGGTGCAGCGAATAAATATTCTTACTATTCTGGTTCTAATTCATTGACAGCTTGGATTACACAATCCTCTTCTGAGTTACCATCGGGTGTATCAAGTACATATGATTTAATCACAACAAACTCTTTGACAAAAGTTAAGACAACAGATAACAATGCCTTTGCAGTCTGTGTGAAATAA
>RHO
ACCAAGGCCAGTATTACTGAAATTACCAATGATAAGTCAACGGCAATAGCAAATGATAAGGACGCCATTACATATACTGTTAAGGTAATGAAAAATGATCAGCCAGTACCGAATCATTTAGTCACATTCACGACGACTTTTGGTAAGTTCAATGGTAAACAGAGCTCAGAAACCGTAACAACAGGGAACGACGGACGCGCTATAGTGACACTGACATCAGGTATAGCCGGAAAAGCAATTGTCAGTGCAAAAGTTAATGAAGTAAATACCGAAGTACAAGCAAAGACTGTTGAGTTTTTCTCAAAACTAGCAGTCGATAGCAATATAGATATTATTGGCACATCAGTGACAGGAGCTTTGCCTGATGTCTGGTTGCAATATGGCCAGTTTAAATTGACTGCTAAAGGTGGAAATGGGAAATATACATGGTATTCTGAAAATCCAAATATTGCATCCGTGGAGTCAACTACAGGACTAGTAACCTTAAATAATAAAGGGAGTGTTAAAATATACGTTTTATCAGATGATAAGCAAACAGCAAGTTATACAATATCTACACCCAAAGATATGATAACAGCTGATATAACAAAACCAGCTATCTATAGTATTGCATTCAAGCAGTGCCAAACAAAAAATGGAAGACTACCATCTTCACAGCGTGAATTAGAAAACGTCTTTAATCTCTGGGGGGCTGCAAATAACTATAGATATTATTCTGCCAAACAATCTATAACAGCATGGATTCAACAATCTGGTTCAGATTCTAGCGGTGGAGTGACTACAACATATGATCTGATAAAGAAGAATCCCCGACAGAATGTTGATGTGAAAAGTCAAAATGTCTATTCGGTTTGCGTTGAATAA
>Epsilon
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAAACAACAGCGAAGGCAAATGGTTCTGATGCGATTACCTATATTGTTAAAGTAATGAAGAATAACCAACCAGAAGCAAACCATTCTGTTACATTCTCAACGAACTTTGGTAATCTGGGGGGAAATTCTAATACCCAAATTGTGAAAACGGATAAAGATGGTAGGGCTACGGTAAAACTGACATCTGGCGTTGCAGGTAATGCTGTTGTTAGTGCAAAAGTCAGCGAAGTTAATACAGAGGTTAAGGCTCCTGAGGTAAAATTCTTCTCAGTTCTGAGCATTGATAGTAATGTGAGTATTATTGGAACCTCCGCTAATGGCGCTTTACCTAATATTTGGTTGCGATATGGTCAGTTTAAGCTGACAGCCAAAGGTGGCGATGGGAAATATCAATGGCGCTCTCAAGATCCAAGTGTTGCATCAGTTGATGCTTTAACTGGTCGAGTTACTTTGCTGAAGAAAGGAACAACAACAATTGAAGTTGTATCGGGTGATAACCAAACAGCAATGTATACAATTAATACACCTACAAAATTTATATCTGTGGAGACACAAAATAAAGTAGTCTATAGTGATGCTGAGGCAACATGTAGAATGAATAATGCACGCTTGCCGTCATCTACGAGTGAGCTAAAGGATGTGTATAATAAATGGGGCGCCGCCAATAGTTATGAAGGCTATAAAGGTAAAAAAACAATAACAGCATGGACACAGCAAACTGAGGATGATAAACAAAAAGGTTGGACTAGTACATTTGACATAGTTACAAAAAATGAAATCCCTAGTAATGGCAGTAATAGTAAAGTCCACGTGAATAAAGCTAACGCTTTTGCCGTCTGTGTAAGATGA
>Kappa
ACCAAGGCCAGCATTACTGAGATTAAGGCTGATAAAACAACTGCAGTAGCAAATGGTCAGGATGCTATTACATACACTGTTAAAGTGATGAAAAATGGTCAGCCATTATCCGGTGAAGAAGTGACTTTCTTTACGGATTTTGGGGCATTGGATAAAACTAAAGTAACGACCGATCAGAGTGGTTATGCTACTGTAAAACTGTCATCCAGTACTTCAGGCAAAGCTATTGTTCGCGCAAAAGTGAGTGATGTTGATACAGAAGTTAAGGCTGCCGCTGTTGAGTTTTTTGCCTCCCTGAGTATTGATAGCAATAAAGTTACTGTAATTGGTACTGGTGTCACTGGTTCTCTGCCAAAGAACTGGTTGCAGTATGGTCAGGCTAAGTTACAGGCAACGGGAGGTAATGGGAAATATACATGGAAATCCAGTGATACTAAAATTGCTTCTGTTGATAGTTCGGGAGTGATAACCCTGAATGGGAAAGGGAGTACCATAATTACGGTCGTATCTGGTGATAATCAGAGTGCGACATATACAATTAGTACTCCTGACAAGATTGTAATAGCTGTGGATAAAATTAATCGGATGACATATTCTGCTGCCGAAAACAAGTGCCGAACAATTAGTGCAAATTTAGCACCGTCAAAGAGTATATTGGCAGACACATATTCTAAATGGGGAGCAGCGAATAAATATTCTTACTATTCTGGTTCTAATTCATTGACAGCTTGGATTACACAATCCTCTTCTGAGTTACCATCGGGTGTATCAAGTACATATGATTTAATCACAACAAACTCTTTGACAAATGTTAAGACAACAGATAACAATGCCTTTGCAGTCTGTGTGAAATAA
>lambda
ACTAAGGCCAGCATTACTGAGATTAAGGCTGATAAGACAACTGCAGTAGCAAATGGTAAGGATCCTATTAAATATACTGTAAAAGTTATGAAAAACGGTCAGCCTGTGAAGGGCTATGATGTAACGTTTTTAACAACTGCAGGTAATCTGAGTAAAACAAAAGAGTTAACAGATAAAGATGGATATGCAACTGTAAACTTGACCTCTAATGCTGCGGGAAAAGCTGTTGTTAGTGCCAAGGTTAGCGATGTTAATTCTGAGGTTAAAGCTTCCGAAGTTGAGTTCTTCACAGAACTGAGTATCAATAAAAATGTAGAGGTTCTTGGTACAAAGGCTAGTGGTGAGTTACCTGATGTCTGGCTACAATATGGTCAGATTAAATTAAATGTTAATGGTGGCAATGATAAATACTCATGGAGCTCCAGTAATCCTAATATTGCGTCGATAGATGCCTCTTCCGGTATTATTACCTTAAAGGAAAAGGGGGAGGCTGTTATTAAGGTTGTATCTGGGGATAAGCAAACGGCGACATACACTATTTCTACACCAAAGAAAATAGTGAGTGTTAATTCGGGTTCTGGGGTTAATTATAATTCGGCCAGTTCTATTTGCGGAAAAATTAATGGCTCGTTACCATCTTCAATTGCAGAGCTGGAAACACTGTATAATAAATGGGGAGCAGCCAATAACTATCAACACTATACTCAATCATCTATAACAGCATGGACGTTGCAAACCAGTGATGATGTGAAGCAAGGGGTCACCAGTACTTACGATTTGGTAAGGAAAAACCCACAATTAAATAAAGTTAATATAAACGACAATAATGCTTATGCTGTCTGTGTGAAATAA
>zeta
ACTCAGGCCAGTATTACTGAGATTAATGCTGATAAAAAAACAGCGAAGGCAAATGGATCTGATGCGATTACCTATACTGTTAAAGTGATGAAAGACGGTAAGCCGTTATCTGCCCAGGATGTGACCTTTACGGCAACCTTAGGTACGTTAAGTAAGTCCACTGAAAAAACGGATGCGAATGGCTATGCTAAAGTAACATTAACATCGAAGACTACAGGAAAATCGCTTGTTAGCGCTCGTATTAGCGGGAGCGCAATTGATGTCAAAGCACCTGAAGTTGAATTCTTCACCCCACTTGCTATTGATGACGGTAATGTTGAAATTGTAGGAACTGGAATCAAAGGGACGTTACCTACGGTATGGTTGCAATATGGACAAGTTAAGTTGAAAGCCAGTGGAGGTGATGGGAAATATACATGGAGTTCGGCAAATACTGGAATTGCTTCGGTGGATTCCACTGGTCAAGTTACCTTAAGGGATAAGGGGAGTACAACTATTACCGTTGTATCAGGTGATAAGCAAACTGCAACTTATATAATTGCAAGACCTAGCAGTATGATTGTCTCTATAAACGAGCGTATGACTTATAATAACGCAATGAGTTCATGCCAGAGCTTGAGCGGGAGGTTGCCATCTTCTCAAAAAGAACTAGCTGATGTATTCGATACATGGGGGGCAGCAAATAAATATGAACATTATGAGACTCGTAACGCCATGATATCGTGGATAAAACAAACAGATCAGGATATGAGGCAGGGAGTTGCAAGTACATACGATCTGGTTAAGAAAAATCCGCTTACTAATAAAGTAGATATAAATAAGCCTAATGCTTATGCAACTTGCGTAAAATAA
>eta
GATAAGACAACAGCGAAGGCAGATGGTTCTGATGCGATTACCTATACTGTTAAAGTAATGAAGAATAACCAACCAGAAGCAAACCATTCTGTTACATTCTCAACGAACTTTGGTGATCTGGGGGGAAATTCTAATACCCAAATTGTGAAAACGGATAAAGATGGTAGGGCTACGGTAAAACTGACATCTGGCGTTGCAGGTAATGCTATTGTTAGTGCAAAAGTCAGCGAAGTTAATACAGAGGTTAAGGCTCCTGAGGCAAAATTCTTCTCAGTTCTGAGCATTGATAATAATGTGAATATTATTGGAACCTCCGCTAATGGTGCTTTACCTAATATTTGGTTGCGATATGGTCAGTTTAAGCTGACAGCCAAAGGTGGCGATGGGAAATATCAATGGCGCTCTCAAGATCCAAGTGTTGCATCAGTTGATGCTTTAACTGGTCGAGTTACTTTGCTGAAGAAAGGAACAACAAAAATTGAAGTTGTGTCTGGCGATAAGCAAACCGCGACATATACAATTAATACACCTGAGAAAATTATAACTGTAGAAACACAAAATAAAGTAATCTATAGCGTTGCTGAAGCAACATGCAGCACGAATAGTGGTCGCTTACCGTCATCTACGAGTGAGCTAAAGGATGTGTATAATAAGTGGGGCGCTGCCAATAGTTATGAAGGTTATAAAGGTAAAAATACAATAACAGCATGGACTCAGCAAACTGCGGCTGATAAACAATCTGGTTGGACTAGTACATTTGACATAGTTACAAAAAATGAAATTCCTAATAATGGAAGTAATAGTAAGGTCAATGTGAACACAGCTAACGCTTTTGCCGTCTGTGTAAGATGA
>xi
TGCGATTACCTATACTGTTAAAGTAATGAAGAATAACCAACCAGAAGCAAACCATTCTGTTACATTCTCAACGAACTTTGGTAATCTGGGGGGGAATTCTAATACCCAAATTGTGAAAACGGATAAAGATGGTAGGGCTACGGTAAAACTGACATCTGGCGTTGCAGGTAATGCTATTGTTAGTGCAAAAGTCAGCGAAGTTAATACAGAGGTTAAGGCTCCTGAGGTAAAATTCTTCTCAGTTCTGAGCATTGATAGTAATGTTAATATTATTGGAACCTCTGCTACTGGCGCTTTGCCTAATATTTGGTTGCAATATGGTCAGTTCAAATTGACTGCTAAAGGTGGTGATGGGAAATATCAATGGCGCTCTCAAGATCCAAAAGTTGCATCAGTTGATGCTTTAACTGGTCGAGTTACTTTGTTGAAGAAAGGAACAACAACAATTGAAGTTGTGTCGGGTGATAACCAAACTGCAACGTATACAATTAATACACCTATAAAAATTATATCTGTGGAGACAAAAAATAAAGTAGTCTATAACGATGCTGAAGCAATATGTAGAACGAATAATGGCCGTTTACCGCTATCTACGAATGAGTTAAAGGACGTGTATAATAAATGGGGAGCGGCCAATAGTTATGAAGGCTATAAAGGTAAAAACACAATAACAGCATGGACTCAGCAAACAGAGGATGATAAAACTAAAGGTTGGACTAGTACATTTGACATAGTTACTAAAAATGAAATCCCTAGTAATGGAAGTAATAATAAGGTCAATGTGACAGCAGCTAATGCCTTTGCTGTCTGTGTAAGATGA
>nu
GCAGTAGCAAATGGTAAGGATGCTATTACATATACTGTTAAAGTGATGAAAAATGGTCAGCCGTTATCCGGTGAAGAAGTGACCTTTACGAAGACCTTAGGTACGTTAAGTAAGTCCATTGAAAAAACGGATGCGAATGGCTATGCTAAAGTAACATTAACATCTGTGACTCAAGGGAAATCCCTCGTTAGTGCCCGTGTTAGCGATGTTGCCGTTGACGTCAAAGCACCTGAAGTTGAATTTTTTACACAGCTTACAATTGATAACAGTAATATCGAAATTGTTGGTACAGGAGTTAAAGGGACGCTGCCTACTGTATGGTTGCAATATGGTCAAGTTAAGTTGAAAGCTAGCGGAGGTAACGGAAGATATACATGGCGTTCAGCAAATACAGCAATTGCTTCGGTGGATGCTTCTTCAGGTCAAATCACCTTAAAAGATAAGGGGACTACAACTATTACCGTTGTATCAGGTGATAATCAAACCGCAACTTATACTATTGCAACACCTAGCAGTATGATTGTCTCTAATATAAACAGTCGTATGACTTATAGTGAGGCAATGAGTTCATGCCAGAGCTTGAGAGGAAGGTTACCATCTTCTCAAAGTGAACTGGCTGATGTATTCGATACATGGGGGGCGGCAAATAAATATGAATATTATAAGACTAGCACAAGCATAATATCGTGGATAAAACAAACCGATGATGATGTGAAGAAGGGCGTTGCAAGCACATACGATCTGGTTAAGAGAAATCCTCTTACTAGTTCAGTAAAAACAGACGCTCGTAATGCTTATGCAACTTGCGTAAAATAA
>omicron
ACTGCAGTAGCAAATGGTAAGGATGCTGTCACATACACTGTTAAAGTAATGAAAGATGGTAAGCCGTTATCCGGTGAAGAAGTGACCTTTACGACGACCTTAGGGACGTTAAGTAAGTCCACTGAAAAAACGAATACGAATGGCTATGCTAAAGTAACGTTAACATCGACGAATCAAGGAAAATCACTCGTTAGTGCCAGCGTTAGCAATGCCGCAGTTGATGTTAAAGCTCCTGAAGTTGAATTTTTTACACAGCTTACAATTGATAACGGTAATGTTGAAATTGTTGGAACCGGAGCTAAGGGGAAACTACCTAATGTATGGTTGCAATATGGCCAAGTTAATCTGAAAGCTAGCGGAGGTAACGGGAAATATACATGGTACTCAGCAAATCCGGCAATTGCTTCGGTGGATCCCTCTTCTGGTCAAGTCACATTAAAAGATAAGGGGGAAACAACAATTACCGTTGTATCTGGGGATAAACAAACGGCCATCTATACAATTGCTATGCCTAATAGCATAGTGAGCGTTAACTCAAGTGGTAGGGTTGATTATAATACAGCCAATAACATTTGCAAAAATATTAAAGGTTCGTTACCTTCTTCAATTAAAGAGCTGAAGGATTTATATGATGATTGGGGGGCTGCAAATAAATATCAACATTATTCCCAGGAATCTATAACTGCCTGGACATTACAAACTAGTGAGAATAAAGTGCAAGGGGTTGCAAGTACTTACGATTTGGTAAGGAAAAATCCACTGATAGATAAAGTTGATATAGCAGGTAATTATGCCTATGCAGTCTGTGTGAGATAA
>pi
ATGCGATTACCTATACTGTTAAAGTAATGAAGAATAACCAACCCGAAGCAAACCATTCTGTTACATTCTCAACGAACTTTGGTAAGTTGGATGGTAATTCTAATACCCAAACTGTGAAAACAGATGAAAATGGTAAAGCTACGGTAAAACTGACATCTGGCGCTGAGGGGAGTGCTGTCGTTAGTGCAAAAGTCAGCGAGATTAATACAGAAGTTAAGGCTCCTGAGGTAAAATTCTTCTCAGTTCTGAGTATTGATAATAATGTGAATATTATTGGAACCTCCGCTAATGCTGCTTTACCTAATATTTGGTTGCAATATGGTCAGTTTAAGCTGACAGCCAAAGGCGGTGATGGAAAATATAAATGGCATTCTAAAGATACAAGTGTTGCATCAGTTGATGCTTCTACTGGTCAGGTTACTTTGCTGAAGAAAGGAACGACAACAATTGAAGTAGTGTCTGGCGATAACCAAACCGCGGCATATACAATTAATACACCTGATAAGATAATATCTGTAGAATCACAAAATACAGAAATCTATAGCTCTGCTGAAAAAACGTGTAGTATTAATAGTGGTCGTTTACCGTTATCTACGAATGAGTTAAAGGATGTGTATAGTAAATGGGGAGCAGCCAATAGTTATGAGAGGTATAAAGATAAAAATACAATAGCAGCATGGACACAGCAAACTGAGGAAGATAAAAAAGCAGGTTGGACTAGTACATTTGACATAGTTACGCAAAATGAAATTCCTAGTAATGGAAGTAATAGCAAGGTCGACGTGAAAAAAGCTAACGGTTTTGCCGTCTGTGTAAGATGA
>UPSILON
CTGATAAGACAACTGCAGTCGCAAATGGTAATGATGCTGTTACATACACTGTTAAAGTGATGAAAAATGGGCTGCCGGAAAAGGGACACACAGTTACATTCTCAACTGATTTAGGTAAATTGAATCAACAGACTGTAGCAACGGATCAAGATGGTTTTGCTTCGGTAACATTGACTTCTGATTCTGTCGGTAAAGCAGTTGTTAGTGCAAAAGTCAGTGAAGCAGGCTCTGTGATTAATGCAGATGTAGTTAACTTCTTCGCTACGCTTAGTATTGATAATAATGTCGAGATCGTTGGTACAAAAGTTCGCGGTGAATTACCTGATATCTGGTTACAATATGGGCAGGTTAAGTTGATTGCAAATGGCGGTAATGCTGGATATTCATGGAGTTCCGATAATCCAGATATTGCATCGATAGATGCTAATACAGGTATTATCACCTTAAATAAAAAAGGAACTACAGTTATTAAAGTCATATCTGGTGATAAACAGATGGCAACGTATACGATTAAGCCACCTAAAGAAATAGTATCTTTAGATAAGAATGTAAAAGTTACATATGACGAAGCCAGTCGTATATGTTTAACCAATGCTGCACATTTATCTGTCTCAGCAGACAGTTTGAAAAAAGTATATGGTCAATGGGGACCAGCAAATAAATATGGACATTATACCCAAGGAACCATAAATGGTTGGATACAACAAACGGATAAGGATAAAAAAGACGGCGTCGCAGCGACATATGATATCGTCACTAATAATACGGTGCCAAATGTTGATTCAACAAAGGTTCACGGTTATGCAGTTTGCATAAAATAA
```
## 1.2 Use Abricate on Pre-classified EHEC & EPEC Assemblies to Calculate the counts and percentages
✅ Overview
You already have:
A table of assemblies classified as EHEC or EPEC based on stx1/stx2/eae presence (`abricate_presence_absence_filtered(ide90cov90).tsv`).
A FASTA reference (intimin_subtypes.fasta) containing eae subtype sequences (γ, α1, α2, β, ε).

We’ll now:
- Separate assemblies into EHEC and EPEC folders.
- Run Abricate using intimin_subtypes database.
- Summarize counts and percentages of each eae subtype within each group.

Step 1. Create subtype database for Abricate

```
mkdir eae_subtype
cd eae_subtype
nano intimin_subtypes.fasta

mkdir -p /home/jing/miniconda3/envs/abricate/db/intimin_subtypes
cp intimin_subtypes.fasta /home/jing/miniconda3/envs/abricate/db/intimin_subtypes/sequences
abricate --setupdb
abricate --list   # should now show intimin_subtypes
```
✅ This makes a custom Abricate database called intimin_subtypes.


Step 2. Split assemblies by classification
Extract the list of EHEC and EPEC assembly IDs from previous analysis.
Save as `classify_ehec_epec.py` in the same folder as your TSV and run `python classify_ehec_epec.py`.
```
#!/usr/bin/env python3
import re
import pandas as pd

IN = "abricate_presence_absence_filtered.tsv"

# 1) Load
df = pd.read_csv(IN, sep="\t")

# 2) Extract a clean sample ID from the Sample path
#    Prefer the first GCA/GCF accession; fallback to basename without .tsv
def extract_id(s):
    if pd.isna(s):
        return None
    m = re.search(r"(GC[AF]_\d+\.\d+)", str(s))
    if m:
        return m.group(1)
    # fallback: basename without extension
    base = str(s).split("/")[-1]
    return re.sub(r"\.tsv$", "", base)

if "Sample" not in df.columns:
    # Some earlier runs used 'sample' lowercase
    sample_col = [c for c in df.columns if c.lower()=="sample"][0]
    df.rename(columns={sample_col: "Sample"}, inplace=True)

df["sample"] = df["Sample"].apply(extract_id)

# 3) Ensure gene columns are numeric 0/1
for col in ["stx1","stx2","eae"]:
    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

# 4) Classify
df["pathovar"] = "Other"
df.loc[(df["eae"]==1) & ((df["stx1"]==1) | (df["stx2"]==1)), "pathovar"] = "EHEC"
df.loc[(df["eae"]==1) & ((df["stx1"]==0) & (df["stx2"]==0)), "pathovar"] = "EPEC"

# 5) Write lists
ehec_ids = df.loc[df["pathovar"]=="EHEC", "sample"].dropna().drop_duplicates()
epec_ids = df.loc[df["pathovar"]=="EPEC", "sample"].dropna().drop_duplicates()

ehec_ids.to_csv("EHEC_list.txt", index=False, header=False)
epec_ids.to_csv("EPEC_list.txt", index=False, header=False)

# 6) Save a small mapping table (useful later)
df[["sample","pathovar","stx1","stx2","eae"]].to_csv(
    "ehec_epec_classification.tsv", sep="\t", index=False
)

print(f"EHEC: {len(ehec_ids)} isolates → EHEC_list.txt")
print(f"EPEC: {len(epec_ids)} isolates → EPEC_list.txt")
print("Saved: ehec_epec_classification.tsv")
```

Step 3 — Separate assemblies by pathovar and run Abricate
I already have assembly folders,`/home/jing/E.coli/ecoli_all/` containing .fna files.
I’ll copy files listed in the pathovar lists:
```
mkdir -p EHEC_assemblies EPEC_assemblies

while read id; do
  find /home/jing/E.coli/ecoli_all -name "${id}*.fna" -exec cp {} EHEC_assemblies/ \;
done < EHEC_list.txt

while read id; do
  find /home/jing/E.coli/ecoli_all -name "${id}*.fna" -exec cp {} EPEC_assemblies/ \;
done < EPEC_list.txt
```
Step 3. Run Abricate on each set
🧬 EHEC
```
mkdir -p abricate_intimin_EHEC
for f in EHEC_assemblies/*.fna; do
  base=$(basename "$f" .fna)
  abricate --db intimin_subtypes --minid 85 --mincov 90 "$f" > abricate_intimin_EHEC/${base}.tsv
done
abricate --summary abricate_intimin_EHEC/*.tsv > abricate_intimin_EHEC_summary.tsv
```
🧬 EPEC
```
mkdir -p abricate_intimin_EPEC
for f in EPEC_assemblies/*.fna; do
  base=$(basename "$f" .fna)
  abricate --db intimin_subtypes --minid 85 --mincov 90 "$f" > abricate_intimin_EPEC/${base}.tsv
done
abricate --summary abricate_intimin_EPEC/*.tsv > abricate_intimin_EPEC_summary.tsv
```

Step 4 — Summarize results (counts & percentages)
Create a new script named `summarize_intimin_subtypes.py`:
```
#!/usr/bin/env python3
import pandas as pd

def summarize_from_summary(file, label):
    # Load abricate summary
    df = pd.read_csv(file, sep="\t")
    
    # Drop the first two metadata columns (#FILE, NUM_FOUND)
    data = df.drop(columns=[df.columns[0], df.columns[1]])

    # Convert "." to NaN and non-numeric values to NaN
    data = data.replace(".", pd.NA)

    # Count non-empty cells per column (i.e. isolates with hits)
    counts = data.notna().sum().reset_index()
    counts.columns = ["Subtype_ID", f"{label}_Count"]

    # Calculate percentage
    total = len(df)
    counts[f"{label}_Percent(%)"] = (counts[f"{label}_Count"] / total * 100).round(2)
    counts["Pathovar"] = label

    return counts

# Summarize EHEC and EPEC datasets
ehec = summarize_from_summary("abricate_intimin_EHEC_summary.tsv", "EHEC")
epec = summarize_from_summary("abricate_intimin_EPEC_summary.tsv", "EPEC")

# Merge both summaries
combined = pd.merge(ehec, epec, on="Subtype_ID", how="outer").fillna(0)
combined = combined[["Subtype_ID", "EHEC_Count", "EHEC_Percent(%)", "EPEC_Count", "EPEC_Percent(%)"]]

# Map reference IDs to meaningful subtype names
subtype_map = {
    "eae_gamma": "Gamma intimin",
    "eae_alpha2": "Alpha 2 intimin",
    "eae_alpha1": "Alpha 1 intimin",
    "eae_beta": "Beta intimin",
    "eae_epsilon": "Epsilon intimin",
    "eae_sigma": "Sigma intimin",
    "eae_mu": "Mu intimin",
    "eae_theta": "Theta intimin",
    "eae_tau": "TAU intimin",
    "eae_RHO": "Rho intimin",
}

combined["Subtype_Name"] = combined["Subtype_ID"].map(subtype_map)
combined = combined[["Subtype_Name", "Subtype_ID", "EHEC_Count", "EHEC_Percent(%)", "EPEC_Count", "EPEC_Percent(%)"]]

# Save output
combined.to_csv("intimin_subtype_distribution_combined.csv", index=False)
print("\n✅ Intimin subtype summary created successfully!\n")
print(combined)

```
Run it:
```
python summarize_intimin_subtypes.py
```
<img width="929" height="127" alt="image" src="https://github.com/user-attachments/assets/0cea76fa-6285-47e0-b0de-8e6cb024cd6c" />

# Serotyping EHEC (1114 isolates) and EPEC (339 isolates) strains
### Option 1: ECTyper (an easy typer) 
(https://github.com/phac-nml/ecoli_serotyping) (https://pmc.ncbi.nlm.nih.gov/articles/PMC8767331/#s5)
1. Install ECTyper
```
conda install -c bioconda ectyper
```

2. Run ECTyper
```
ectyper -i EHEC_assemblies/*.fna -o EHEC_Serotype
ectyper -i EPEC_assemblies/*.fna -o EPEC_Serotype
```

3. Results Summary
serotype_distribution_summary.py
```


