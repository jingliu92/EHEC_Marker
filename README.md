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
python make_presence_matrix.py --glob "/home/jing/E.coli/blast_hits_only/*.tsv"
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
🌀 Step 3: ABRicate with a custom DB (threshold: --minid 80 --mincov 80)
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
🌀 Step 6 Check how many assemblies have hits? (2390)
```
wc -l abricate_presence_absence_filtered.tsv

# 2390
```
<img width="837" height="41" alt="image" src="https://github.com/user-attachments/assets/3b084a46-cd6c-44ad-972a-039052be807c" />
