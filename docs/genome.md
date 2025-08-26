# Querying Genomes

## Summary of this installation

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti installed -i data/apes-114
```
> **Note**
> :material-download: [Download all the data](ensembl_tui_data.zip) (zip, ~196 MB).

## Summary for a species {#summary-species}

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti species-summary -i data/apes-114 --species human
```

## Export gene meta-data for a species {#export-genes}

> **Note**
> The list of data from this query only covers human chromosome 22 because we are using a custom subset of the original Ensembl data.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti dump-genes -i data/apes-114 --species human -od human_data
```

```console exec="1" source="console" result="ansi" workdir="./docs"
$ head human_data/homo_sapiens-114-gene_metadata.tsv
```

## Defining intergenic regions

In order to utilize `ensembl-tui` for sampling non-genic regions you need to write code that will produce a coordinate file. We're going to do that here using `cogent3`. In brief, the algorithmic steps are

1. Load the metadata file into a `cogent3` table
2. Sort the table by the genomic coordinate columns seqid, start, stop
3. For each seqid (e.g. "22" for chromosome 22)
4. Get the start, stop coordinates for all genes and merge overlapping
5. Defining intergenic as last gene stop and current gene start
6. Write these out to a tab delimited file with the correct column headings

```python exec="on" result="ansi" workdir="./docs" source="above"
from cogent3 import load_table, make_table
from cogent3.util.misc import get_merged_overlapping_coords

table = load_table("human_data/homo_sapiens-114-gene_metadata.tsv")
# make sure the seqid column is a string type
table.columns["seqid"] = table.columns["seqid"].astype(str)
table = table.sorted(columns=["seqid", "start", "stop"])

# if we were doing this for real, we would work on each unique
# seqid at a time
seqids = table.distinct_values("seqid")

# but we just do one seqid here, chrom 22
seqid = "22"
chrom22 = table.filtered(lambda x: x == seqid, columns="seqid")
start_stop = chrom22.to_list(["start", "stop"])

# we use the utility function to merge overlapping gene coordinates
start_stop = get_merged_overlapping_coords(start_stop)

# define the remaining constants to be output into the tsv file
# required by eti
species = "homo_sapiens"
strand = 1
inter_genic = [(species, seqid, 0, start_stop[0][0], strand)]
last_end = start_stop[0][1]
for start, stop in start_stop:
   inter_genic.append((species, seqid, last_end, start, strand))
   last_end = stop

intergen_tab = make_table(header=["species", "seqid", "start", "stop", "strand"], data=inter_genic)
intergen_tab.write("data/chrom22-intergenic.tsv")
```

```bash exec="1"
rm -rf human_data # markdown-exec: hide
```