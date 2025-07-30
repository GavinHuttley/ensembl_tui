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
4. Defining intergenic as last gene stop and current gene start
5. Write these out to a tab delimited file with the correwct column headings

```python exec="on" result="ansi" workdir="./docs" source="above"
import cogent3

table = cogent3.load_table("human_data/homo_sapiens-114-gene_metadata.tsv")
table = table.sorted(columns=["seqid", "start", "stop"])
# make sure the seqid column is a string type
table.columns["seqid"] = table.columns["seqid"].astype(str)
seqids = table.distinct_values("seqid")
# let's just do chrom 22 ... because that's all we have!
chrom22 = table.filtered(lambda x: x == "22", columns="seqid")
start_stop = chrom22.to_list(["start", "stop"])
inter_genic = [("22", 0, start_stop[0][0], 1)]
last_end = start_stop[0][1]
for start, stop in start_stop:
   inter_genic.append(("22", last_end, start, 1))
   last_end = stop

intergen_tab = cogent3.make_table(header=["seqid", "start", "stop", "strand"], data=inter_genic)
intergen_tab.write("data/chrom22-intergenic.tsv")
```

```bash exec="1"
rm -rf human_data # markdown-exec: hide
```

> **Warning**
> The above does not handle the case where genes overlap!