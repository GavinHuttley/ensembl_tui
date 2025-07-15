# Download data from Ensembl

This command downloads the data from the Ensembl server saving it into `staging_path` as specified in the config file. This is the only step that requires Internet access.

In addition to writing the original files from Ensembl, an expanded version of the config file is also written into the download directory. This is required for the subsequent `install` step.

```
$ eti download -c example/sample.cfg
```

> **Note**
> Downloading takes advantage of multiple threads and can be interrupted and resumed.
