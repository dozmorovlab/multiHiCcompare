- During installation, several microarray annotation packages get installed. The most annoying is the 300 Mb KEGG data, "KEGGdzPathwaysGEO_1.18.0.tar.gz". Why are they installed? Can one get rid of them?

```
> devtools::install_github('jstansfield0/multiHiCcompare', build_vignettes = TRUE)
Downloading GitHub repo jstansfield0/multiHiCcompare@master
Installing 3 packages: hgu133a.db, hgu133plus2.db, KEGGdzPathwaysGEO
Installing packages into â€?/Users/mdozmorov/Library/R/3.5/libraryâ€™
(as â€?libâ€™ is unspecified)
installing the source packages â€?hgu133a.dbâ€™, â€?hgu133plus2.dbâ€™, â€?KEGGdzPathwaysGEOâ€™

trying URL 'https://bioconductor.org/packages/3.7/data/annotation/src/contrib/hgu133a.db_3.2.3.tar.gz'
Content type 'application/x-gzip' length 903473 bytes (882 KB)
==================================================
downloaded 882 KB

trying URL 'https://bioconductor.org/packages/3.7/data/annotation/src/contrib/hgu133plus2.db_3.2.3.tar.gz'
Content type 'application/x-gzip' length 2139642 bytes (2.0 MB)
==================================================
downloaded 2.0 MB

trying URL 'https://bioconductor.org/packages/3.7/data/experiment/src/contrib/KEGGdzPathwaysGEO_1.18.0.tar.gz'
Content type 'application/x-gzip' length 314643015 bytes (300.1 MB)
==================================================
downloaded 300.1 MB
```

- When `multiHiCcompare` is loaded, `KEGG.db` is also loaded. Why? It is not a part of package's functionality.

```
> library(multiHiCcompare)

KEGG.db contains mappings based on older data because the original resource was removed from the the public domain before the most recent update was
  produced. This package should now be considered deprecated and future versions of Bioconductor may not have it available.  Users who want more current
  data are encouraged to look at the KEGGREST or reactome.db packages
```

- Add description/functionality how to use vignettes. Currently, they are not working.

```
> browseVignettes(package = "multiHiCcompare")
No vignettes found by browseVignettes(package = "multiHiCcompare")
```

# `juiceboxVisualization.Rmd`

- `exportJuicebox(rao2017, logfc_cutoff = 2, logcpm_cutoff = 1, p.adj_cutoff = 0.001, file_name = "rao2017Annotations.txt")` - Add "rao2017Annotations.txt" file as the example data to the package, so the user won't need to search for it.