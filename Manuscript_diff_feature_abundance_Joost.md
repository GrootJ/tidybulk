tidybulk Manuscript code - differential feature abundance
================

add renv for tracking package management - done after version controlled
notebook initiated and run only execute once

note renv snapshot looks for package dependencies throughout the
project - I added several test tidybulk notebooks to the project (from
tidybulk paper etc) - they differ in packages called which is detected
by renv. Renv auto-detects all packages called across notebooks and I
added a few notebooks from other (tidy-omics) publications to test
tidy-bulk here. So unlikely all packages detected by renv are used
across these notebooks are true requirements for tidy-bulk to function.

after snapshot choose option 2: Install the packages, then snapshot.
this throws error “Error: package ‘nanny’ is not available” - which
looks to be stemangiola’s nanny package
<https://github.com/stemangiola/nanny> which is not longer maintained
and superseeded by tidybulk

note that renv sets the libpaths to project local check some renv
settings

install bunch of packages not captured by renv::snapshot() installation
(execute only once)

notes package ‘nanny’ is not available - but that is depricated
precursor to tidybulk - ignore here

also install required bioconductor packages (execute only once)

check bioconductor packages install

note that the code below is from notebook Manuscript code -
transcriptional signature identification
<http://127.0.0.1:16527/library/tidybulk/doc/manuscript_transcriptional_signatures.html>

Load the installed libraries

# Load data from pasilla dataset

``` r
pasCts = system.file("extdata",
                                         "pasilla_gene_counts.tsv",
                                         package = "pasilla",
                                         mustWork = TRUE)
pasAnno = system.file(
    "extdata",
    "pasilla_sample_annotation.csv",
    package = "pasilla",
    mustWork = TRUE
)
cts = as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))
coldata = read.csv(pasAnno, row.names = 1)
coldata = coldata[, c("condition", "type")]

# Create tidybulk object
counts =
    cts %>%
    as_tibble(rownames = "feature") %>%
    pivot_longer(names_to = "sample",
                             values_to = "count",
                             cols = -feature) %>%
    left_join(
        coldata %>%
        as_tibble(rownames = "sample") %>%
        mutate(sample = gsub("fb", "", sample))
    ) %>%
    mutate_if(is.character, as.factor)
```

alternative (NOT DONE) load dataset from tidybulk README.Rmd (does not
work?)

``` r
# loading this summarized experiment object doesn't work for me - needs to load from dev branch? only have main forked.. 
# counts_SE = here("dev/counts_SE.rda") |> load()
# counts_SE = here("counts_SE.rda") |> load()
# tibble_counts = counts_SE %>% tidybulk() %>% as_tibble()

# can load se_mini which is used in introduction.Rmd - proceed running functions with se_mini 
# data(se_mini)
# tibble_counts.se_mini = tidybulk::se_mini |> tidybulk() |> as_tibble()
```

# Create a tt object with unique raw and normalised counts

``` r
tt_scaled <- 
    tidybulk(counts, sample, feature, count) %>%
    aggregate_duplicates() %>%
    identify_abundant() %>%
    scale_abundance()

# Plot count densities
tt_scaled %>%
    pivot_longer(
        c(count, count_scaled),
        values_to = "count", 
        names_to = "Normalisation"
    ) %>%
    ggplot(aes(count + 1, group=sample, color=type)) +
    facet_grid(~Normalisation) +
    geom_density() +
    scale_x_log10()
```

![](Manuscript_diff_feature_abundance_Joost_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

paired end has lower density is logical (sequence from other end of
insert on flowcell clusters is more likely to yield less (complete)
sequences) depending on flowcell(s) design single end alignment (or SE
k-mer counting) could provide a technical solution for the single vs
paired end differences

``` r
# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- tt_scaled %>% reduce_dimensions(method="MDS", .dims = 3)

# Plot all-vs-all MDS dimensions 
tt_mds %>%
    pivot_sample() %>%
    GGally::ggpairs(columns = 7:9, ggplot2::aes(colour=condition))
```

![](Manuscript_diff_feature_abundance_Joost_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

# Adjust for visualisation

``` r
# Adjust for visualisation
tt_adj <- tt_mds %>% adjust_abundance(~ condition + type)

# Visualise the association between reduced dimensions and factors
tt_mds_adj_mds <- 
    tt_adj %>%
    filter( count_scaled_adjusted %>% is.na %>% `!`) %>%

    # Calculate reduced dimensions on the adjusted counts as well
    reduce_dimensions(
       .abundance = count_scaled_adjusted, 
       method="MDS", .dim = 3
    )
```

OLD NOTES from running first version of notebook (first time running
notebook) when earlier chunks to set up package environment were not
added yet:

NOTE tidybulk said it required combat in the adjust_abundance function -
then installed sva (package w Combat) and used it for adjustment for
unwanted variation this means function adjust_abundance executes batch
correction - in this case over paired-end vs single-end - which may not
be a very common use case in the documentation a good example of batch
correction (giving idea of when to use it) is probably important - batch
correction can distort true biological signal a first idea for an
example would be multi-site data (maybe TCGA lung squamous vs lung
adenocarcinoma? normal lung controls across batches should map on top of
each other? look further into this..) UPDATED: stemangiola now added
reference to combat in the adjust_abundance documentation

Data manipulation and visualisation

``` r
tt_mds_adj_mds %>%
    pivot_sample() %>%

    # First level reshaping
    pivot_longer(contains("Dim"), names_to = "Dim", values_to = ".value")   %>%
    separate(Dim, c("Dim", "Adj"), sep="\\.") %>%
    mutate(Adj = ifelse(Adj == "y", "non", "adj") %>% factor(c("scaled", "adj"))) %>%

    # Second level reshaping
    pivot_longer(c(type, condition), names_to = "covar", values_to = "which") %>%

    # Visualise the integrative plot
    ggplot(aes(y = .value, x = covar, fill = `which`)) +
    geom_boxplot() +
    facet_grid(Adj ~ Dim)
```

![](Manuscript_diff_feature_abundance_Joost_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

should the bottom “NA” be the adjusted data? there differences between
paired and single end seem to have been “corrected” and signal of
treated vs untreated on PC1 shows more clear separation..

``` r
tt_test <- tt_adj %>% test_differential_abundance(~ condition + type)
```

type (read type) is used as covariate - as test_differential_abundance
takes in raw count data (does its own correction, normalization)

Couple of notes: 1) Looking into test_differential_abundance function
documentation it is not clear whether count, count_scaled, or
count_scaled_adjusted columns are used as input defaults 2)
test_differential_abundance outputs the design columns used and the
design column tested - here: (Intercept) - this limited info on what
contrast the outputted fold change represents (and also a bit confusing
for a novice user)

look at subset of data relevant to differential abundance test

``` r
tt_test_diff_abundance_data <- tt_test %>% dplyr::select(feature, sample, count, count_scaled, count_scaled_adjusted, logFC, logCPM, F, PValue, FDR) %>% distinct()

tt_test_diff_abundance_data %>% head(n=20)
```

then look at fold change differential expression data only

``` r
tt_test_diff_abundance_data_FC_only <- tt_test %>% dplyr::select(feature, logFC, logCPM, F, PValue, FDR) %>% distinct()
tt_test_diff_abundance_data_FC_only  %>% head(n=20)
```

default diff testing method is edgeR - which outputs logCPM - a per gene
expression value it is not very clear what the outputted logCPM
represents here - it looks like it is the “Intercept” value of the
linear modeling

``` r
# MA plot
tt_test %>%
        keep_abundant() %>%
      pivot_transcript() %>%

    # Subset data
    mutate(significant = FDR<0.05 & abs(logFC) >=2) %>%
    mutate(feature = ifelse(significant, as.character(feature), NA)) %>%

    # Plot
    ggplot(aes(x = logCPM, y = logFC, label=feature)) +
    geom_point(aes(color = significant, size = significant, alpha=significant)) +
    geom_text_repel() +
    scale_color_manual(values=c("black", "#e11f28")) +
    scale_size_discrete(range = c(0, 2))
```

![](Manuscript_diff_feature_abundance_Joost_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

``` r
 tt_test %>%

     # Select top genes and reshape data
     inner_join( arrange((.), PValue) %>% distinct(feature) %>% head(6)) %>%

     # High level reshaping of the data.
     # All three count columns are shaped as two columns:
     # (i) the columns name and (ii) the value of those columns
     pivot_longer(
        c(count, count_scaled, count_scaled_adjusted),
        names_to = "Stage", values_to = "count"
     ) %>%

     # This allows the faceted plot
     ggplot(aes(x = Stage, y = count + 1, fill = condition)) +
        geom_boxplot() +
     facet_wrap(~feature) +
     scale_y_log10()
```

``` r
 # Heatmap
 tt_test %>%
        as_tibble() %>%

     # Select differentially abundant
     filter(FDR < 0.05 & abs(logFC) > 2) %>%

     # Plot
     heatmap( feature, sample, count_scaled_adjusted) %>%
     add_tile(condition) %>%
     add_tile(type)
```

heatmap does not seem to provide right visual for DE genes - perhaps a
scaling error? error reported = “tidyHeatmap says: (once per session)
from release 1.7.0 the scaling is set to”none” by default. Please use
scale = “row”, “column” or “both” to apply scaling”

snapshot libraries installed and display sessionInfo for
renv::snapshot() chose for option 1: Snapshot, just using the currently
installed packages. snapshot extensively looks for required packages -
not all are necessary for this notebook..

``` r
renv::snapshot()
sessionInfo()
```

render notebook to html for local rendered output check and to
github_document to check into github

``` r
rmarkdown::render("Manuscript_diff_feature_abundance_Joost.Rmd")
```
