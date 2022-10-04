Supp. fig. S3 Distribution of the intensities of peptides or protein
groups from 100 ng of tryptic HeLa and overlapping between 1 ng and 100
ng tryptic HeLa
================
Samuel Grégoire

Supp. fig. S3 Distribution of the intensities of peptides (left) or
protein groups (right) from 100 ng of tryptic HeLa and overlapping of
identified peptides and protein groups between 1 ng and 100 ng tryptic
HeLa.

Comparison of the number of identified peptides on DIA-LIT-based method
on Normal scanning mode for WhisperTM 20 SPD and Rapid scanning mode for
WhisperTM 40 SPD from 1 ng of tryptic HeLa lysate with different numbers
of windows and IT at fixed cycle time.

# Load packages

``` r
library("tidyverse")
library("ggVennDiagram")
library("patchwork")
```

# Load data

The path to every folder used is stored to keep everything compact. LIT
data is stored in the folder located in lit_path. To reproduce analysis
the path must be changed to corresponding path on the local computer.

``` r
#lit_path <- "C:/Data/LIT/"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in RStudio with the function
`read.csv2()`.

``` r
LITDIA_1 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_100 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

# Replace protein accession by consensus

The same peptide sequence can sometime be attributed to a slightly
different Protein Group between samples. Since we will not focus on
proteoforms, we keep the most used Protein Group for the same peptide
sequence as a consensus.

``` r
LITDIA_1 <- LITDIA_1 %>% 
  group_by(PEP.StrippedSequence) %>% 
  mutate(consensus_PG = names(sort(table(PG.ProteinGroups),
                              decreasing = TRUE))[1])
LITDIA_100 <- LITDIA_100 %>% 
  group_by(PEP.StrippedSequence) %>% 
  mutate(consensus_PG = names(sort(table(PG.ProteinGroups),
                              decreasing = TRUE))[1])
```

# Clean data tables

Imported tables contain a lot of information which is not useful for our
analysis. The function below only keep the file names (replicates),
quantitatives values for each file name and the peptide sequence. It
also calculates the mean between replicates. Spectronaut does some
imputation by default but we want to keep missing values. Imputed values
are replaced by NA. Peptides missing in more than 1 replicates are
highlighted (TRUE) in the “filtered” column.

``` r
filter_pep <- function(x) {
  x$FG.Quantity[as.logical(x$EG.IsImputed)] <- NA
  y <- x %>% 
  group_by(R.FileName, PEP.StrippedSequence, consensus_PG) %>% 
  summarise(mean_quantity = mean(FG.Quantity, na.rm = TRUE)) %>%  
  pivot_wider(names_from = R.FileName, 
              values_from = mean_quantity) %>%
    suppressMessages()
  z <- y %>% 
  as.data.frame() %>% 
  mutate(na_nbr = rowSums(is.na(y[-c(1,2)]))) %>%
  mutate(mean_quantity = rowMeans(y[3:ncol(y)], na.rm = TRUE)) %>%
  mutate(filtered = !na_nbr <= 1) %>% 
  select(-na_nbr) %>% 
  suppressMessages("`summarise()` has grouped output by 'R.FileName'. You can override using the `.groups` argument.")
  message("Highlighted ", sum(z$filtered), " peptide(s) found in less than ", ncol(y) -3, " replicates")
  
 return(z) }
```

The filter_pep function is run on every dataset

``` r
LITDIA_1_filter <- filter_pep(LITDIA_1) %>% 
  filter(filtered == FALSE)
```

    ## Highlighted 12 peptide(s) found in less than 3 replicates

``` r
LITDIA_100_filter <- filter_pep(LITDIA_100) %>% 
  filter(filtered == FALSE)
```

    ## Highlighted 282 peptide(s) found in less than 3 replicates

# Highlight peptides found in 1 ng

``` r
LITDIA_100_filter$in_1 <- 
  LITDIA_100_filter$PEP.StrippedSequence %in%
  LITDIA_1_filter$PEP.StrippedSequence

LITDIA_100_filter$in_1 %>% table
```

    ## .
    ## FALSE  TRUE 
    ##  6319  3269

# Plot peptides

``` r
pep_overlap <- LITDIA_100_filter %>% 
  ggplot(aes(x = log2(mean_quantity), fill = in_1))+
  geom_histogram(bins = 40,
                 position = "identity",
                 alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(values = c("grey20", "#cd3700"),
                    labels = c("100 ng only", "100 ng and 1 ng"))+
  labs(x = "log2(FG.Quantity)",
       y = "Identified peptides")+
  theme(legend.position = c(0.783, 0.927),
        legend.background = element_blank(),
        axis.title = element_text(size = 27),
        axis.text = element_text(size = 27),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        strip.text.x = element_text(size = 27),
        strip.text.y = element_text(size = 27, vjust = 1.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.12))

pep_overlap
```

![](intensity_overlap_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Plot proteins

``` r
LITDIA_1_filter_prot <- LITDIA_1_filter %>% 
  group_by(consensus_PG) %>%
  summarise(mean_prot_quantity = mean(mean(mean_quantity)))

LITDIA_100_filter_prot <- LITDIA_100_filter %>% 
  group_by(consensus_PG) %>%
  summarise(mean_prot_quantity = mean(mean(mean_quantity)))

LITDIA_100_filter_prot$in_1 <- 
  LITDIA_100_filter_prot$consensus_PG %in%
  LITDIA_1_filter_prot$consensus_PG

LITDIA_100_filter_prot$in_1 %>% table
```

    ## .
    ## FALSE  TRUE 
    ##  1672  1100

``` r
prot_overlap <- LITDIA_100_filter_prot %>% 
  ggplot(aes(x = log2(mean_prot_quantity), fill = in_1))+
  geom_histogram(bins = 40,
                 position = "identity",
                 alpha = 0.4)+
  theme_bw()+
  scale_x_continuous(breaks = c(5, 10, 15, 20))+
  scale_fill_manual(values = c("grey20", "#cd3700"),
                    labels = c("100 ng only", "100 ng and 1 ng"))+
  labs(x = "log2(FG.Quantity)",
       y = "Identified proteins")+
  theme(legend.position = c(0.783, 0.927),
        legend.background = element_blank(),
        axis.title = element_text(size = 27),
        axis.text = element_text(size = 27),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        strip.text.x = element_text(size = 27),
        strip.text.y = element_text(size = 27, vjust = 1.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.12))
prot_overlap
```

![](intensity_overlap_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Join peptides and proteins

``` r
fig_sup3 <- pep_overlap + prot_overlap

#save as pdf
#cairo_pdf(filename = "fig_sup3.pdf", width = 20.83, height = 7.81)
#fig_sup3
#dev.off()
```
