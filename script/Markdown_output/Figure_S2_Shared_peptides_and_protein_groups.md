Supplementary figure S2 Shared peptides and protein groups between 1
ng-, 5 ng-, 10 ng-, and 100 ng tryptic HeLa
================
Samuel Grégoire

Supp figure S2 Shared peptides and protein groups between 1 ng-, 5 ng-,
10 ng-, and 100 ng tryptic HeLa.

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
LITDIA_5 <- read.csv2(paste0(lit_path, "20220502_LIT_Rapid_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_10 <- read.csv2(paste0(lit_path, "20220502_LIT_Rapid_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_100 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

# Replace protein accession by consensus

The same peptide sequence can sometime be attributed to a slightly
different Protein Group between samples. Since we will not focus on
proteoforms, we keep the most used Protein Group for the same peptide
sequence as a consensus.

``` r
LITDIA <- rbind(LITDIA_1, LITDIA_5, LITDIA_10, LITDIA_100)

LITDIA <- LITDIA[LITDIA$EG.IsImputed == "False",]

LITDIA <- LITDIA %>% 
  group_by(PEP.StrippedSequence) %>% 
  mutate(consensus_PG = names(sort(table(PG.ProteinGroups),
                              decreasing = TRUE))[1])
```

# Peptides

## Pull peptides from replicates

``` r
# 1 ng
pep_1_1 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_1_1)
```

    ## [1] 4533

``` r
pep_1_2 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s02") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_1_2)
```

    ## [1] 4543

``` r
pep_1_3 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s03") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_1_3)
```

    ## [1] 4540

``` r
pep_1_4 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s04") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_1_4)
```

    ## [1] 4533

``` r
# 5 ng
pep_5_1 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_5_1)
```

    ## [1] 7908

``` r
pep_5_2 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s02") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_5_2)
```

    ## [1] 7894

``` r
pep_5_3 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s03") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_5_3)
```

    ## [1] 7908

``` r
pep_5_4 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s04") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_5_4)
```

    ## [1] 7895

``` r
# 10 ng
pep_10_1 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01" ) %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_10_1)
```

    ## [1] 8728

``` r
pep_10_2 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01" ) %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_10_2)
```

    ## [1] 8728

``` r
pep_10_3 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_10_3)
```

    ## [1] 8728

``` r
pep_10_4 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_10_4)
```

    ## [1] 8728

``` r
# 100 ng
pep_100_1 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_100_1)
```

    ## [1] 9556

``` r
pep_100_2 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_100_2)
```

    ## [1] 9556

``` r
pep_100_3 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_100_3)
```

    ## [1] 9556

``` r
pep_100_4 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(PEP.StrippedSequence) %>% 
  unique()
length(pep_100_4)
```

    ## [1] 9556

## Filter peptides found in less than 3 replicates

``` r
pep_1 <- c(pep_1_1, pep_1_2, pep_1_3, pep_1_4)
pep_5 <- c(pep_5_1, pep_5_2, pep_5_3, pep_5_4)
pep_10 <- c(pep_10_1, pep_10_2, pep_10_3, pep_10_4)
pep_100 <- c(pep_100_1, pep_100_2, pep_100_3, pep_100_4)

pep_1_filter <- pep_1 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()

pep_5_filter <- pep_5 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()

pep_10_filter <-pep_10 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()

pep_100_filter <-pep_100 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()
```

## Venn diagram

``` r
pep_list <- list('1 ng' = unique(pep_1_filter),
     '5 ng' = unique(pep_5_filter),
     '10 ng' = unique(pep_10_filter),
     '100 ng' = unique(pep_100_filter))
str(pep_list)
```

    ## List of 4
    ##  $ 1 ng  : chr [1:4547] "AAAAELSLLEK" "AAAEVAGQFVIK" "AAAGEFADDPCSSVK" "AAAIGIDLGTTYSCVGVFQHGK" ...
    ##  $ 5 ng  : chr [1:7918] "AAAAELSLLEK" "AAAEAGGDDAR" "AAAEVAGQFVIK" "AAAEVNQDYGLDPK" ...
    ##  $ 10 ng : chr [1:8728] "AAAAELSLLEK" "AAADGDDSLYPIAVLIDELR" "AAAEAGGDDAR" "AAAEVAGQFVIK" ...
    ##  $ 100 ng: chr [1:9556] "AAAAAWEEPSSGNGTAR" "AAAASAAEAGIATTGTEDSDDALLK" "AAADTLQGPMQAAYR" "AAAEAGGDDAR" ...

``` r
pep_venn <- Venn(pep_list)
pep_data <- process_data(pep_venn)

pep_venn_plot <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(pep_data)) +
  geom_sf(aes(color = name), data = venn_setedge(pep_data), show.legend = FALSE, size = 1) +
  geom_sf_text(aes(label = c("1 ng", "5 ng", "10 ng", "100 ng")),
               nudge_x = c(0.06, 0.06, -0.06, -0.06),
               nudge_y = 0.02,
               data = venn_setlabel(pep_data),
               size = 5) +
  geom_sf_text(aes(label = count), 
                data = venn_region(pep_data),
                size = 4)+
  scale_fill_gradient(low = "#FFFFFF", high = "#2574A2")+
  scale_color_manual(values = rep("black", 4))+
  theme_void()+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  labs(fill = "Peptides")

pep_venn_plot
```

![](Figure_S2_Shared_peptides_and_protein_groups_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Proteins

## Pull proteins from replicates

``` r
# 1 ng
prot_1_1 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_1_1)
```

    ## [1] 1439

``` r
prot_1_2 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s02") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_1_2)
```

    ## [1] 1439

``` r
prot_1_3 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s03") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_1_3)
```

    ## [1] 1439

``` r
prot_1_4 <- LITDIA %>% 
  filter(R.FileName == "20220219_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Rapid_auto_s04") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_1_4)
```

    ## [1] 1441

``` r
# 5 ng
prot_5_1 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_5_1)
```

    ## [1] 2321

``` r
prot_5_2 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s02") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_5_2)
```

    ## [1] 2316

``` r
prot_5_3 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s03") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_5_3)
```

    ## [1] 2321

``` r
prot_5_4 <- LITDIA %>% 
  filter(R.FileName == "20220220_EV_Evo_Whisper100_40SPD_HeLa_5ng_1CV45_LITDIA_Rapid_auto_s04") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_5_4)
```

    ## [1] 2315

``` r
# 10 ng
prot_10_1 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01" ) %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_10_1)
```

    ## [1] 2429

``` r
prot_10_2 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01" ) %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_10_2)
```

    ## [1] 2429

``` r
prot_10_3 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_10_3)
```

    ## [1] 2429

``` r
prot_10_4 <- LITDIA %>% 
  filter(R.FileName == "20220221_EV_Evo_Whisper100_40SPD_HeLa_10ng_1CV45_LITDIA_Rapid_auto_s01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_10_4)
```

    ## [1] 2429

``` r
# 100 ng
prot_100_1 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_100_1)
```

    ## [1] 2957

``` r
prot_100_2 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_100_2)
```

    ## [1] 2957

``` r
prot_100_3 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_100_3)
```

    ## [1] 2957

``` r
prot_100_4 <- LITDIA %>% 
  filter(R.FileName == "20220415_HB_Evo_Whisper100_40SPD_LITDIA_HeLa_100ng_Rapid_1CV_auto_40w_S01") %>% 
  pull(consensus_PG) %>% 
  unique()
length(prot_100_4)
```

    ## [1] 2957

## Filter proteins found in less than 3 replicates

``` r
prot_1 <- c(prot_1_1, prot_1_2, prot_1_3, prot_1_4)
prot_5 <- c(prot_5_1, prot_5_2, prot_5_3, prot_5_4)
prot_10 <- c(prot_10_1, prot_10_2, prot_10_3, prot_10_4)
prot_100 <- c(prot_100_1, prot_100_2, prot_100_3, prot_100_4)

prot_1_filter <- prot_1 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()

prot_5_filter <- prot_5 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()

prot_10_filter <-prot_10 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()

prot_100_filter <- prot_100 %>% 
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 3) %>% 
  pull('.') %>%
  as.character()
```

## Venn diagram

``` r
prot_list <- list('1 ng' = unique(prot_1_filter),
                  '5 ng' = unique(prot_5_filter),
                  '10 ng' = unique(prot_10_filter),
                  '100 ng' = unique(prot_100_filter))

prot_venn <- Venn(prot_list)
prot_data <- process_data(prot_venn)

prot_venn_plot <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(prot_data)) +
  geom_sf(aes(color = name), data = venn_setedge(prot_data), show.legend = FALSE, size = 1) +
  geom_sf_text(aes(label = c("1 ng", "5 ng", "10 ng", "100 ng")),
               nudge_x = c(0.06, 0.06, -0.06, -0.06),
               nudge_y = 0.02,
               data = venn_setlabel(prot_data),
               size = 5) +
  geom_sf_text(aes(label = count), 
                data = venn_region(prot_data),
                size = 4)+
  scale_fill_gradient(low = "#FFFFFF", high = "#2574A2")+
  scale_color_manual(values = rep("black", 4))+
  theme_void()+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  labs(fill = "Proteins")

prot_venn_plot
```

![](Figure_S2_Shared_peptides_and_protein_groups_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# fig sup 2

``` r
fig_sup2 <- pep_venn_plot + prot_venn_plot

#save as pdf
#cairo_pdf(filename = "fig_sup2.pdf", width = 12.5, height = 8.33)
#fig_sup2
#dev.off()
```
