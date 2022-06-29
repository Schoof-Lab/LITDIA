Fig. 1 Comparison between an Orbitrap Mass analyzer and Linear Ion Trap
Mass Analyzer on low-input proteomics
================
Lukas R. Woltereck & Samuel Grégoire

# Load packages

``` r
library("tidyverse")
library("ggpubr")
library("patchwork")
library("ggpointdensity")
#set bw theme as default theme
theme_set(theme_bw())
```

# Load data

The path to every folder is stored to keep everything compact. LIT data
is stored in the folder located in lit_path and OT data in ot_path. To
reproduce analysis the path must be changed to corresponding path on the
local computer.

``` r
lit_path <- "C://data/lit" 
ot_path <- "C://data/ot"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in R with the function `read.csv2()`.

LIT Turbo DIA

``` r
LITDIA_Turbo_1 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_5 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_10 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_100 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

LIT Rapid DIA

``` r
LITDIA_Rapid_1 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_5 <- read.csv2(paste0(lit_path, "20220502_LIT_Rapid_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_10 <- read.csv2(paste0(lit_path, "20220502_LIT_Rapid_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_100 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

LIT Normal DIA

``` r
LITDIA_Normal_1 <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_5 <- read.csv2(paste0(lit_path, "20220502_LIT_Normal_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_10 <- read.csv2(paste0(lit_path, "20220502_LIT_Normal_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_100 <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

OT 7.5k DIA

``` r
OTDIA_7k_1 <- read.csv2(paste0(ot_path, "20220421_OT_7k_DIA_1ng_40SPD_whisper100_1CV.csv"))
OTDIA_7k_5 <- read.csv2(paste0(ot_path, "20220502_OT_7k_DIA_5ng_40SPD_whisper100_1CV.csv"))
OTDIA_7k_10 <- read.csv2(paste0(ot_path, "20220502_OT_7k_DIA_10ng_40SPD_whisper100_1CV.csv"))
OTDIA_7k_100 <- read.csv2(paste0(ot_path, "20220421_OT_7k_DIA_100ng_40SPD_whisper100_1CV.csv"))
```

OT 15k DIA

``` r
OTDIA_15k_1 <- read.csv2(paste0(ot_path, "20220421_OT_15k_DIA_1ng_40SPD_whisper100_1CV.csv"))
OTDIA_15k_5 <- read.csv2(paste0(ot_path, "20220509_OT_15k_DIA_5ng_40SPD_whisper100_1CV.csv"))
OTDIA_15k_10 <- read.csv2(paste0(ot_path, "20220502_OT_15k_DIA_10ng_40SPD_whisper100_1CV.csv"))
OTDIA_15k_100 <- read.csv2(paste0(ot_path, "20220425_OT_15k_DIA_100ng_40SPD_whisper100_1CV.csv"))
```

OT 30k DIA

``` r
OTDIA_30k_1 <- read.csv2(paste0(ot_path, "20220421_OT_30k_DIA_1ng_40SPD_whisper100_1CV.csv"))
OTDIA_30k_5 <- read.csv2(paste0(ot_path, "20220509_OT_30k_DIA_5ng_40SPD_whisper100_1CV.csv"))
OTDIA_30k_10 <- read.csv2(paste0(ot_path, "20220502_OT_30k_DIA_10ng_40SPD_whisper100_1CV.csv"))
OTDIA_30k_100 <- read.csv2(paste0(ot_path, "20220425_OT_30k_DIA_100ng_40SPD_whisper100_1CV.csv"))
```

# Clean data tables

Imported tables contain a lot of information, which are not useful for
our analysis. The function `filter_pep()` only keeps the file names
(replicates), quantitative values for each file name, and the peptide
sequence. It also calculates the mean between the four replicates.
Spectronaut does some imputation by default, but we want to keep missing
values. Imputed values are replaced by NA. Peptides missing in more than
one replicate are highlighted as TRUE.

``` r
filter_pep <- function(x) {
  x$FG.Quantity[as.logical(x$EG.IsImputed)] <- NA
  y <- x %>% 
  group_by(R.FileName, PEP.StrippedSequence) %>% 
  summarise(mean_quantity = mean(FG.Quantity, na.rm = TRUE)) %>%  
  pivot_wider(names_from = R.FileName, 
              values_from = mean_quantity) %>%
    suppressMessages()
  z <- y %>% 
  mutate(na_nbr = rowSums(is.na(y[-1]))) %>%
  mutate(mean_quantity = rowMeans(y[2:ncol(y)], na.rm = TRUE)) %>%
  mutate(filtered = !na_nbr <= 1) %>% 
  select(-na_nbr) %>% 
  suppressMessages("`summarise()` has grouped output by 'R.FileName'. You can override using the `.groups` argument.")
  message("
Highlighted ", sum(z$filtered) , " peptide(s) found in less than ", ncol(y) -2, " replicates")
  
  return(z)
}
```

The `filter_pep()` function is used for every data-set. Furthermore a
“filter” data-set is created by removing all peptides missing in more
than one replicate.

1 ng

``` r
OTDIA_7k_1 <-  filter_pep(OTDIA_7k_1)
```

    ## 
    ## Highlighted 5 peptide(s) found in less than 3 replicates

``` r
OTDIA_7k_1_filter <- OTDIA_7k_1[OTDIA_7k_1$filtered == FALSE,]
OTDIA_15k_1 <-  filter_pep(OTDIA_15k_1)
```

    ## 
    ## Highlighted 7 peptide(s) found in less than 3 replicates

``` r
OTDIA_15k_1_filter <- OTDIA_15k_1[OTDIA_15k_1$filtered == FALSE,]
OTDIA_30k_1 <-  filter_pep(OTDIA_30k_1)
```

    ## 
    ## Highlighted 38 peptide(s) found in less than 3 replicates

``` r
OTDIA_30k_1_filter <- OTDIA_30k_1[OTDIA_30k_1$filtered == FALSE,]
LITDIA_Turbo_1 <- filter_pep(LITDIA_Turbo_1)
```

    ## 
    ## Highlighted 14 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_1_filter <- LITDIA_Turbo_1[LITDIA_Turbo_1$filtered == FALSE,]
LITDIA_Rapid_1 <- filter_pep(LITDIA_Rapid_1)
```

    ## 
    ## Highlighted 12 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_1_filter <- LITDIA_Rapid_1[LITDIA_Rapid_1$filtered == FALSE,]
LITDIA_Normal_1 <- filter_pep(LITDIA_Normal_1)
```

    ## 
    ## Highlighted 41 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_1_filter <- LITDIA_Normal_1[LITDIA_Normal_1$filtered == FALSE,]
```

5 ng

``` r
OTDIA_7k_5 <-  filter_pep(OTDIA_7k_5)
```

    ## 
    ## Highlighted 22 peptide(s) found in less than 3 replicates

``` r
OTDIA_7k_5_filter <- OTDIA_7k_5[OTDIA_7k_5$filtered == FALSE,]
OTDIA_15k_5 <-  filter_pep(OTDIA_15k_5)
```

    ## 
    ## Highlighted 88 peptide(s) found in less than 3 replicates

``` r
OTDIA_15k_5_filter <- OTDIA_15k_5[OTDIA_15k_5$filtered == FALSE,]
OTDIA_30k_5 <-  filter_pep(OTDIA_30k_5)
```

    ## 
    ## Highlighted 591 peptide(s) found in less than 3 replicates

``` r
OTDIA_30k_5_filter <- OTDIA_30k_5[OTDIA_30k_5$filtered == FALSE,]
LITDIA_Turbo_5 <- filter_pep(LITDIA_Turbo_5)
```

    ## 
    ## Highlighted 27 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_5_filter <- LITDIA_Turbo_5[LITDIA_Turbo_5$filtered == FALSE,]
LITDIA_Rapid_5 <- filter_pep(LITDIA_Rapid_5)
```

    ## 
    ## Highlighted 13 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_5_filter <- LITDIA_Rapid_5[LITDIA_Rapid_5$filtered == FALSE,]
LITDIA_Normal_5 <- filter_pep(LITDIA_Normal_5)
```

    ## 
    ## Highlighted 30 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_5_filter <- LITDIA_Normal_5[LITDIA_Normal_5$filtered == FALSE,]
```

10 ng

``` r
OTDIA_7k_10 <-  filter_pep(OTDIA_7k_10)
```

    ## 
    ## Highlighted 50 peptide(s) found in less than 3 replicates

``` r
OTDIA_7k_10_filter <- OTDIA_7k_10[OTDIA_7k_10$filtered == FALSE,]
OTDIA_15k_10 <-  filter_pep(OTDIA_15k_10)
```

    ## 
    ## Highlighted 46 peptide(s) found in less than 3 replicates

``` r
OTDIA_15k_10_filter <- OTDIA_15k_10[OTDIA_15k_10$filtered == FALSE,]
OTDIA_30k_10 <-  filter_pep(OTDIA_30k_10)
```

    ## 
    ## Highlighted 194 peptide(s) found in less than 3 replicates

``` r
OTDIA_30k_10_filter <- OTDIA_30k_10[OTDIA_30k_10$filtered == FALSE,]
LITDIA_Turbo_10 <- filter_pep(LITDIA_Turbo_10)
```

    ## 
    ## Highlighted 14 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_10_filter <- LITDIA_Turbo_10[LITDIA_Turbo_10$filtered == FALSE,]
LITDIA_Rapid_10 <- filter_pep(LITDIA_Rapid_10)
```

    ## 
    ## Highlighted 21 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_10_filter <- LITDIA_Rapid_10[LITDIA_Rapid_10$filtered == FALSE,]
LITDIA_Normal_10 <- filter_pep(LITDIA_Normal_10)
```

    ## 
    ## Highlighted 28 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_10_filter <- LITDIA_Normal_10[LITDIA_Normal_10$filtered == FALSE,]
```

100 ng

``` r
OTDIA_7k_100 <-  filter_pep(OTDIA_7k_100)
```

    ## 
    ## Highlighted 182 peptide(s) found in less than 3 replicates

``` r
OTDIA_7k_100_filter <- OTDIA_7k_100[OTDIA_7k_100$filtered == FALSE,]
OTDIA_15k_100 <-  filter_pep(OTDIA_15k_100)
```

    ## 
    ## Highlighted 132 peptide(s) found in less than 3 replicates

``` r
OTDIA_15k_100_filter <- OTDIA_15k_100[OTDIA_15k_100$filtered == FALSE,]
OTDIA_30k_100 <-  filter_pep(OTDIA_30k_100)
```

    ## 
    ## Highlighted 271 peptide(s) found in less than 3 replicates

``` r
OTDIA_30k_100_filter <- OTDIA_30k_100[OTDIA_30k_100$filtered == FALSE,]
LITDIA_Turbo_100 <- filter_pep(LITDIA_Turbo_100)
```

    ## 
    ## Highlighted 324 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_filter <- LITDIA_Turbo_100[LITDIA_Turbo_100$filtered == FALSE,]
LITDIA_Rapid_100 <- filter_pep(LITDIA_Rapid_100)
```

    ## 
    ## Highlighted 282 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_100_filter <- LITDIA_Rapid_100[LITDIA_Rapid_100$filtered == FALSE,]
LITDIA_Normal_100 <- filter_pep(LITDIA_Normal_100)
```

    ## 
    ## Highlighted 316 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_100_filter <- LITDIA_Normal_100[LITDIA_Normal_100$filtered == FALSE,]
```

# Figure A

Comparison of the number of identified peptides in serial dilution
(1-,5-,10-, and 100 ng) of HeLa tryptic digested between DIA-OT-based
methods with different resolution scans (7.5k, 15k, and 30k) and
DIA-LIT-based methods with different scanning modes (Turbo, Rapid,
Normal on 40 SPD).

## Number of identified peptides

Empty vectors are created to store the number of identified peptides in
each replicates for each data-set.

``` r
id_LITDIA_Turbo_1 <- id_LITDIA_Turbo_5 <- 
  id_LITDIA_Turbo_10 <- id_LITDIA_Turbo_100 <- 
  id_LITDIA_Rapid_1 <- id_LITDIA_Rapid_5 <- 
  id_LITDIA_Rapid_10 <- id_LITDIA_Rapid_100 <-
  id_LITDIA_Normal_1 <- id_LITDIA_Normal_5 <- 
  id_LITDIA_Normal_10 <- id_LITDIA_Normal_100 <-  
  id_OTDIA_7k_1 <- id_OTDIA_7k_5 <- 
  id_OTDIA_7k_10 <- id_OTDIA_7k_100 <- 
  id_OTDIA_15k_1 <- id_OTDIA_15k_5 <- 
  id_OTDIA_15k_10 <- id_OTDIA_15k_100 <- 
  id_OTDIA_30k_1 <- id_OTDIA_30k_5 <- 
  id_OTDIA_30k_10 <- id_OTDIA_30k_100 <- rep(NA, 4)
```

The numbers of identified peptides are counted by summing non-missing
values for each replicate and stored in the vectors created above.

### 1 ng

OT

``` r
for(i in 2:5){
  id_OTDIA_7k_1[i-1] <- sum(!is.na(OTDIA_7k_1[, i]))}
for(i in 2:5){
  id_OTDIA_15k_1[i-1] <- sum(!is.na(OTDIA_15k_1[, i]))}
for(i in 2:5){
  id_OTDIA_30k_1[i-1] <- sum(!is.na(OTDIA_30k_1[, i]))}
```

LIT

``` r
for(i in 2:5){
id_LITDIA_Turbo_1[i-1] <- sum(!is.na(LITDIA_Turbo_1[, i]))}
for(i in 2:5){
id_LITDIA_Rapid_1[i-1] <- sum(!is.na(LITDIA_Rapid_1[, i]))}
for(i in 2:5){
id_LITDIA_Normal_1[i-1] <- sum(!is.na(LITDIA_Normal_1[, i]))}
```

### 5 ng

OT

``` r
for(i in 2:5){
  id_OTDIA_7k_5[i-1] <- sum(!is.na(OTDIA_7k_5[, i]))}
for(i in 2:5){
  id_OTDIA_15k_5[i-1] <- sum(!is.na(OTDIA_15k_5[, i]))}
for(i in 2:5){
  id_OTDIA_30k_5[i-1] <- sum(!is.na(OTDIA_30k_5[, i]))}
```

LIT

``` r
for(i in 2:5){
id_LITDIA_Turbo_5[i-1] <- sum(!is.na(LITDIA_Turbo_5[, i]))}
for(i in 2:5){
id_LITDIA_Rapid_5[i-1] <- sum(!is.na(LITDIA_Rapid_5[, i]))}
for(i in 2:5){
id_LITDIA_Normal_5[i-1] <- sum(!is.na(LITDIA_Normal_5[, i]))}
```

### 10 ng

OT

``` r
for(i in 2:5){
  id_OTDIA_7k_10[i-1] <- sum(!is.na(OTDIA_7k_10[, i]))}
for(i in 2:5){
  id_OTDIA_15k_10[i-1] <- sum(!is.na(OTDIA_15k_10[, i]))}
for(i in 2:5){
  id_OTDIA_30k_10[i-1] <- sum(!is.na(OTDIA_30k_10[, i]))}
```

LIT Turbo DIA auto

``` r
for(i in 2:5){
id_LITDIA_Turbo_10[i-1] <- sum(!is.na(LITDIA_Turbo_10[, i]))}
for(i in 2:5){
id_LITDIA_Rapid_10[i-1] <- sum(!is.na(LITDIA_Rapid_10[, i]))}
for(i in 2:5){
id_LITDIA_Normal_10[i-1] <- sum(!is.na(LITDIA_Normal_10[, i]))}
```

### 100 ng

OT

``` r
for(i in 2:5){
  id_OTDIA_7k_100[i-1] <- sum(!is.na(OTDIA_7k_100[, i]))}
for(i in 2:5){
  id_OTDIA_15k_100[i-1] <- sum(!is.na(OTDIA_15k_100[, i]))}
for(i in 2:5){
  id_OTDIA_30k_100[i-1] <- sum(!is.na(OTDIA_30k_100[, i]))}
```

LIT

``` r
for(i in 2:5){
id_LITDIA_Turbo_100[i-1] <- sum(!is.na(LITDIA_Turbo_100[, i]))}
for(i in 2:5){
id_LITDIA_Rapid_100[i-1] <- sum(!is.na(LITDIA_Rapid_100[, i]))}
for(i in 2:5){
id_LITDIA_Normal_100[i-1] <- sum(!is.na(LITDIA_Normal_100[, i]))}
```

## Join A

Number of identified peptides for each files are combined into a data
frame containing information about the method (OT resolution and LIT
scanning mode), the replicate and the input. Method and input variables
are factorized to set a relevant order in the plot instead of the
alphabetical order.

``` r
id_pep <- data.frame(method = rep(c(rep("OT 7.5k", 4), rep("OT 15k", 4), rep("OT 30k", 4),
                                    rep("LIT Turbo", 4), rep("LIT Rapid", 4), rep("LIT Normal", 4)), 4),
           replicate = rep(1:4, 24),
           id = c(id_OTDIA_7k_1, id_OTDIA_15k_1, id_OTDIA_30k_1, 
                  id_LITDIA_Turbo_1, id_LITDIA_Rapid_1, id_LITDIA_Normal_1,
                  id_OTDIA_7k_5, id_OTDIA_15k_5, id_OTDIA_30k_5, 
                  id_LITDIA_Turbo_5, id_LITDIA_Rapid_5, id_LITDIA_Normal_5,
                  id_OTDIA_7k_10, id_OTDIA_15k_10, id_OTDIA_30k_10, 
                  id_LITDIA_Turbo_10, id_LITDIA_Rapid_10, id_LITDIA_Normal_10,
                  id_OTDIA_7k_100, id_OTDIA_15k_100, id_OTDIA_30k_100, 
                  id_LITDIA_Turbo_100, id_LITDIA_Rapid_100, id_LITDIA_Normal_100),
           input = c(rep("1 ng", 24), rep("5 ng", 24), rep("10 ng", 24), rep("100 ng", 24)))
# Factorise to set order
id_pep$method <- factor(id_pep$method, levels = unique(id_pep$method))
id_pep$input <- factor(id_pep$input, levels = unique(id_pep$input))
```

## Plot A

The mean and standard deviation is calculated for the replicates of each
method. Both are plotted in a line plot using `ggplot2`. Error bars
representing the standard deviation are displayed using
`geom_errorbar()`.

``` r
figure1a <- id_pep %>% 
  group_by(method, input) %>% 
  summarise(mean_id = mean(id),
            sd_id = sd(id))%>%
  ggplot(aes(x = input, y = mean_id, group = method, color = method))+
  geom_line(size = 2.5)+
  geom_point(size = 4)+
  scale_color_manual(values = c("#ADC2E8", "#7F9FDB", "#527DCE", "#EEB9CF", "#DF7FA7", "#D35087"))+
  geom_errorbar(aes(ymin = mean_id - sd_id, 
                         ymax = mean_id + sd_id),
                    width = 0.1,
                    size = 1.5)+
  labs(color = "Methods",
       x = "Input",
       y = "Identified peptides",
       title = "A")+
  theme(legend.position = c(0.23, 0.85),
        legend.background = element_blank(),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 22),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 27),
        legend.key.size = unit(0.9, "cm"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.2),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
figure1a
```

![](Fig.-1-Comparison-between-an-Orbitrap-Mass-analyzer-and-Linear-Ion-Trap-Mass-Analyzer-on-low-input-proteomics_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

# Figure B

Comparison of identified peptides between the DIA-OT-based methods and
the DIA-LIT-based methods. Identified peptides with a coefficient of
variation (CV) between 10% and 15% are colored with light red and those
with a CV below 10 % with dark red.

## Compute CV

The function `compute_CV()` below takes the data-sets, already filtered
by `filter_pep()`, as an entry and compute the coefficient of variation
(CV) for the replicates of each method. Most of our data-sets include
four replicates but some contain only three. In this case, an empty
column is added as a fourth replicate to keep the same format for every
data-set.

``` r
compute_CV <- function(x){
  y <- x %>% 
    mutate(sd = apply(x[, 2:(ncol(x)-2)],
                      FUN = sd, MARGIN =  1, na.rm = TRUE),
           CV = sd/mean_quantity)
  colnames(y)[2:(ncol(x)-2)] <- paste0("S", 1:(ncol(x)-3))
  if(ncol(x) <= 6){
  y$S4 <- NA
 }
  return(y)
}
```

Information about the method and the input are included afterwards.

### 1 ng

LIT

``` r
LITDIA_Turbo_1_CV <- compute_CV(LITDIA_Turbo_1)
LITDIA_Turbo_1_CV$method <- "LIT Turbo"
LITDIA_Turbo_1_CV$input<- "1 ng"
LITDIA_Rapid_1_CV <- compute_CV(LITDIA_Rapid_1)
LITDIA_Rapid_1_CV$method <- "LIT Rapid"
LITDIA_Rapid_1_CV$input<- "1 ng"
LITDIA_Normal_1_CV <- compute_CV(LITDIA_Normal_1)
LITDIA_Normal_1_CV$method <- "LIT Normal"
LITDIA_Normal_1_CV$input<- "1 ng"
```

OT

``` r
OTDIA_7k_1_CV <- compute_CV(OTDIA_7k_1)
OTDIA_7k_1_CV$method <- "OT 7k"
OTDIA_7k_1_CV$input<- "1 ng"
OTDIA_15k_1_CV <- compute_CV(OTDIA_15k_1)
OTDIA_15k_1_CV$method <- "OT 15k"
OTDIA_15k_1_CV$input<- "1 ng"
OTDIA_30k_1_CV <- compute_CV(OTDIA_30k_1)
OTDIA_30k_1_CV$method <- "OT 30k"
OTDIA_30k_1_CV$input<- "1 ng"
```

### 5 ng

LIT

``` r
LITDIA_Turbo_5_CV <- compute_CV(LITDIA_Turbo_5)
LITDIA_Turbo_5_CV$method <- "LIT Turbo"
LITDIA_Turbo_5_CV$input<- "5 ng"
LITDIA_Rapid_5_CV <- compute_CV(LITDIA_Rapid_5)
LITDIA_Rapid_5_CV$method <- "LIT Rapid"
LITDIA_Rapid_5_CV$input<- "5 ng"
LITDIA_Normal_5_CV <- compute_CV(LITDIA_Normal_5)
LITDIA_Normal_5_CV$method <- "LIT Normal"
LITDIA_Normal_5_CV$input<- "5 ng"
```

OT

``` r
OTDIA_7k_5_CV <- compute_CV(OTDIA_7k_5)
OTDIA_7k_5_CV$method <- "OT 7k"
OTDIA_7k_5_CV$input<- "5 ng"
OTDIA_15k_5_CV <- compute_CV(OTDIA_15k_5)
OTDIA_15k_5_CV$method <- "OT 15k"
OTDIA_15k_5_CV$input<- "5 ng"
OTDIA_30k_5_CV <- compute_CV(OTDIA_30k_5)
OTDIA_30k_5_CV$method <- "OT 30k"
OTDIA_30k_5_CV$input<- "5 ng"
```

### 10 ng

LIT

``` r
LITDIA_Turbo_10_CV <- compute_CV(LITDIA_Turbo_10)
LITDIA_Turbo_10_CV$method <- "LIT Turbo"
LITDIA_Turbo_10_CV$input<- "10 ng"
LITDIA_Rapid_10_CV <- compute_CV(LITDIA_Rapid_10)
LITDIA_Rapid_10_CV$method <- "LIT Rapid"
LITDIA_Rapid_10_CV$input<- "10 ng"
LITDIA_Normal_10_CV <- compute_CV(LITDIA_Normal_10)
LITDIA_Normal_10_CV$method <- "LIT Normal"
LITDIA_Normal_10_CV$input<- "10 ng"
```

OT

``` r
OTDIA_7k_10_CV <- compute_CV(OTDIA_7k_10)
OTDIA_7k_10_CV$method <- "OT 7k"
OTDIA_7k_10_CV$input<- "10 ng"
OTDIA_15k_10_CV <- compute_CV(OTDIA_15k_10)
OTDIA_15k_10_CV$method <- "OT 15k"
OTDIA_15k_10_CV$input<- "10 ng"
OTDIA_30k_10_CV <- compute_CV(OTDIA_30k_10)
OTDIA_30k_10_CV$method <- "OT 30k"
OTDIA_30k_10_CV$input<- "10 ng"
```

### 100 ng

LIT

``` r
LITDIA_Turbo_100_CV <- compute_CV(LITDIA_Turbo_100)
LITDIA_Turbo_100_CV$method <- "LIT Turbo"
LITDIA_Turbo_100_CV$input<- "100 ng"
LITDIA_Rapid_100_CV <- compute_CV(LITDIA_Rapid_100)
LITDIA_Rapid_100_CV$method <- "LIT Rapid"
LITDIA_Rapid_100_CV$input<- "100 ng"
LITDIA_Normal_100_CV <- compute_CV(LITDIA_Normal_100)
LITDIA_Normal_100_CV$method <- "LIT Normal"
LITDIA_Normal_100_CV$input<- "100 ng"
```

OT

``` r
OTDIA_7k_100_CV <- compute_CV(OTDIA_7k_100)
OTDIA_7k_100_CV$method <- "OT 7k"
OTDIA_7k_100_CV$input<- "100 ng"
OTDIA_15k_100_CV <- compute_CV(OTDIA_15k_100)
OTDIA_15k_100_CV$method <- "OT 15k"
OTDIA_15k_100_CV$input<- "100 ng"
OTDIA_30k_100_CV <- compute_CV(OTDIA_30k_100)
OTDIA_30k_100_CV$method <- "OT 30k"
OTDIA_30k_100_CV$input<- "100 ng"
```

## Join B

Data frames containing CV are combined using `rbind()`.

``` r
CV_joined <- 
  rbind(LITDIA_Turbo_1_CV, LITDIA_Rapid_1_CV, LITDIA_Normal_1_CV,
      OTDIA_7k_1_CV, OTDIA_15k_1_CV, OTDIA_30k_1_CV,
      LITDIA_Turbo_5_CV, LITDIA_Rapid_5_CV, LITDIA_Normal_5_CV,
      OTDIA_7k_5_CV, OTDIA_15k_5_CV, OTDIA_30k_5_CV,
      LITDIA_Turbo_10_CV, LITDIA_Rapid_10_CV, LITDIA_Normal_10_CV,
      OTDIA_7k_10_CV, OTDIA_15k_10_CV, OTDIA_30k_10_CV,
      LITDIA_Turbo_100_CV, LITDIA_Rapid_100_CV, LITDIA_Normal_100_CV,
      OTDIA_7k_100_CV, OTDIA_15k_100_CV, OTDIA_30k_100_CV)
CV_joined$method %>% table()
```

    ## .
    ## LIT Normal  LIT Rapid  LIT Turbo     OT 15k     OT 30k      OT 7k 
    ##      34020      31156      27124      19463      33600      12130

``` r
CV_joined$input %>% table()
```

    ## .
    ##   1 ng  10 ng 100 ng   5 ng 
    ##  15351  39034  69104  34004

``` r
# Factorise to set order
CV_joined$method <- factor(CV_joined$method, 
                          levels = unique(CV_joined$method))
CV_joined$input <- factor(CV_joined$input, 
                          levels = unique(CV_joined$input))
```

A new column called “CV_filter” is added to identify peptides with a CV
over 0.15 (1), between 0.1 and 0.15 (2) and under 0.1 (3).

``` r
CV_joined$CV_filter <- NA
CV_joined$CV_filter[CV_joined$CV < 0.1] <- "3"
CV_joined$CV_filter[CV_joined$CV >= 0.1 & CV_joined$CV <= 0.15] <- "2"
CV_joined$CV_filter[CV_joined$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out. The
number of identified peptides for each replicate and each CV_filter is
calculated by summing the number of non-missing peptides. The mean and
the standard deviation of the replicates are also calculated for each
method.

``` r
CV_joined <- CV_joined %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(method, input, CV_filter) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
CV_joined$mean_id <- rowMeans(CV_joined[4:7], na.rm = TRUE)
CV_joined$sd_id <- apply(CV_joined[4:7], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

## Plot B

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

``` r
CV_joined$y_errorbar <- CV_joined$mean_id
CV_joined$y_errorbar[CV_joined$CV_filter == "1"] <- 
  CV_joined$y_errorbar[CV_joined$CV_filter == "1"] +
  CV_joined$y_errorbar[CV_joined$CV_filter == "2"]+
  CV_joined$y_errorbar[CV_joined$CV_filter == "3"]
CV_joined$y_errorbar[CV_joined$CV_filter == "2"] <-
  CV_joined$y_errorbar[CV_joined$CV_filter == "2"]+
  CV_joined$y_errorbar[CV_joined$CV_filter == "3"]
figure1b <- CV_joined %>% 
  ggplot(aes(x = input, y = mean_id, fill = CV_filter))+
    geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)+
  geom_col(width = 0.7) +
  facet_wrap(~method) +
  labs(x = "Input",
       y = "Identified peptides",
       fill = "CV",
       title = "B") +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 27),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.title.y = element_text(vjust = 2.2),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 22),
        legend.key.size = unit(0.9, "cm"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.08),
        plot.margin = margin(l = 25),
        legend.position = c(0.1, 0.285),
        legend.background = element_blank())+
  scale_fill_manual(values = c("grey82", "#ff9090", "#ff3030"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 18000))
 
figure1b
```

![](Fig.-1-Comparison-between-an-Orbitrap-Mass-analyzer-and-Linear-Ion-Trap-Mass-Analyzer-on-low-input-proteomics_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

# Figure C

Showing the pearson correlation of identified peptides between
DIA-OT-based methods and DIA-LIT-based methods, with different input (1
ng, 5 ng, 10 ng and 100 ng). The HeLa lysate was analyzed in quadruple.

## Join LIT and OT data frames

New data frames are created with LIT and OT quantification data for
peptides found in both methods. LIT Turbo is compared with OT 7.5k,
Rapid with 15k and Normal with 30k. Only peptides found in at least 3 of
the 4 replicates are kept. To complete this data frame, information
about the employed method and the input are included.

### 1 ng

``` r
Turbo_vs_7k_1 <- 
  inner_join(OTDIA_7k_1_filter[,c(1,6)],
             LITDIA_Turbo_1_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Turbo_vs_7k_1$method <- "LIT Turbo vs OT 7.5k"
Turbo_vs_7k_1$input <- "1 ng"
Rapid_vs_15k_1 <-  
  inner_join(OTDIA_15k_1_filter[,c(1,6)],
             LITDIA_Rapid_1_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Rapid_vs_15k_1$method <- "LIT Rapid vs OT 15k"
Rapid_vs_15k_1$input <- "1 ng"
Normal_vs_30k_1 <-  
  inner_join(OTDIA_30k_1_filter[,c(1,6)],
             LITDIA_Normal_1_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Normal_vs_30k_1$method <- "LIT Normal vs OT 30k"
Normal_vs_30k_1$input <- "1 ng"
```

### 5 ng

``` r
Turbo_vs_7k_5 <- 
  inner_join(OTDIA_7k_5_filter[,c(1,6)],
             LITDIA_Turbo_5_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Turbo_vs_7k_5$method <- "LIT Turbo vs OT 7.5k"
Turbo_vs_7k_5$input <- "5 ng"
Rapid_vs_15k_5 <-  
  inner_join(OTDIA_15k_5_filter[,c(1,6)],
             LITDIA_Rapid_5_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Rapid_vs_15k_5$method <- "LIT Rapid vs OT 15k"
Rapid_vs_15k_5$input <- "5 ng"
Normal_vs_30k_5 <-  
  inner_join(OTDIA_30k_5_filter[,c(1,6)],
             LITDIA_Normal_5_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Normal_vs_30k_5$method <- "LIT Normal vs OT 30k"
Normal_vs_30k_5$input <- "5 ng"
```

### 10 ng

``` r
Turbo_vs_7k_10 <- 
  inner_join(OTDIA_7k_10_filter[,c(1,6)],
             LITDIA_Turbo_10_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Turbo_vs_7k_10$method <- "LIT Turbo vs OT 7.5k"
Turbo_vs_7k_10$input <- "10 ng"
Rapid_vs_15k_10 <-  
  inner_join(OTDIA_15k_10_filter[,c(1,6)],
             LITDIA_Rapid_10_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Rapid_vs_15k_10$method <- "LIT Rapid vs OT 15k"
Rapid_vs_15k_10$input <- "10 ng"
Normal_vs_30k_10 <-  
  inner_join(OTDIA_30k_10_filter[,c(1,6)],
             LITDIA_Normal_10_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Normal_vs_30k_10$method <- "LIT Normal vs OT 30k"
Normal_vs_30k_10$input <- "10 ng"
```

### 100 ng

``` r
Turbo_vs_7k_100 <-  
  inner_join(OTDIA_7k_100_filter[,c(1,6)],
             LITDIA_Turbo_100_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Turbo_vs_7k_100$method <- "LIT Turbo vs OT 7.5k"
Turbo_vs_7k_100$input <- "100 ng"
Rapid_vs_15k_100 <- 
  inner_join(OTDIA_15k_100_filter[,c(1,6)],
             LITDIA_Rapid_100_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Rapid_vs_15k_100$method <- "LIT Rapid vs OT 15k"
Rapid_vs_15k_100$input <- "100 ng"
Normal_vs_30k_100 <- 
  inner_join(OTDIA_30k_100_filter[,c(1,6)],
             LITDIA_Normal_100_filter[,c(1,6)],
             by = "PEP.StrippedSequence",
             suffix = c("_OT", "_LIT"))
Normal_vs_30k_100$method <- "LIT Normal vs OT 30k"
Normal_vs_30k_100$input <- "100 ng"
```

## Join C

The data frames are combined with the function `rbind()`.

``` r
figure1c_join <- rbind(Turbo_vs_7k_1, Rapid_vs_15k_1, Normal_vs_30k_1,
               Turbo_vs_7k_5, Rapid_vs_15k_5, Normal_vs_30k_5,
               Turbo_vs_7k_10, Rapid_vs_15k_10, Normal_vs_30k_10,
               Turbo_vs_7k_100, Rapid_vs_15k_100, Normal_vs_30k_100)
# Factorise to set order
figure1c_join$method <- factor(figure1c_join$method,
                       levels = c("LIT Turbo vs OT 7.5k",
                                  "LIT Rapid vs OT 15k",
                                  "LIT Normal vs OT 30k"))
figure1c_join$input <- factor(figure1c_join$input, 
                              levels = unique(figure1c_join$input))
```

## Plot C

LIT and OT abundances are plotted against each other in a scatter plot
using `ggplot2`. The functions `geom_smooth()`and `stat_cor()` were used
to add a regression line and its coefficient of correlation to the plot.
Numbers of peptides found in both LIT and OT are also displayed using
`geom_text`. Logarithmic scale was used for both x and y axis.

``` r
figure1c <- figure1c_join %>% 
  ggplot(aes(x = mean_quantity_OT, y = mean_quantity_LIT)) +
  geom_pointdensity(size = 1.5)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.3),
        axis.title = element_text(size = 27),
        axis.text = element_text(size = 24),
        strip.text.x = element_text(size = 27),
        strip.text.y = element_text(size = 27, vjust = 1.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.12)) +
  geom_smooth(method = lm,
              se = FALSE,
              color = "red",
              size = 1.2)+
  stat_cor(method = "pearson",
           aes(label = ..r.label..),
           r.digits = 3, color = "red",
           size = 7,
           label.x.npc = 0,
           label.y.npc = 0.9)+
  labs(x = "OT Abundance (log10)",
       y = "LIT Abundance (log10)",
       title = "C") +
  scale_y_log10(breaks = c(1000, 100000)) +
  scale_x_log10() +
  facet_grid(input~method)
labels <- 
  paste0(c(sum(figure1c_join$method == "LIT Turbo vs OT 7.5k" &
                 figure1c_join$input == "1 ng"),
           sum(figure1c_join$method == "LIT Rapid vs OT 15k" &
                 figure1c_join$input == "1 ng"),
           sum(figure1c_join$method == "LIT Normal vs OT 30k" &
                 figure1c_join$input == "1 ng"),
           sum(figure1c_join$method == "LIT Turbo vs OT 7.5k" &
                 figure1c_join$input == "5 ng"),
           sum(figure1c_join$method == "LIT Rapid vs OT 15k" & 
                 figure1c_join$input == "5 ng"),
           sum(figure1c_join$method == "LIT Normal vs OT 30k" &
                 figure1c_join$input == "5 ng"),
           sum(figure1c_join$method == "LIT Turbo vs OT 7.5k" &
                 figure1c_join$input == "10 ng"),
           sum(figure1c_join$method == "LIT Rapid vs OT 15k" &
                 figure1c_join$input == "10 ng"),
           sum(figure1c_join$method == "LIT Normal vs OT 30k" &
                 figure1c_join$input == "10 ng"),
           sum(figure1c_join$method == "LIT Turbo vs OT 7.5k" & 
                 figure1c_join$input == "100 ng"),
           sum(figure1c_join$method == "LIT Rapid vs OT 15k" &
                 figure1c_join$input == "100 ng"),
           sum(figure1c_join$method == "LIT Normal vs OT 30k" &
                 figure1c_join$input == "100 ng")),
         " pep.")
data_labels <- 
  data.frame(method = as.factor(rep(c("LIT Turbo vs OT 7.5k",
                                      "LIT Rapid vs OT 15k",
                                      "LIT Normal vs OT 30k"),4)),
             input = as.factor(c(rep("1 ng",3), rep("5 ng",3),rep("10 ng",3), rep("100 ng",3))),
             label = labels)
figure1c <- figure1c +
  geom_text(x = -Inf, y = Inf,
            hjust = -0.15, vjust = 3.2,
            aes(label = labels, group = NULL), data = data_labels,
            size = 7,
            face = "bold")
figure1c
```

![](Fig.-1-Comparison-between-an-Orbitrap-Mass-analyzer-and-Linear-Ion-Trap-Mass-Analyzer-on-low-input-proteomics_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

# Figure D

Distribution of the range of detection between the DIA-OT-based method
and DIA-LIT-based methods and overlapping of identified peptides based
on their intensities.

## Highlight peptides found in OT

A new column is added to each LIT data-set to identify peptides which
were also found in their corresponding OT data-sets. Information about
the scanning mode and input are added.

### 1 ng

``` r
#Turbo
LITDIA_Turbo_1_filter$in_OT <- 
  LITDIA_Turbo_1_filter$PEP.StrippedSequence %in%
  OTDIA_7k_1_filter$PEP.StrippedSequence
LITDIA_Turbo_1_filter$method <- "Turbo"
#Rapid
LITDIA_Rapid_1_filter$in_OT <- 
  LITDIA_Rapid_1_filter$PEP.StrippedSequence %in%
  OTDIA_15k_1_filter$PEP.StrippedSequence
LITDIA_Rapid_1_filter$method <- "Rapid"
#Normal
LITDIA_Normal_1_filter$in_OT <- 
  LITDIA_Normal_1_filter$PEP.StrippedSequence %in%
  OTDIA_30k_1_filter$PEP.StrippedSequence
LITDIA_Normal_1_filter$method <- "Normal"
#Input
LITDIA_Turbo_1_filter$input <-
  LITDIA_Rapid_1_filter$input <-
  LITDIA_Normal_1_filter$input <- "1 ng"
```

### 5 ng

``` r
#Turbo
LITDIA_Turbo_5_filter$in_OT <- 
  LITDIA_Turbo_5_filter$PEP.StrippedSequence %in%
  OTDIA_7k_5_filter$PEP.StrippedSequence
LITDIA_Turbo_5_filter$method <- "Turbo"
#Rapid
LITDIA_Rapid_5_filter$in_OT <- 
  LITDIA_Rapid_5_filter$PEP.StrippedSequence %in%
  OTDIA_15k_5_filter$PEP.StrippedSequence
LITDIA_Rapid_5_filter$method <- "Rapid"
#Normal
LITDIA_Normal_5_filter$in_OT <- 
  LITDIA_Normal_5_filter$PEP.StrippedSequence %in%
  OTDIA_30k_5_filter$PEP.StrippedSequence
LITDIA_Normal_5_filter$method <- "Normal"
#Input
LITDIA_Turbo_5_filter$input <-
  LITDIA_Rapid_5_filter$input <-
  LITDIA_Normal_5_filter$input <- "5 ng"
```

### 10 ng

``` r
#Turbo
LITDIA_Turbo_10_filter$in_OT <- 
  LITDIA_Turbo_10_filter$PEP.StrippedSequence %in%
  OTDIA_7k_10_filter$PEP.StrippedSequence
LITDIA_Turbo_10_filter$method <- "Turbo"
#Rapid
LITDIA_Rapid_10_filter$in_OT <- 
  LITDIA_Rapid_10_filter$PEP.StrippedSequence %in%
  OTDIA_15k_10_filter$PEP.StrippedSequence
LITDIA_Rapid_10_filter$method <- "Rapid"
#Normal
LITDIA_Normal_10_filter$in_OT <- 
  LITDIA_Normal_10_filter$PEP.StrippedSequence %in%
  OTDIA_30k_10_filter$PEP.StrippedSequence
LITDIA_Normal_10_filter$method <- "Normal"
#Input
LITDIA_Turbo_10_filter$input <-
  LITDIA_Rapid_10_filter$input <-
  LITDIA_Normal_10_filter$input <- "10 ng"
```

### 100 ng

``` r
#Turbo
LITDIA_Turbo_100_filter$in_OT <- 
  LITDIA_Turbo_100_filter$PEP.StrippedSequence %in%
  OTDIA_7k_100_filter$PEP.StrippedSequence
LITDIA_Turbo_100_filter$method <- "Turbo"
#Rapid
LITDIA_Rapid_100_filter$in_OT <- 
  LITDIA_Rapid_100_filter$PEP.StrippedSequence %in%
  OTDIA_15k_100_filter$PEP.StrippedSequence
LITDIA_Rapid_100_filter$method <- "Rapid"
#Normal
LITDIA_Normal_100_filter$in_OT <- 
  LITDIA_Normal_100_filter$PEP.StrippedSequence %in%
  OTDIA_30k_100_filter$PEP.StrippedSequence
LITDIA_Normal_100_filter$method <- "Normal"
#Input
LITDIA_Turbo_100_filter$input <-
  LITDIA_Rapid_100_filter$input <- 
  LITDIA_Normal_100_filter$input <- "100 ng"
```

## Join D

The data frames are combined with the function `rbind()`. Only the
relevant columns are kept.

``` r
figure1d_join <- 
  rbind(LITDIA_Turbo_1_filter[,c(1,6:10)], LITDIA_Turbo_5_filter[,c(1,6:10)],
        LITDIA_Turbo_10_filter[,c(1,6:10)], LITDIA_Turbo_100_filter[,c(1,6:10)],
        LITDIA_Rapid_1_filter[,c(1,6:10)], LITDIA_Rapid_5_filter[,c(1,6:10)],
        LITDIA_Rapid_10_filter[,c(1,6:10)], LITDIA_Rapid_100_filter[,c(1,6:10)],
        LITDIA_Normal_1_filter[,c(1,6:10)], LITDIA_Normal_5_filter[,c(1,6:10)],
        LITDIA_Normal_10_filter[,c(1,6:10)], LITDIA_Normal_100_filter[,c(1,6:10)]) 
# Factorise to set order
figure1d_join$method <- factor(figure1d_join$method,
                       levels = unique(figure1d_join$method))
figure1d_join$input <- factor(figure1d_join$input,
                       levels = unique(figure1d_join$input))
```

## Plot D

LIT abundances for peptides found in LIT only (grey), or in both LIT and
OT (red), are plotted in a histogram using `ggplot2`. Numbers of
peptides found in LIT only (grey) and in both (LIT and OT, red) are
displayed using `geom_text`. Logarithmic scale was used for x axis.

``` r
figure1d <- figure1d_join %>% 
  ggplot(aes(x = mean_quantity, fill = in_OT))+
  geom_histogram(bins = 40,
                 position = "identity",
                 alpha = 0.4)+
  scale_x_log10()+
  scale_y_continuous(breaks = c(0, 500, 1000)) +
  scale_fill_manual(values = c("grey20", "red"),
                    labels = c("LIT only", "LIT and OT"))+
  labs(x = "DIA-LIT Abundance (log10)",
       y = "Identified peptides",
       title = "D")+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.3),
        legend.position = c(0.915, 0.927),
        legend.background = element_blank(),
        axis.title = element_text(size = 27),
        axis.text = element_text(size = 27),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        strip.text.x = element_text(size = 27),
        strip.text.y = element_text(size = 27, vjust = 1.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.12)) +
facet_grid(input ~ method)
labels1 <- 
  paste0(c(sum(figure1d_join$method == "Turbo" & 
                 figure1d_join$input == "1 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Rapid" & 
                 figure1d_join$input == "1 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Normal" & 
                 figure1d_join$input == "1 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Turbo" & 
                 figure1d_join$input == "5 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Rapid" & 
                 figure1d_join$input == "5 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Normal" & 
                 figure1d_join$input == "5 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Turbo" & 
                 figure1d_join$input == "10 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Rapid" & 
                 figure1d_join$input == "10 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Normal" & 
                 figure1d_join$input == "10 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Turbo" & 
                 figure1d_join$input == "100 ng" & 
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Rapid" &
                 figure1d_join$input == "100 ng" &
                 figure1d_join$in_OT == FALSE),
           sum(figure1d_join$method == "Normal" &
                 figure1d_join$input == "100 ng" &
                 figure1d_join$in_OT == FALSE)),
         " pep.")
                       
data_labels1 <- 
  data.frame(method = as.factor(rep(c("Turbo", "Rapid", "Normal"),4)),
             input = as.factor(c(rep("1 ng",3), rep("5 ng",3),
                                 rep("10 ng",3), rep("100 ng",3))),
             label = labels1,
             in_OT = rep(FALSE, 12))
labels2<- 
  paste0(c(sum(figure1d_join$method == "Turbo" &
                 figure1d_join$input == "1 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Rapid" &
                 figure1d_join$input == "1 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Normal" &
                 figure1d_join$input == "1 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Turbo" &
                 figure1d_join$input == "5 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Rapid" &
                 figure1d_join$input == "5 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Normal" &
                 figure1d_join$input == "5 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Turbo" &
                 figure1d_join$input == "10 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Rapid" &
                 figure1d_join$input == "10 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Normal" &
                 figure1d_join$input == "10 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Turbo" &
                 figure1d_join$input == "100 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Rapid" &
                 figure1d_join$input == "100 ng" &
                 figure1d_join$in_OT == TRUE),
           sum(figure1d_join$method == "Normal" &
                 figure1d_join$input == "100 ng" &
                 figure1d_join$in_OT == TRUE)),
         " pep.")
                       
data_labels2 <- 
  data.frame(method = as.factor(rep(c("Turbo", "Rapid", "Normal"),4)),
             input = as.factor(c(rep("1 ng",3), rep("5 ng",3),
                                 rep("10 ng",3), rep("100 ng",3))),
             label = labels2,
             in_OT = rep(TRUE, 12))
figure1d <- figure1d +
  geom_text(x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.82,
            aes(label = label, group = NULL), data = data_labels1,
            color = "grey50",
            size = 7,
            face = "bold")+
  geom_text(x = -Inf, y = Inf,
            hjust = -0.1, vjust = 3.32,
            aes(label = label, group = NULL), data = data_labels2,
            color = "red",
            size = 7,
            face = "bold")
figure1d
```

![](Fig.-1-Comparison-between-an-Orbitrap-Mass-analyzer-and-Linear-Ion-Trap-Mass-Analyzer-on-low-input-proteomics_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

# Join A,B,C and D

Plot 1A, 1B, 1C and 1D are combined using patchwork to create figure 1.

``` r
temp <- figure1a +  figure1b + plot_layout(widths = c(1, 2))
temp2 <- plot_spacer() + figure1c + plot_spacer() + plot_layout(widths = c(1, 10, 1))
temp3 <- plot_spacer() + figure1d + plot_spacer() + plot_layout(widths = c(1, 10, 1))
figure1 <- temp/temp2/temp3 + plot_layout(heights = c(0.8, 1, 1))
figure1
```

![](Fig.-1-Comparison-between-an-Orbitrap-Mass-analyzer-and-Linear-Ion-Trap-Mass-Analyzer-on-low-input-proteomics_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

Save it as PDF.

``` r
cairo_pdf(filename = "fig1.pdf", width = 20.83, height = 27.08)
figure1
dev.off()
```

    ## png 
    ##   2

# Session information

    ## R version 4.2.0 (2022-04-22 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8   
    ## [3] LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    ## [5] LC_TIME=Danish_Denmark.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggpointdensity_0.1.0 patchwork_1.1.1      ggpubr_0.4.0        
    ##  [4] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.9         
    ##  [7] purrr_0.3.4          readr_2.1.2          tidyr_1.2.0         
    ## [10] tibble_3.1.7         ggplot2_3.3.6        tidyverse_1.3.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lattice_0.20-45  lubridate_1.8.0  assertthat_0.2.1 digest_0.6.29   
    ##  [5] utf8_1.2.2       R6_2.5.1         cellranger_1.1.0 backports_1.4.1 
    ##  [9] reprex_2.0.1     evaluate_0.15    highr_0.9        httr_1.4.3      
    ## [13] pillar_1.7.0     rlang_1.0.2      readxl_1.4.0     rstudioapi_0.13 
    ## [17] car_3.0-13       Matrix_1.4-1     rmarkdown_2.14   splines_4.2.0   
    ## [21] labeling_0.4.2   munsell_0.5.0    broom_0.8.0      compiler_4.2.0  
    ## [25] modelr_0.1.8     xfun_0.31        pkgconfig_2.0.3  mgcv_1.8-40     
    ## [29] htmltools_0.5.2  tidyselect_1.1.2 fansi_1.0.3      crayon_1.5.1    
    ## [33] tzdb_0.3.0       dbplyr_2.1.1     withr_2.5.0      grid_4.2.0      
    ## [37] nlme_3.1-157     jsonlite_1.8.0   gtable_0.3.0     lifecycle_1.0.1 
    ## [41] DBI_1.1.2        magrittr_2.0.3   scales_1.2.0     cli_3.3.0       
    ## [45] stringi_1.7.6    carData_3.0-5    farver_2.1.0     ggsignif_0.6.3  
    ## [49] fs_1.5.2         xml2_1.3.3       ellipsis_0.3.2   generics_0.1.2  
    ## [53] vctrs_0.4.1      tools_4.2.0      glue_1.6.2       hms_1.1.1       
    ## [57] abind_1.4-5      fastmap_1.1.0    yaml_2.3.5       colorspace_2.0-3
    ## [61] rstatix_0.7.0    rvest_1.0.2      knitr_1.39       haven_2.5.0
