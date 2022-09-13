Fig. 6 Comparison of Injection time on linear ion trap mass analyzer.
================
Samuel Grégoire & Lukas R. Woltereck

Fig. 6 Comparison of Injection time on linear ion trap mass analyzer.

Comparison of the number of identified peptides on DIA-LIT-based method
on Turbo, Rapid, and Normal scanning mode for WhisperTM 100 40 SPD from
1, 5, and 100 ng of tryptic HeLa lysate with different ITs at fixed 40
isolation windows. For Turbo scanning mode auto IT (16 ms), half of its
auto-IT (8 ms), and the auto IT of Rapid (23 ms) and Normal (38 ms) were
applied to determine the compromise between ions filling time and
scanning speed of LIT on this mode. For Rapid scanning mode, the ITs of
23 ms and 38 ms suffice. On Normal scanning mode, only its auto-IT was
applied for this comparison.

# Load packages

``` r
library("tidyverse")
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library("patchwork")
#set bw theme as default theme
theme_set(theme_bw())
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

## 8ms

``` r
#Turbo
LITDIA_Turbo_1_8ms <- read.csv2(paste0(lit_path, "20220420_LIT_Turbo_DIA_1ng_40SPD_whisper100_1CV_8ms.csv"))
LITDIA_Turbo_5_8ms <- read.csv2(paste0(lit_path, "20220420_LIT_Turbo_DIA_5ng_40SPD_whisper100_1CV_8ms.csv"))
LITDIA_Turbo_100_8ms <- read.csv2(paste0(lit_path, "20220420_LIT_Turbo_DIA_100ng_40SPD_whisper100_1CV_8ms.csv"))
```

## 16ms

``` r
#Turbo
LITDIA_Turbo_1_16ms <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_5_16ms <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_100_16ms <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

## 23ms

``` r
#Turbo
LITDIA_Turbo_1_23ms <- read.csv2(paste0(lit_path, "20220524_LIT_Turbo_DIA_1ng_40SPD_whisper100_1CV_23ms.csv"))
LITDIA_Turbo_5_23ms <- read.csv2(paste0(lit_path, "20220524_LIT_Turbo_DIA_5ng_40SPD_whisper100_1CV_23ms.csv"))
LITDIA_Turbo_100_23ms <- read.csv2(paste0(lit_path, "20220525_LIT_Turbo_DIA_100ng_40SPD_whisper100_1CV_23ms.csv"))

#Rapid
LITDIA_Rapid_1_23ms <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_5_23ms <- read.csv2(paste0(lit_path, "20220502_LIT_Rapid_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_100_23ms <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

## 38ms

``` r
#Turbo
LITDIA_Turbo_1_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Turbo_DIA_1ng_40SPD_whisper100_1CV_38ms.csv"))
LITDIA_Turbo_5_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Turbo_DIA_5ng_40SPD_whisper100_1CV_38ms.csv"))
LITDIA_Turbo_100_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Turbo_DIA_100ng_40SPD_whisper100_1CV_38ms.csv"))

#Rapid
LITDIA_Rapid_1_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_38ms.csv"))
LITDIA_Rapid_5_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Rapid_DIA_5ng_40SPD_whisper100_1CV_38ms_40SPD_File1.csv"))
LITDIA_Rapid_100_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Rapid_DIA_100ng_40SPD_whisper100_1CV_38ms.csv"))

#Normal
LITDIA_Normal_1_38ms <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_5_38ms <- read.csv2(paste0(lit_path, "20220502_LIT_Normal_DIA_5ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_100_38ms <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

# Clean data tables

Imported tables contain a lot of information, which are not useful for
our analysis. The function `filter_pep()` only keeps the file names
(replicates), quantitative values for each file name, and the peptide
sequence. It is also calculating the mean between the four replicates.
Spectronaut does some imputation by default , but we want to keep
missing values. Imputed values are replaced by NA. Peptides missing in
more than one replicate are highlighted as *TRUE*.

``` r
filter_pep <- function(x) {
  x$FG.Quantity[as.logical(x$EG.IsImputed)] <- NA
  x_clean_with_mean <- x %>% 
  group_by(R.FileName, PEP.StrippedSequence) %>% 
  summarise(mean_quantity = mean(FG.Quantity, na.rm = TRUE)) %>%  
  pivot_wider(names_from = R.FileName, 
              values_from = mean_quantity) %>%
    suppressMessages()
  x_filter_highlight <- x_clean_with_mean %>% 
  mutate(na_nbr = rowSums(is.na(x_clean_with_mean[-1]))) %>%
  mutate(mean_quantity = 
           rowMeans(x_clean_with_mean[2:ncol(x_clean_with_mean)], 
                    na.rm = TRUE)) %>%
  mutate(filtered = !na_nbr <= 1) %>% 
  select(-na_nbr) %>% 
  suppressMessages("`summarise()` has grouped output by 'R.FileName'. You can override using the `.groups` argument.")
  message("
Highlighted ", sum(x_filter_highlight$filtered) , " peptide(s) found in less than ", ncol(x_clean_with_mean) -2, " replicates")
  
  return(x_filter_highlight)
}
```

The `filter_pep()` function is used for every data-set.

## 8ms

``` r
#Turbo
LITDIA_Turbo_1_8ms <-  filter_pep(LITDIA_Turbo_1_8ms)
```

    ## 
    ## Highlighted 13 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_5_8ms <-  filter_pep(LITDIA_Turbo_5_8ms)
```

    ## 
    ## Highlighted 13 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_8ms <-  filter_pep(LITDIA_Turbo_100_8ms)
```

    ## 
    ## Highlighted 219 peptide(s) found in less than 3 replicates

## 16ms

``` r
#Turbo
LITDIA_Turbo_1_16ms <-  filter_pep(LITDIA_Turbo_1_16ms)
```

    ## 
    ## Highlighted 14 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_5_16ms <-  filter_pep(LITDIA_Turbo_5_16ms)
```

    ## 
    ## Highlighted 27 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_16ms <-  filter_pep(LITDIA_Turbo_100_16ms)
```

    ## 
    ## Highlighted 324 peptide(s) found in less than 3 replicates

## 23ms

``` r
#Turbo
LITDIA_Turbo_1_23ms <-  filter_pep(LITDIA_Turbo_1_23ms)
```

    ## 
    ## Highlighted 44 peptide(s) found in less than 2 replicates

``` r
LITDIA_Turbo_5_23ms <-  filter_pep(LITDIA_Turbo_5_23ms)
```

    ## 
    ## Highlighted 22 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_23ms <-  filter_pep(LITDIA_Turbo_100_23ms)
```

    ## 
    ## Highlighted 258 peptide(s) found in less than 3 replicates

``` r
#Rapid
LITDIA_Rapid_1_23ms <-  filter_pep(LITDIA_Rapid_1_23ms)
```

    ## 
    ## Highlighted 12 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_5_23ms <-  filter_pep(LITDIA_Rapid_5_23ms)
```

    ## 
    ## Highlighted 13 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_100_23ms <-  filter_pep(LITDIA_Rapid_100_23ms)
```

    ## 
    ## Highlighted 282 peptide(s) found in less than 3 replicates

## 38ms

``` r
#Turbo
LITDIA_Turbo_1_38ms <-  filter_pep(LITDIA_Turbo_1_38ms)
```

    ## 
    ## Highlighted 17 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_5_38ms <-  filter_pep(LITDIA_Turbo_5_38ms)
```

    ## 
    ## Highlighted 20 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_38ms <-  filter_pep(LITDIA_Turbo_100_38ms)
```

    ## 
    ## Highlighted 346 peptide(s) found in less than 3 replicates

``` r
#Rapid
LITDIA_Rapid_1_38ms <-  filter_pep(LITDIA_Rapid_1_38ms)
```

    ## 
    ## Highlighted 19 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_5_38ms <-  filter_pep(LITDIA_Rapid_5_38ms)
```

    ## 
    ## Highlighted 31 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_100_38ms <-  filter_pep(LITDIA_Rapid_100_38ms)
```

    ## 
    ## Highlighted 287 peptide(s) found in less than 3 replicates

``` r
#Normal
LITDIA_Normal_1_38ms <-  filter_pep(LITDIA_Normal_1_38ms)
```

    ## 
    ## Highlighted 41 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_5_38ms <-  filter_pep(LITDIA_Normal_5_38ms)
```

    ## 
    ## Highlighted 30 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_100_38ms <-  filter_pep(LITDIA_Normal_100_38ms)
```

    ## 
    ## Highlighted 316 peptide(s) found in less than 3 replicates

# Count number of identified peptides

Empty vectors to store the number of identified peptides in each
replicates for each data-set are created.

``` r
id_LITDIA_Turbo_1_8ms <- id_LITDIA_Turbo_5_8ms <- 
  id_LITDIA_Turbo_100_8ms <- id_LITDIA_Turbo_1_16ms <- 
  id_LITDIA_Turbo_5_16ms <- id_LITDIA_Turbo_100_16ms <- 
  id_LITDIA_Turbo_1_23ms <- id_LITDIA_Turbo_5_23ms <-
  id_LITDIA_Turbo_100_23ms <- id_LITDIA_Rapid_1_23ms <- 
  id_LITDIA_Rapid_5_23ms <- id_LITDIA_Rapid_100_23ms <-
  id_LITDIA_Turbo_1_38ms <- id_LITDIA_Turbo_5_38ms <-
  id_LITDIA_Turbo_100_38ms <- id_LITDIA_Rapid_1_38ms <-
  id_LITDIA_Rapid_5_38ms <- id_LITDIA_Rapid_100_38ms <-
  id_LITDIA_Normal_1_38ms <- id_LITDIA_Normal_5_38ms <-
  id_LITDIA_Normal_100_38ms <- rep(NA, 4)
```

The identified peptides, that are not *NA*, are summed for each
replicate and replace the empty vectors.

## 8ms

``` r
#Turbo
for(i in 2:5){
  id_LITDIA_Turbo_1_8ms[i-1] <- sum(!is.na(LITDIA_Turbo_1_8ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_5_8ms[i-1] <- sum(!is.na(LITDIA_Turbo_5_8ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_100_8ms[i-1] <- sum(!is.na(LITDIA_Turbo_100_8ms[, i]))}
```

## 16ms

``` r
#Turbo
for(i in 2:5){
  id_LITDIA_Turbo_1_16ms[i-1] <- sum(!is.na(LITDIA_Turbo_1_16ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_5_16ms[i-1] <- sum(!is.na(LITDIA_Turbo_5_16ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_100_16ms[i-1] <- sum(!is.na(LITDIA_Turbo_100_16ms[, i]))}
```

## 23ms

``` r
#Turbo
for(i in 2:5){
  id_LITDIA_Turbo_1_23ms[i-1] <- sum(!is.na(LITDIA_Turbo_1_23ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_5_23ms[i-1] <- sum(!is.na(LITDIA_Turbo_5_23ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_100_23ms[i-1] <- sum(!is.na(LITDIA_Turbo_100_23ms[, i]))}

#Rapid
for(i in 2:5){
  id_LITDIA_Rapid_1_23ms[i-1] <- sum(!is.na(LITDIA_Rapid_1_23ms[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_5_23ms[i-1] <- sum(!is.na(LITDIA_Rapid_5_23ms[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_100_23ms[i-1] <- sum(!is.na(LITDIA_Rapid_100_23ms[, i]))}
```

## 38ms

``` r
#Turbo
for(i in 2:5){
  id_LITDIA_Turbo_1_38ms[i-1] <- sum(!is.na(LITDIA_Turbo_1_38ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_5_38ms[i-1] <- sum(!is.na(LITDIA_Turbo_5_38ms[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_100_38ms[i-1] <- sum(!is.na(LITDIA_Turbo_100_38ms[, i]))}

#Rapid
for(i in 2:5){
  id_LITDIA_Rapid_1_38ms[i-1] <- sum(!is.na(LITDIA_Rapid_1_38ms[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_5_38ms[i-1] <- sum(!is.na(LITDIA_Rapid_5_38ms[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_100_38ms[i-1] <- sum(!is.na(LITDIA_Rapid_100_38ms[, i]))}

#Normal
for(i in 2:5){
  id_LITDIA_Normal_1_38ms[i-1] <- sum(!is.na(LITDIA_Normal_1_38ms[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_5_38ms[i-1] <- sum(!is.na(LITDIA_Normal_5_38ms[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_100_38ms[i-1] <- sum(!is.na(LITDIA_Normal_100_38ms[, i]))}
```

# Compute CV

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

Information about the scanning mode, the injection time and the input
are included afterwards.

## 8 ms

``` r
#Turbo
LITDIA_Turbo_1_8ms_cv <- compute_CV(LITDIA_Turbo_1_8ms)
LITDIA_Turbo_1_8ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_1_8ms_cv$injection_time<- "8 ms"
LITDIA_Turbo_1_8ms_cv$input<- 1

LITDIA_Turbo_5_8ms_cv <- compute_CV(LITDIA_Turbo_5_8ms)
LITDIA_Turbo_5_8ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_5_8ms_cv$injection_time<- "8 ms"
LITDIA_Turbo_5_8ms_cv$input<- 5

LITDIA_Turbo_100_8ms_cv <- compute_CV(LITDIA_Turbo_100_8ms)
LITDIA_Turbo_100_8ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_100_8ms_cv$injection_time<- "8 ms"
LITDIA_Turbo_100_8ms_cv$input<- 100
```

## 16 ms

``` r
#Turbo
LITDIA_Turbo_1_16ms_cv <- compute_CV(LITDIA_Turbo_1_16ms)
LITDIA_Turbo_1_16ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_1_16ms_cv$injection_time<- "16 ms"
LITDIA_Turbo_1_16ms_cv$input<- 1

LITDIA_Turbo_5_16ms_cv <- compute_CV(LITDIA_Turbo_5_16ms)
LITDIA_Turbo_5_16ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_5_16ms_cv$injection_time<- "16 ms"
LITDIA_Turbo_5_16ms_cv$input<- 5

LITDIA_Turbo_100_16ms_cv <- compute_CV(LITDIA_Turbo_100_16ms)
LITDIA_Turbo_100_16ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_100_16ms_cv$injection_time<- "16 ms"
LITDIA_Turbo_100_16ms_cv$input<- 100
```

## 23 ms

``` r
#Turbo
LITDIA_Turbo_1_23ms_cv <- compute_CV(LITDIA_Turbo_1_23ms)
LITDIA_Turbo_1_23ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_1_23ms_cv$injection_time<- "23 ms"
LITDIA_Turbo_1_23ms_cv$input<- 1

LITDIA_Turbo_5_23ms_cv <- compute_CV(LITDIA_Turbo_5_23ms)
LITDIA_Turbo_5_23ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_5_23ms_cv$injection_time<- "23 ms"
LITDIA_Turbo_5_23ms_cv$input<- 5

LITDIA_Turbo_100_23ms_cv <- compute_CV(LITDIA_Turbo_100_23ms)
LITDIA_Turbo_100_23ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_100_23ms_cv$injection_time<- "23 ms"
LITDIA_Turbo_100_23ms_cv$input<- 100

#Rapid
LITDIA_Rapid_1_23ms_cv <- compute_CV(LITDIA_Rapid_1_23ms)
LITDIA_Rapid_1_23ms_cv$scanning_mode <- "Rapid"
LITDIA_Rapid_1_23ms_cv$injection_time<- "23 ms"
LITDIA_Rapid_1_23ms_cv$input<- 1

LITDIA_Rapid_5_23ms_cv <- compute_CV(LITDIA_Rapid_5_23ms)
LITDIA_Rapid_5_23ms_cv$scanning_mode <- "Rapid"
LITDIA_Rapid_5_23ms_cv$injection_time<- "23 ms"
LITDIA_Rapid_5_23ms_cv$input<- 5

LITDIA_Rapid_100_23ms_cv <- compute_CV(LITDIA_Rapid_100_23ms)
LITDIA_Rapid_100_23ms_cv$scanning_mode <- "Rapid"
LITDIA_Rapid_100_23ms_cv$injection_time<- "23 ms"
LITDIA_Rapid_100_23ms_cv$input<- 100
```

## 38 ms

``` r
#Turbo
LITDIA_Turbo_1_38ms_cv <- compute_CV(LITDIA_Turbo_1_38ms)
LITDIA_Turbo_1_38ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_1_38ms_cv$injection_time<- "38 ms"
LITDIA_Turbo_1_38ms_cv$input<- 1

LITDIA_Turbo_5_38ms_cv <- compute_CV(LITDIA_Turbo_5_38ms)
LITDIA_Turbo_5_38ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_5_38ms_cv$injection_time<- "38 ms"
LITDIA_Turbo_5_38ms_cv$input<- 5

LITDIA_Turbo_100_38ms_cv <- compute_CV(LITDIA_Turbo_100_38ms)
LITDIA_Turbo_100_38ms_cv$scanning_mode <- "Turbo"
LITDIA_Turbo_100_38ms_cv$injection_time<- "38 ms"
LITDIA_Turbo_100_38ms_cv$input<- 100

#Rapid
LITDIA_Rapid_1_38ms_cv <- compute_CV(LITDIA_Rapid_1_38ms)
LITDIA_Rapid_1_38ms_cv$scanning_mode <- "Rapid"
LITDIA_Rapid_1_38ms_cv$injection_time<- "38 ms"
LITDIA_Rapid_1_38ms_cv$input<- 1

LITDIA_Rapid_5_38ms_cv <- compute_CV(LITDIA_Rapid_5_38ms)
LITDIA_Rapid_5_38ms_cv$scanning_mode <- "Rapid"
LITDIA_Rapid_5_38ms_cv$injection_time<- "38 ms"
LITDIA_Rapid_5_38ms_cv$input<- 5

LITDIA_Rapid_100_38ms_cv <- compute_CV(LITDIA_Rapid_100_38ms)
LITDIA_Rapid_100_38ms_cv$scanning_mode <- "Rapid"
LITDIA_Rapid_100_38ms_cv$injection_time<- "38 ms"
LITDIA_Rapid_100_38ms_cv$input<- 100

#Normal
LITDIA_Normal_1_38ms_cv <- compute_CV(LITDIA_Normal_1_38ms)
LITDIA_Normal_1_38ms_cv$scanning_mode <- "Normal"
LITDIA_Normal_1_38ms_cv$injection_time<- "38 ms"
LITDIA_Normal_1_38ms_cv$input<- 1

LITDIA_Normal_5_38ms_cv <- compute_CV(LITDIA_Normal_5_38ms)
LITDIA_Normal_5_38ms_cv$scanning_mode <- "Normal"
LITDIA_Normal_5_38ms_cv$injection_time<- "38 ms"
LITDIA_Normal_5_38ms_cv$input<- 5

LITDIA_Normal_100_38ms_cv <- compute_CV(LITDIA_Normal_100_38ms)
LITDIA_Normal_100_38ms_cv$scanning_mode <- "Normal"
LITDIA_Normal_100_38ms_cv$injection_time<- "38 ms"
LITDIA_Normal_100_38ms_cv$input<- 100
```

# Join

Data frames containing CV are combined with the function `rbind()` and
arranged with the function `factor()`.

``` r
id_figure6 <- 
  rbind(LITDIA_Turbo_1_8ms_cv, LITDIA_Turbo_5_8ms_cv,
        LITDIA_Turbo_100_8ms_cv, LITDIA_Turbo_1_16ms_cv,
        LITDIA_Turbo_5_16ms_cv, LITDIA_Turbo_100_16ms_cv, 
        LITDIA_Turbo_1_23ms_cv, LITDIA_Turbo_5_23ms_cv,
        LITDIA_Turbo_100_23ms_cv, LITDIA_Rapid_1_23ms_cv,
        LITDIA_Rapid_5_23ms_cv, LITDIA_Rapid_100_23ms_cv, 
        LITDIA_Turbo_1_38ms_cv, LITDIA_Turbo_5_38ms_cv,
        LITDIA_Turbo_100_38ms_cv, LITDIA_Rapid_1_38ms_cv,
        LITDIA_Rapid_5_38ms_cv, LITDIA_Rapid_100_38ms_cv, 
        LITDIA_Normal_1_38ms_cv, LITDIA_Normal_5_38ms_cv,
        LITDIA_Normal_100_38ms_cv)

id_figure6$scanning_mode %>% table()
```

    ## .
    ## Normal  Rapid  Turbo 
    ##  24643  45842  78492

``` r
id_figure6$injection_time %>% table()
```

    ## .
    ## 16 ms 23 ms 38 ms  8 ms 
    ## 19681 43604 71541 14151

``` r
id_figure6$input %>% table()
```

    ## .
    ##     1     5   100 
    ## 29025 53013 66939

``` r
id_figure6$scanning_mode <- factor(id_figure6$scanning_mode, 
                          levels = unique(id_figure6$scanning_mode))

id_figure6$input <- factor(id_figure6$input, 
                          levels = unique(id_figure6$input))

id_figure6$injection_time <- factor(id_figure6$injection_time, 
                          levels = unique(id_figure6$injection_time))
```

A new column called “CV_filter” is added to identify peptides with a CV
over 0.15 (1), CV between 0.1 and 0.15 (2) and CV under 0.1 (3).

``` r
id_figure6$CV_filter <- NA
id_figure6$CV_filter[id_figure6$CV < 0.1] <- "3"
id_figure6$CV_filter[id_figure6$CV >= 0.1 & id_figure6$CV <= 0.15] <- "2"
id_figure6$CV_filter[id_figure6$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out and
the number of identified peptides for each replicate and each CV_filter
is calculated by summing the number of non-missing peptides. The mean
and the standard deviation of the replicates are also calculated for
each method.

``` r
id_figure6 <- id_figure6 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(scanning_mode, injection_time, input, CV_filter) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
```

    ## `summarise()` has grouped output by 'scanning_mode', 'injection_time', 'input'.
    ## You can override using the `.groups` argument.

``` r
id_figure6$id_S4[id_figure6$id_S4 == 0] <- NA 

id_figure6$mean_id <- rowMeans(id_figure6[5:8], na.rm = TRUE)
id_figure6$sd_id <- apply(id_figure6[5:8], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

# Plot

For the figure the data from Turbo, Rapid and Normal were combined with
patchwork.

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

## Turbo

``` r
id_figure6$y_errorbar <- id_figure6$mean_id
id_figure6$y_errorbar[id_figure6$CV_filter == "1"] <-
  id_figure6$y_errorbar[id_figure6$CV_filter == "1"] +
  id_figure6$y_errorbar[id_figure6$CV_filter == "2"]+
  id_figure6$y_errorbar[id_figure6$CV_filter == "3"]
id_figure6$y_errorbar[id_figure6$CV_filter == "2"] <-
  id_figure6$y_errorbar[id_figure6$CV_filter == "2"]+
  id_figure6$y_errorbar[id_figure6$CV_filter == "3"]

figure6_turbo <- id_figure6 %>% filter(scanning_mode == "Turbo") %>% 
  ggplot(aes(x = input, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
 facet_grid(scanning_mode~injection_time)+
  labs(x = "",
       y = "Identified peptides",
       fill = "CV",
       title = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.09),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 27),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.title.y = element_text(vjust = 2.2),
        legend.title = element_text(size = 27),
        axis.title.x = element_text(vjust = 0.2),
        legend.text = element_text(size = 27),
        legend.key.size = unit(0.9, "cm"),
        legend.position = "top")+
  scale_fill_manual(values = c("grey82", "#ffaa8b", "#ff713e"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 16000))+
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure6_turbo
```

![](Figure6_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

## Rapid

``` r
id_figure6$y_errorbar <- id_figure6$mean_id
id_figure6$y_errorbar[id_figure6$CV_filter == "1"] <-
  id_figure6$y_errorbar[id_figure6$CV_filter == "1"] +
  id_figure6$y_errorbar[id_figure6$CV_filter == "2"]+
  id_figure6$y_errorbar[id_figure6$CV_filter == "3"]
id_figure6$y_errorbar[id_figure6$CV_filter == "2"] <-
  id_figure6$y_errorbar[id_figure6$CV_filter == "2"]+
  id_figure6$y_errorbar[id_figure6$CV_filter == "3"]

figure6_rapid <- id_figure6 %>% filter(scanning_mode == "Rapid") %>% 
  ggplot(aes(x = input, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
 facet_grid(scanning_mode~injection_time)+
  labs(x = "Input (ng)",
       y = "Identified peptides",
       fill = "CV",
       title = "") +
  theme_bw() +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.09),
    axis.text = element_text(size = 24),
        axis.title = element_text(size = 27),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.title.y = element_text(vjust = 2.2),
        legend.title = element_text(size = 27),
        axis.title.x = element_text(vjust = 0.2),
        legend.text = element_text(size = 27),
        legend.key.size = unit(0.9, "cm"))+
scale_fill_manual(values = c("grey82", "#ffaa8b", "#ff713e"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 16000))+
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure6_rapid
```

![](Figure6_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

## Normal

``` r
id_figure6$y_errorbar <- id_figure6$mean_id
id_figure6$y_errorbar[id_figure6$CV_filter == "1"] <-
  id_figure6$y_errorbar[id_figure6$CV_filter == "1"] +
  id_figure6$y_errorbar[id_figure6$CV_filter == "2"]+
  id_figure6$y_errorbar[id_figure6$CV_filter == "3"]
id_figure6$y_errorbar[id_figure6$CV_filter == "2"] <-
  id_figure6$y_errorbar[id_figure6$CV_filter == "2"]+
  id_figure6$y_errorbar[id_figure6$CV_filter == "3"]

figure6_normal <- id_figure6 %>% filter(scanning_mode == "Normal") %>% 
  ggplot(aes(x = input, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
 facet_grid(scanning_mode~injection_time)+
  labs(x = "Input (ng)",
       y = "",
       fill = "CV",
       title = "") +
  theme_bw() +
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.09),
    axis.text = element_text(size = 24),
        axis.title = element_text(size = 27),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.title.y = element_text(vjust = 2.2),
        legend.title = element_text(size = 27),
        axis.title.x = element_text(vjust = 0.2),
        legend.text = element_text(size = 27),
        legend.key.size = unit(0.9, "cm"))+
scale_fill_manual(values = c("grey82", "#ffaa8b", "#ff713e"),
                  labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 16000))+
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),  
                  width = 0.2,
                    size = 1)
figure6_normal
```

![](Figure6_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

# Join Turbo, Rapid, Normal

Figure6_turbo, Figure6_rapid and Figure6_normal are combined using
patchwork to create figure 4.

``` r
figure6 <- figure6_turbo/((figure6_rapid+figure6_normal)+plot_layout(widths = c(2,1)))
figure6
```

![](Figure6_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
#save as pdf
#cairo_pdf(filename = "fig6.pdf", width = 20.83, height = 15.625)
#figure6
#dev.off()
```

# Session information

    ## R version 4.2.1 (2022-06-23 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=French_Belgium.utf8  LC_CTYPE=French_Belgium.utf8   
    ## [3] LC_MONETARY=French_Belgium.utf8 LC_NUMERIC=C                   
    ## [5] LC_TIME=French_Belgium.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.1.1 forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9    
    ##  [5] purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.7   
    ##  [9] ggplot2_3.3.6   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.31        haven_2.5.0      colorspace_2.0-3
    ##  [5] vctrs_0.4.1      generics_0.1.2   htmltools_0.5.2  yaml_2.3.5      
    ##  [9] utf8_1.2.2       rlang_1.0.2      pillar_1.7.0     glue_1.6.2      
    ## [13] withr_2.5.0      DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8    
    ## [17] readxl_1.4.0     lifecycle_1.0.1  munsell_0.5.0    gtable_0.3.0    
    ## [21] cellranger_1.1.0 rvest_1.0.2      evaluate_0.15    labeling_0.4.2  
    ## [25] knitr_1.39       tzdb_0.3.0       fastmap_1.1.0    fansi_1.0.3     
    ## [29] highr_0.9        broom_0.8.0      backports_1.4.1  scales_1.2.0    
    ## [33] jsonlite_1.8.0   farver_2.1.0     fs_1.5.2         hms_1.1.1       
    ## [37] digest_0.6.29    stringi_1.7.6    grid_4.2.1       cli_3.3.0       
    ## [41] tools_4.2.1      magrittr_2.0.3   crayon_1.5.1     pkgconfig_2.0.3 
    ## [45] ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0 
    ## [49] assertthat_0.2.1 rmarkdown_2.14   httr_1.4.3       rstudioapi_0.13 
    ## [53] R6_2.5.1         compiler_4.2.1
