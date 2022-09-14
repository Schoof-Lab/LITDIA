Fig. 8 Comparison of low-input protein identifications between
SpectronautTM version 15 and 16
================
Lukas R. Woltereck & Samuel Grégoire

Fig. 8 Comparison of low-input protein identifications between
SpectronautTM version 15 and 16

Comparison of the number of identified peptides in serial dilution
(1,10, and 100 ng) of HeLa tryptic digested analyzed by DIA-LIT-based
methods with Normal scanning mode on 40 SPD between Spectronaut version
15 and 16.

# Load packages

``` r
library("tidyverse")
library("patchwork")

#set bw theme as default theme
theme_set(theme_bw())
```

# Load data

The path to every folder used is stored to keep everything compact. LIT
data is stored in the folder located in lit\_path. To reproduce analysis
the path must be changed to corresponding path on the local computer.The
path to every folder used is stored to keep everything compact. LIT data
is stored in the folder located in lit\_path. To reproduce analysis the
path must be changed to corresponding path on the local computer.

``` r
#lit_path <- "C:/Data/LIT/"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in RStudio with the function
`read.csv2()`.Data tables were exported from Spectronaut (version
15.5.211111.50606 or 16.0.220606.53000) in csv files and need to be
loaded in RStudio with the function `read.csv2()`.

## Spectronaut V.15

``` r
#Turbo
LITDIA_Turbo_1_SN15 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_10_SN15 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Turbo_100_SN15 <- read.csv2(paste0(lit_path, "20220502_LIT_Turbo_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
#Rapid
LITDIA_Rapid_1_SN15 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_10_SN15 <- read.csv2(paste0(lit_path, "20220502_LIT_Rapid_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Rapid_100_SN15 <- read.csv2(paste0(lit_path, "20220421_LIT_Rapid_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
#Normal
LITDIA_Normal_1_SN15 <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_10_SN15 <- read.csv2(paste0(lit_path, "20220502_LIT_Normal_DIA_10ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_100_SN15 <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
```

## Spectronaut V.16

``` r
#Turbo
LITDIA_Turbo_1_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Turbo_1CV_1ng_auto_40w.csv"))
LITDIA_Turbo_10_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Turbo_1CV_10ng_auto_40w.csv"))
LITDIA_Turbo_100_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Turbo_1CV_100ng_auto_40w.csv"))
#Rapid
LITDIA_Rapid_1_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Rapid_1CV_1ng_auto_40w.csv"))
LITDIA_Rapid_10_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Rapid_1CV_10ng_auto_40w.csv"))
LITDIA_Rapid_100_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Rapid_1CV_100ng_auto_40w.csv"))
#Normal
LITDIA_Normal_1_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Normal_1CV_1ng_auto_40w.csv"))
LITDIA_Normal_10_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Normal_1CV_10ng_auto_40w.csv"))
LITDIA_Normal_100_SN16 <- read.csv2(paste0(lit_path, "20220615_SN16_40SPD_LIT_Normal_1CV_100ng_auto_40w.csv"))
```

# Clean data tables

Imported tables contain a lot of information, which are not useful for
our analysis. The function `filter_pep()` only keeps the file names
(replicates), quantitative values for each file name, and the peptide
sequence. It is also calculating the mean between the four replicates.
Spectronaut does imputation by default, but we want to keep missing
values. Imputed values are replaced by *NA*. Peptides missing in more
than one replicate are highlighted as *TRUE*.

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

## Spectronaut V.15

``` r
LITDIA_Turbo_1_SN15 <-  filter_pep(LITDIA_Turbo_1_SN15)
```

    ## 
    ## Highlighted 14 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_10_SN15 <-  filter_pep(LITDIA_Turbo_10_SN15)
```

    ## 
    ## Highlighted 14 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_SN15 <-  filter_pep(LITDIA_Turbo_100_SN15)
```

    ## 
    ## Highlighted 324 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_1_SN15 <-  filter_pep(LITDIA_Rapid_1_SN15)
```

    ## 
    ## Highlighted 12 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_10_SN15 <-  filter_pep(LITDIA_Rapid_10_SN15)
```

    ## 
    ## Highlighted 21 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_100_SN15 <-  filter_pep(LITDIA_Rapid_100_SN15)
```

    ## 
    ## Highlighted 282 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_1_SN15 <-  filter_pep(LITDIA_Normal_1_SN15)
```

    ## 
    ## Highlighted 41 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_10_SN15 <-  filter_pep(LITDIA_Normal_10_SN15)
```

    ## 
    ## Highlighted 28 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_100_SN15 <-  filter_pep(LITDIA_Normal_100_SN15)
```

    ## 
    ## Highlighted 316 peptide(s) found in less than 3 replicates

## Spectronaut V.16

``` r
LITDIA_Turbo_1_SN16 <-  filter_pep(LITDIA_Turbo_1_SN16)
```

    ## 
    ## Highlighted 10 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_10_SN16 <-  filter_pep(LITDIA_Turbo_10_SN16)
```

    ## 
    ## Highlighted 20 peptide(s) found in less than 3 replicates

``` r
LITDIA_Turbo_100_SN16 <-  filter_pep(LITDIA_Turbo_100_SN16)
```

    ## 
    ## Highlighted 267 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_1_SN16 <-  filter_pep(LITDIA_Rapid_1_SN16)
```

    ## 
    ## Highlighted 8 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_10_SN16 <-  filter_pep(LITDIA_Rapid_10_SN16)
```

    ## 
    ## Highlighted 27 peptide(s) found in less than 3 replicates

``` r
LITDIA_Rapid_100_SN16 <-  filter_pep(LITDIA_Rapid_100_SN16)
```

    ## 
    ## Highlighted 300 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_1_SN16 <-  filter_pep(LITDIA_Normal_1_SN16)
```

    ## 
    ## Highlighted 16 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_10_SN16 <-  filter_pep(LITDIA_Normal_10_SN16)
```

    ## 
    ## Highlighted 29 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_100_SN16 <-  filter_pep(LITDIA_Normal_100_SN16)
```

    ## 
    ## Highlighted 264 peptide(s) found in less than 3 replicates

# Count number of identified peptides

Empty vectors are created to store the number of identified peptides in
each replicates for each data-set.

``` r
id_LITDIA_Turbo_1_SN15 <- id_LITDIA_Turbo_10_SN15 <- id_LITDIA_Turbo_100_SN15 <-
id_LITDIA_Rapid_1_SN15 <- id_LITDIA_Rapid_10_SN15 <- id_LITDIA_Rapid_100_SN15 <-
id_LITDIA_Normal_1_SN15 <- id_LITDIA_Normal_10_SN15 <- id_LITDIA_Normal_100_SN15 <-
id_LITDIA_Turbo_1_SN16 <- id_LITDIA_Turbo_10_SN16 <- id_LITDIA_Turbo_100_SN16 <-
id_LITDIA_Rapid_1_SN16 <- id_LITDIA_Rapid_10_SN16 <- id_LITDIA_Rapid_100_SN16 <-
id_LITDIA_Normal_1_SN16 <- id_LITDIA_Normal_10_SN16 <- id_LITDIA_Normal_100_SN16 <- rep(NA, 4)
```

The numbers of identified peptides are counted by summing non-missing
values for each replicate and stored in the vectors created above.

## Spectronaut V.15

``` r
#Turbo
for(i in 2:5){
  id_LITDIA_Turbo_1_SN15[i-1] <- sum(!is.na(LITDIA_Turbo_1_SN15[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_10_SN15[i-1] <- sum(!is.na(LITDIA_Turbo_10_SN15[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_100_SN15[i-1] <- sum(!is.na(LITDIA_Turbo_100_SN15[, i]))}

#Rapid
for(i in 2:5){
  id_LITDIA_Rapid_1_SN15[i-1] <- sum(!is.na(LITDIA_Rapid_1_SN15[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_10_SN15[i-1] <- sum(!is.na(LITDIA_Rapid_10_SN15[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_100_SN15[i-1] <- sum(!is.na(LITDIA_Rapid_100_SN15[, i]))}

#Normal
for(i in 2:5){
  id_LITDIA_Normal_1_SN15[i-1] <- sum(!is.na(LITDIA_Normal_1_SN15[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_10_SN15[i-1] <- sum(!is.na(LITDIA_Normal_10_SN15[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_100_SN15[i-1] <- sum(!is.na(LITDIA_Normal_100_SN15[, i]))}
```

## Spectronaut V.16

``` r
#Turbo
for(i in 2:5){
  id_LITDIA_Turbo_1_SN16[i-1] <- sum(!is.na(LITDIA_Turbo_1_SN16[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_10_SN16[i-1] <- sum(!is.na(LITDIA_Turbo_10_SN16[, i]))}
for(i in 2:5){
  id_LITDIA_Turbo_100_SN16[i-1] <- sum(!is.na(LITDIA_Turbo_100_SN16[, i]))}

#Rapid
for(i in 2:5){
  id_LITDIA_Rapid_1_SN16[i-1] <- sum(!is.na(LITDIA_Rapid_1_SN16[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_10_SN16[i-1] <- sum(!is.na(LITDIA_Rapid_10_SN16[, i]))}
for(i in 2:5){
  id_LITDIA_Rapid_100_SN16[i-1] <- sum(!is.na(LITDIA_Rapid_100_SN16[, i]))}

#Normal
for(i in 2:5){
  id_LITDIA_Normal_1_SN16[i-1] <- sum(!is.na(LITDIA_Normal_1_SN16[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_10_SN16[i-1] <- sum(!is.na(LITDIA_Normal_10_SN16[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_100_SN16[i-1] <- sum(!is.na(LITDIA_Normal_100_SN16[, i]))}
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

Information about the scanning mode, the Scectronaut version and the
input are included afterwards.

## Spectronaut V.15

``` r
#Turbo
LITDIA_Turbo_1_SN15_cv <- compute_CV(LITDIA_Turbo_1_SN15)
LITDIA_Turbo_1_SN15_cv$input<- "1 ng"
LITDIA_Turbo_1_SN15_cv$version<- 15
LITDIA_Turbo_1_SN15_cv$scanning_mode<- "Turbo"
LITDIA_Turbo_10_SN15_cv <- compute_CV(LITDIA_Turbo_10_SN15)
LITDIA_Turbo_10_SN15_cv$input<- "10 ng"
LITDIA_Turbo_10_SN15_cv$version<- 15
LITDIA_Turbo_10_SN15_cv$scanning_mode<- "Turbo"
LITDIA_Turbo_100_SN15_cv <- compute_CV(LITDIA_Turbo_100_SN15)
LITDIA_Turbo_100_SN15_cv$input<- "100 ng"
LITDIA_Turbo_100_SN15_cv$version<- 15
LITDIA_Turbo_100_SN15_cv$scanning_mode<- "Turbo"

#Rapid
LITDIA_Rapid_1_SN15_cv <- compute_CV(LITDIA_Rapid_1_SN15)
LITDIA_Rapid_1_SN15_cv$input<- "1 ng"
LITDIA_Rapid_1_SN15_cv$version<- 15
LITDIA_Rapid_1_SN15_cv$scanning_mode<- "Rapid"
LITDIA_Rapid_10_SN15_cv <- compute_CV(LITDIA_Rapid_10_SN15)
LITDIA_Rapid_10_SN15_cv$input<- "10 ng"
LITDIA_Rapid_10_SN15_cv$version<- 15
LITDIA_Rapid_10_SN15_cv$scanning_mode<- "Rapid"
LITDIA_Rapid_100_SN15_cv <- compute_CV(LITDIA_Rapid_100_SN15)
LITDIA_Rapid_100_SN15_cv$input<- "100 ng"
LITDIA_Rapid_100_SN15_cv$version<- 15
LITDIA_Rapid_100_SN15_cv$scanning_mode<- "Rapid"

#Normal
LITDIA_Normal_1_SN15_cv <- compute_CV(LITDIA_Normal_1_SN15)
LITDIA_Normal_1_SN15_cv$input<- "1 ng"
LITDIA_Normal_1_SN15_cv$version<- 15
LITDIA_Normal_1_SN15_cv$scanning_mode<- "Normal"
LITDIA_Normal_10_SN15_cv <- compute_CV(LITDIA_Normal_10_SN15)
LITDIA_Normal_10_SN15_cv$input<- "10 ng"
LITDIA_Normal_10_SN15_cv$version<- 15
LITDIA_Normal_10_SN15_cv$scanning_mode<- "Normal"
LITDIA_Normal_100_SN15_cv <- compute_CV(LITDIA_Normal_100_SN15)
LITDIA_Normal_100_SN15_cv$input<- "100 ng"
LITDIA_Normal_100_SN15_cv$version<- 15
LITDIA_Normal_100_SN15_cv$scanning_mode<- "Normal"
```

## Spectronaut V.16

``` r
#Turbo
LITDIA_Turbo_1_SN16_cv <- compute_CV(LITDIA_Turbo_1_SN16)
LITDIA_Turbo_1_SN16_cv$input<- "1 ng"
LITDIA_Turbo_1_SN16_cv$version<- 16
LITDIA_Turbo_1_SN16_cv$scanning_mode<- "Turbo"
LITDIA_Turbo_10_SN16_cv <- compute_CV(LITDIA_Turbo_10_SN16)
LITDIA_Turbo_10_SN16_cv$input<- "10 ng"
LITDIA_Turbo_10_SN16_cv$version<- 16
LITDIA_Turbo_10_SN16_cv$scanning_mode<- "Turbo"
LITDIA_Turbo_100_SN16_cv <- compute_CV(LITDIA_Turbo_100_SN16)
LITDIA_Turbo_100_SN16_cv$input<- "100 ng"
LITDIA_Turbo_100_SN16_cv$version<- 16
LITDIA_Turbo_100_SN16_cv$scanning_mode<- "Turbo"

#Rapid
LITDIA_Rapid_1_SN16_cv <- compute_CV(LITDIA_Rapid_1_SN16)
LITDIA_Rapid_1_SN16_cv$input<- "1 ng"
LITDIA_Rapid_1_SN16_cv$version<- 16
LITDIA_Rapid_1_SN16_cv$scanning_mode<- "Rapid"
LITDIA_Rapid_5_SN16_cv <- compute_CV(LITDIA_Rapid_1_SN16)
LITDIA_Rapid_10_SN16_cv <- compute_CV(LITDIA_Rapid_10_SN16)
LITDIA_Rapid_10_SN16_cv$input<- "10 ng"
LITDIA_Rapid_10_SN16_cv$version<- 16
LITDIA_Rapid_10_SN16_cv$scanning_mode<- "Rapid"
LITDIA_Rapid_100_SN16_cv <- compute_CV(LITDIA_Rapid_100_SN16)
LITDIA_Rapid_100_SN16_cv$input<- "100 ng"
LITDIA_Rapid_100_SN16_cv$version<- 16
LITDIA_Rapid_100_SN16_cv$scanning_mode<- "Rapid"

#Normal
LITDIA_Normal_1_SN16_cv <- compute_CV(LITDIA_Normal_1_SN16)
LITDIA_Normal_1_SN16_cv$input<- "1 ng"
LITDIA_Normal_1_SN16_cv$version<- 16
LITDIA_Normal_1_SN16_cv$scanning_mode<- "Normal"
LITDIA_Normal_10_SN16_cv <- compute_CV(LITDIA_Normal_10_SN16)
LITDIA_Normal_10_SN16_cv$input<- "10 ng"
LITDIA_Normal_10_SN16_cv$version<- 16
LITDIA_Normal_10_SN16_cv$scanning_mode<- "Normal"
LITDIA_Normal_100_SN16_cv <- compute_CV(LITDIA_Normal_100_SN16)
LITDIA_Normal_100_SN16_cv$input<- "100 ng"
LITDIA_Normal_100_SN16_cv$version<- 16
LITDIA_Normal_100_SN16_cv$scanning_mode<- "Normal"
```

# Join

Data frames containing CV are combined with the function `rbind()`.

``` r
id_figure8 <- 
rbind(LITDIA_Turbo_1_SN15_cv, LITDIA_Turbo_10_SN15_cv, LITDIA_Turbo_100_SN15_cv,
LITDIA_Rapid_1_SN15_cv, LITDIA_Rapid_10_SN15_cv, LITDIA_Rapid_100_SN15_cv,
LITDIA_Normal_1_SN15_cv, LITDIA_Normal_10_SN15_cv, LITDIA_Normal_100_SN15_cv,
LITDIA_Turbo_1_SN16_cv,LITDIA_Turbo_10_SN16_cv, LITDIA_Turbo_100_SN16_cv,
LITDIA_Rapid_1_SN16_cv, LITDIA_Rapid_10_SN16_cv, LITDIA_Rapid_100_SN16_cv,
LITDIA_Normal_1_SN16_cv, LITDIA_Normal_10_SN16_cv, LITDIA_Normal_100_SN16_cv)
# Factorise to set order
id_figure8$input <- factor(id_figure8$input, 
                          levels = unique(id_figure8$input))
id_figure8$scanning_mode <- factor(id_figure8$scanning_mode, 
                          levels = unique(id_figure8$scanning_mode))
id_figure8$version <- factor(id_figure8$version, 
                          levels = unique(id_figure8$version))
```

A new column called “CV\_filter” is added to identify peptides with a CV
over 0.15 (1), between 0.1 and 0.15 (2) and under 0.1 (3).

``` r
id_figure8$CV_filter <- NA
id_figure8$CV_filter[id_figure8$CV < 0.1] <- "3"
id_figure8$CV_filter[id_figure8$CV >= 0.1 & id_figure8$CV <= 0.15] <- "2"
id_figure8$CV_filter[id_figure8$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out. The
number of identified peptides for each replicate and each CV\_filter is
calculated by summing the number of non-missing peptides. The mean and
the standard deviation of the replicates are also calculated for each
method.

``` r
id_figure8 <- id_figure8 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(input, CV_filter, scanning_mode, version) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
```

    ## `summarise()` has grouped output by 'input', 'CV_filter', 'scanning_mode'. You
    ## can override using the `.groups` argument.

``` r
id_figure8$id_S4[id_figure8$id_S4 == 0] <- NA

#Calculating mean and sd from the replicates 
id_figure8$mean_id <- rowMeans(id_figure8[5:8], na.rm = TRUE)
id_figure8$sd_id <- apply(id_figure8[5:8], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

# Plot

The mean and the standard deviation of identified peptides for each
CV\_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y\_errorbar variable instead of
the mean\_id.

``` r
id_figure8$y_errorbar <- id_figure8$mean_id
id_figure8$y_errorbar[id_figure8$CV_filter == "1"] <-
  id_figure8$y_errorbar[id_figure8$CV_filter == "1"] +
  id_figure8$y_errorbar[id_figure8$CV_filter == "2"]+
  id_figure8$y_errorbar[id_figure8$CV_filter == "3"]
id_figure8$y_errorbar[id_figure8$CV_filter == "2"] <-
  id_figure8$y_errorbar[id_figure8$CV_filter == "2"]+
  id_figure8$y_errorbar[id_figure8$CV_filter == "3"]

figure8 <- id_figure8 %>%
  ggplot(aes(x = version, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
  labs(x = "Spectronaut version",
       y = "Identified peptides",
       fill = "CV") +
  facet_grid(input~scanning_mode)+
  theme_bw() +
  theme(axis.text = element_text(size = 24),
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
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure8
```

![](Figure8_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Figure 8 was saved as a PDF.

``` r
cairo_pdf(filename = "fig8", width = 20.83, height = 25)
figure8
dev.off()
```

    ## png 
    ##   2

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
    ##  [1] patchwork_1.1.2 forcats_0.5.2   stringr_1.4.1   dplyr_1.0.10   
    ##  [5] purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.8   
    ##  [9] ggplot2_3.3.6   tidyverse_1.3.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2    xfun_0.32           haven_2.5.1        
    ##  [4] gargle_1.2.1        colorspace_2.0-3    vctrs_0.4.1        
    ##  [7] generics_0.1.3      htmltools_0.5.3     yaml_2.3.5         
    ## [10] utf8_1.2.2          rlang_1.0.5         pillar_1.8.1       
    ## [13] withr_2.5.0         glue_1.6.2          DBI_1.1.3          
    ## [16] dbplyr_2.2.1        modelr_0.1.9        readxl_1.4.1       
    ## [19] lifecycle_1.0.1     munsell_0.5.0       gtable_0.3.1       
    ## [22] cellranger_1.1.0    rvest_1.0.3         evaluate_0.16      
    ## [25] labeling_0.4.2      knitr_1.40          tzdb_0.3.0         
    ## [28] fastmap_1.1.0       fansi_1.0.3         highr_0.9          
    ## [31] broom_1.0.1         scales_1.2.1        backports_1.4.1    
    ## [34] googlesheets4_1.0.1 jsonlite_1.8.0      farver_2.1.1       
    ## [37] fs_1.5.2            hms_1.1.2           digest_0.6.29      
    ## [40] stringi_1.7.8       grid_4.2.1          cli_3.3.0          
    ## [43] tools_4.2.1         magrittr_2.0.3      crayon_1.5.1       
    ## [46] pkgconfig_2.0.3     ellipsis_0.3.2      xml2_1.3.3         
    ## [49] reprex_2.0.2        googledrive_2.0.0   lubridate_1.8.0    
    ## [52] assertthat_0.2.1    rmarkdown_2.16      httr_1.4.4         
    ## [55] rstudioapi_0.14     R6_2.5.1            compiler_4.2.1
