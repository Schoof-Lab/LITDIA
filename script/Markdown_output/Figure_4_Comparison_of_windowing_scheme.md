Fig. 4 Comparison of windowing scheme.
================
Samuel Grégoire & Lukas R. Woltereck

Fig. 4 Comparison of windowing scheme.

Comparison of the number of identified peptides on LIT-DIA-based method
on Normal scanning mode with constant auto-IIT and different numbers of
windows (34, 40, and 45 windows) and input material (low, and high-input
tryptic HeLa digest) on WhisperTM 100 for 40 SPD and 20 SPD. Identified
peptides with a coefficient of variation (CV) between 10% and 15% are
coloured with light red and those with a CV below 10% with dark red.

# Load packages

``` r
library("tidyverse")

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

## 20 SPD

``` r
#1 ng
LITDIA_Normal_20SPD_34w_1 <- read.csv2(paste0(lit_path, "20220218_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV_auto_34w.csv"))
LITDIA_Normal_20SPD_40w_1 <- read.csv2(paste0(lit_path, "20220502_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_20SPD_45w_1 <- read.csv2(paste0(lit_path, "20220218_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV_auto_45w.csv"))

#100 ng
LITDIA_Normal_20SPD_34w_100 <- read.csv2(paste0(lit_path, "20220601_LIT_Normal_DIA_100ng_20SPD_whisper100_1CV_auto_34w.csv"))
LITDIA_Normal_20SPD_40w_100 <- read.csv2(paste0(lit_path, "20220413_LIT_Normal_DIA_100ng_20SPD_whisper100_1CV_auto_40w.csv"))
LITDIA_Normal_20SPD_45w_100 <- read.csv2(paste0(lit_path, "20220418_LIT_Normal_DIA_100ng_20SPD_whisper100_1CV_auto_45w.csv"))
```

## 40SPD

``` r
#1 ng
LITDIA_Normal_34w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_1ng_40SPD_whisper100_1CV_auto_34w.csv"))
LITDIA_Normal_40w_1 <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_1ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_45w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_1ng_40SPD_whisper100_1CV_auto_45w.csv"))

#100 ng
LITDIA_Normal_34w_100 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_100ng_40SPD_whisper100_1CV_auto_34w.csv"))
LITDIA_Normal_40w_100 <- read.csv2(paste0(lit_path, "20220421_LIT_Normal_DIA_100ng_40SPD_whisper100_1CV_auto.csv"))
LITDIA_Normal_45w_100 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_100ng_40SPD_whisper100_1CV_auto_45w.csv"))
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

## 20 SPD

``` r
#1 ng
LITDIA_Normal_20SPD_34w_1 <-  filter_pep(LITDIA_Normal_20SPD_34w_1)
```

    ## 
    ## Highlighted 15 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_20SPD_40w_1 <-  filter_pep(LITDIA_Normal_20SPD_40w_1)
```

    ## 
    ## Highlighted 39 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_20SPD_45w_1 <-  filter_pep(LITDIA_Normal_20SPD_45w_1)
```

    ## 
    ## Highlighted 79 peptide(s) found in less than 3 replicates

``` r
#100 ng
LITDIA_Normal_20SPD_34w_100 <-  filter_pep(LITDIA_Normal_20SPD_34w_100)
```

    ## 
    ## Highlighted 39 peptide(s) found in less than 2 replicates

``` r
LITDIA_Normal_20SPD_40w_100 <-  filter_pep(LITDIA_Normal_20SPD_40w_100)
```

    ## 
    ## Highlighted 99 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_20SPD_45w_100 <-  filter_pep(LITDIA_Normal_20SPD_45w_100)
```

    ## 
    ## Highlighted 56 peptide(s) found in less than 2 replicates

## 40 SPD

``` r
#1 ng
LITDIA_Normal_34w_1 <-  filter_pep(LITDIA_Normal_34w_1)
```

    ## 
    ## Highlighted 11 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_40w_1 <-  filter_pep(LITDIA_Normal_40w_1)
```

    ## 
    ## Highlighted 41 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_45w_1 <-  filter_pep(LITDIA_Normal_45w_1)
```

    ## 
    ## Highlighted 19 peptide(s) found in less than 3 replicates

``` r
#100 ng
LITDIA_Normal_34w_100 <-  filter_pep(LITDIA_Normal_34w_100)
```

    ## 
    ## Highlighted 204 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_40w_100 <-  filter_pep(LITDIA_Normal_40w_100)
```

    ## 
    ## Highlighted 316 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_45w_100 <-  filter_pep(LITDIA_Normal_45w_100)
```

    ## 
    ## Highlighted 405 peptide(s) found in less than 3 replicates

# Count number of identified peptides

Empty vectors are created to store the number of identified peptides in
each replicates for each data-set.

``` r
id_LITDIA_Normal_20SPD_34w_1 <- id_LITDIA_Normal_20SPD_40w_1 <-
  id_LITDIA_Normal_20SPD_45w_1 <- id_LITDIA_Normal_20SPD_34w_100 <- 
  id_LITDIA_Normal_20SPD_40w_100 <- id_LITDIA_Normal_20SPD_45w_100 <- 
  id_LITDIA_Normal_34w_1 <- id_LITDIA_Normal_40w_1 <-
  id_LITDIA_Normal_45w_1 <- id_LITDIA_Normal_34w_100 <-
  id_LITDIA_Normal_40w_100 <- id_LITDIA_Normal_45w_100 <- rep(NA, 4)
```

The numbers of identified peptides are counted by summing non-missing
values for each replicate and stored in the vectors created above.

## 20 SPD

``` r
#1 ng
for(i in 2:5){
  id_LITDIA_Normal_20SPD_34w_1[i-1] <- 
    sum(!is.na(LITDIA_Normal_20SPD_34w_1[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_20SPD_40w_1[i-1] <- 
    sum(!is.na(LITDIA_Normal_20SPD_40w_1[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_20SPD_45w_1[i-1] <- 
    sum(!is.na(LITDIA_Normal_20SPD_45w_1[, i]))}

#100 ng
for(i in 2:5){
  id_LITDIA_Normal_20SPD_34w_100[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_34w_100[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_20SPD_40w_100[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_40w_100[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_20SPD_45w_100[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_45w_100[, i]))}
```

## 40 SPD

``` r
#1 ng
for(i in 2:5){
  id_LITDIA_Normal_34w_1[i-1] <- sum(!is.na(LITDIA_Normal_34w_1[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_40w_1[i-1] <- sum(!is.na(LITDIA_Normal_40w_1[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_45w_1[i-1] <- sum(!is.na(LITDIA_Normal_45w_1[, i]))}

#100 ng
for(i in 2:5){
  id_LITDIA_Normal_34w_100[i-1] <- sum(!is.na(LITDIA_Normal_34w_100[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_40w_100[i-1] <- sum(!is.na(LITDIA_Normal_40w_100[, i]))}
for(i in 2:5){
  id_LITDIA_Normal_45w_100[i-1] <- sum(!is.na(LITDIA_Normal_45w_100[, i]))}
```

# Compute CV

The function `compute_CV()` below takes the data-sets, already filtered
by `filter_pep()`, as an entry and compute the coefficients of variation
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

Information about the method, the number of samples per day and the
input are included afterwards.

## 20 SPD

``` r
LITDIA_Normal_20SPD_34w_1_cv <- compute_CV(LITDIA_Normal_20SPD_34w_1)
LITDIA_Normal_20SPD_34w_1_cv$method <- 34
LITDIA_Normal_20SPD_34w_1_cv$spd<- "20 SPD"
LITDIA_Normal_20SPD_34w_1_cv$input<- "1 ng"
LITDIA_Normal_20SPD_34w_100_cv <- compute_CV(LITDIA_Normal_20SPD_34w_100)
LITDIA_Normal_20SPD_34w_100_cv$method <- 34
LITDIA_Normal_20SPD_34w_100_cv$spd<- "20 SPD"
LITDIA_Normal_20SPD_34w_100_cv$input<- "100 ng"
LITDIA_Normal_20SPD_40w_1_cv <- compute_CV(LITDIA_Normal_20SPD_40w_1)
LITDIA_Normal_20SPD_40w_1_cv$method <- 40
LITDIA_Normal_20SPD_40w_1_cv$spd<- "20 SPD"
LITDIA_Normal_20SPD_40w_1_cv$input<- "1 ng"
LITDIA_Normal_20SPD_40w_100_cv <- compute_CV(LITDIA_Normal_20SPD_40w_100)
LITDIA_Normal_20SPD_40w_100_cv$method <- 40
LITDIA_Normal_20SPD_40w_100_cv$spd<- "20 SPD"
LITDIA_Normal_20SPD_40w_100_cv$input<- "100 ng"
LITDIA_Normal_20SPD_45w_1_cv <- compute_CV(LITDIA_Normal_20SPD_45w_1)
LITDIA_Normal_20SPD_45w_1_cv$method <- 45
LITDIA_Normal_20SPD_45w_1_cv$spd<- "20 SPD"
LITDIA_Normal_20SPD_45w_1_cv$input<- "1 ng"
LITDIA_Normal_20SPD_45w_100_cv <- compute_CV(LITDIA_Normal_20SPD_45w_100)
LITDIA_Normal_20SPD_45w_100_cv$method <- 45
LITDIA_Normal_20SPD_45w_100_cv$spd<- "20 SPD"
LITDIA_Normal_20SPD_45w_100_cv$input<- "100 ng"
```

## 40 SPD

``` r
LITDIA_Normal_34w_1_cv <- compute_CV(LITDIA_Normal_34w_1)
LITDIA_Normal_34w_1_cv$method <- 34
LITDIA_Normal_34w_1_cv$spd<- "40 SPD"
LITDIA_Normal_34w_1_cv$input<- "1 ng"
LITDIA_Normal_34w_100_cv <- compute_CV(LITDIA_Normal_34w_100)
LITDIA_Normal_34w_100_cv$method <- 34
LITDIA_Normal_34w_100_cv$spd<- "40 SPD"
LITDIA_Normal_34w_100_cv$input<- "100 ng"
LITDIA_Normal_40w_1_cv <- compute_CV(LITDIA_Normal_40w_1)
LITDIA_Normal_40w_1_cv$method <- 40
LITDIA_Normal_40w_1_cv$spd<- "40 SPD"
LITDIA_Normal_40w_1_cv$input<- "1 ng"
LITDIA_Normal_40w_100_cv <- compute_CV(LITDIA_Normal_40w_100)
LITDIA_Normal_40w_100_cv$method <- 40
LITDIA_Normal_40w_100_cv$spd<- "40 SPD"
LITDIA_Normal_40w_100_cv$input<- "100 ng"
LITDIA_Normal_45w_1_cv <- compute_CV(LITDIA_Normal_45w_1)
LITDIA_Normal_45w_1_cv$method <- 45
LITDIA_Normal_45w_1_cv$spd<- "40 SPD"
LITDIA_Normal_45w_1_cv$input<- "1 ng"
LITDIA_Normal_45w_100_cv <- compute_CV(LITDIA_Normal_34w_100)
LITDIA_Normal_45w_100_cv$method <- 45
LITDIA_Normal_45w_100_cv$spd<- "40 SPD"
LITDIA_Normal_45w_100_cv$input<- "100 ng"
```

# Join

Data frames containing CV are combined with the function `rbind()`.

``` r
join_fig4 <- 
  rbind(LITDIA_Normal_20SPD_34w_1_cv, LITDIA_Normal_20SPD_34w_100_cv, 
        LITDIA_Normal_20SPD_40w_1_cv, LITDIA_Normal_20SPD_40w_100_cv, 
        LITDIA_Normal_20SPD_45w_1_cv, LITDIA_Normal_20SPD_45w_100_cv, 
        LITDIA_Normal_34w_1_cv, LITDIA_Normal_34w_100_cv, 
        LITDIA_Normal_40w_1_cv, LITDIA_Normal_40w_100_cv, 
        LITDIA_Normal_45w_1_cv, LITDIA_Normal_45w_100_cv)

join_fig4$method %>% table()
```

    ## .
    ##    34    40    45 
    ## 28249 35745 31095

``` r
join_fig4$input %>% table()
```

    ## .
    ##   1 ng 100 ng 
    ##  27863  67226

``` r
join_fig4$spd %>% table()
```

    ## .
    ## 20 SPD 40 SPD 
    ##  55878  39211

``` r
# Factorise to set order
join_fig4$method <- factor(join_fig4$method, 
                          levels = unique(join_fig4$method))
join_fig4$input <- factor(join_fig4$input, 
                          levels = unique(join_fig4$input))
join_fig4$spd <- factor(join_fig4$spd, 
                          levels = c("40 SPD", "20 SPD"))
```

A new column called “CV_filter” is added to identify peptides with a CV
over 0.15 (1), between 0.1 and 0.15 (2) and under 0.1 (3).

``` r
join_fig4$CV_filter <- NA
join_fig4$CV_filter[join_fig4$CV < 0.1] <- "3"
join_fig4$CV_filter[join_fig4$CV >= 0.1 & join_fig4$CV <= 0.15] <- "2"
join_fig4$CV_filter[join_fig4$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out. The
number of identified peptides for each replicate and each CV_filter is
calculated by summing the number of non-missing peptides. The mean and
the standard deviation of the replicates are also calculated for each
method.

``` r
join_fig4 <- join_fig4 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(method, spd, input, CV_filter) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
```

    ## `summarise()` has grouped output by 'method', 'spd', 'input'. You can override
    ## using the `.groups` argument.

``` r
join_fig4$id_S4[join_fig4$id_S4 == 0] <- NA 
join_fig4$mean_id <- rowMeans(join_fig4[5:8], na.rm = TRUE)
join_fig4$sd_id <- apply(join_fig4[5:8], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

# Plot

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

``` r
join_fig4$y_errorbar <- join_fig4$mean_id
join_fig4$y_errorbar[join_fig4$CV_filter == "1"] <-
  join_fig4$y_errorbar[join_fig4$CV_filter == "1"] +
  join_fig4$y_errorbar[join_fig4$CV_filter == "2"]+
  join_fig4$y_errorbar[join_fig4$CV_filter == "3"]
join_fig4$y_errorbar[join_fig4$CV_filter == "2"] <-
  join_fig4$y_errorbar[join_fig4$CV_filter == "2"]+
  join_fig4$y_errorbar[join_fig4$CV_filter == "3"]
figure4 <- join_fig4 %>% 
  ggplot(aes(x = method, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
 facet_grid(input~spd)+
  labs(x = "Number of Windows",
       y = "Identified peptides",
       fill = "CV") +
  theme_bw() +
  theme(axis.text = element_text(size = 24),
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
figure4
```

![](Figure4_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
#Save as pdf
#cairo_pdf(filename = "fig4.pdf", width = 18.75, height = 12.5)
#figure4
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
    ## [1] forcats_0.5.2   stringr_1.4.1   dplyr_1.0.10    purrr_0.3.4    
    ## [5] readr_2.1.2     tidyr_1.2.0     tibble_3.1.8    ggplot2_3.3.6  
    ## [9] tidyverse_1.3.2
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
