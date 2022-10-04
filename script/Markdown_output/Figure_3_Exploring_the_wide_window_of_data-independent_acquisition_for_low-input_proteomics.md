Fig. 3 Exploring the wide window of data-independent acquisition for
low-input proteomics.
================
Samuel Grégoire & Lukas R. Woltereck

Fig. 3 Exploring the wide window of data-independent acquisition for
low-input proteomics.

Comparison between LIT-DIA with 10 m/z isolation windows and OT-DIA
methods (10 m/z isolation windows on 15k resolution, 20 m/z isolation
windows on 30k resolution and various windows sizes (4 windows with
120-, 120-, 200-, and 580 m/z) on 120k resolution)

# Load packages

``` r
library("tidyverse")
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

## LIT Normal 10mz

``` r
LIT_Normal_10mz_1 <- 
  read.csv2(paste0(lit_path, "20220924_HB_Evo_Whisper100_40SPD_LITDIA_Normal_40w_10mz_auto_HeLa_1ng_directedDIA_Report.csv"))

LIT_Normal_10mz_5 <- 
  read.csv2(paste0(lit_path, "20220924_EV_Evo_Whisper100_40SPD_LITDIA_Normal_40w_10mz_auto_HeLa_5ng_directedDIA_Report.csv"))

LIT_Normal_10mz_10 <- 
  read.csv2(paste0(lit_path, "20220924_EV_Evo_Whisper100_40SPD_LITDIA_Normal_40w_10mz_auto_HeLa_10ng_directedDIA_Report.csv"))

LIT_Normal_10mz_100 <- 
  read.csv2(paste0(lit_path, "20220924_HB_Evo_Whisper100_40SPD_LITDIA_Normal_40w_10mz_auto_HeLa_100ng_directedDIA_Report.csv"))
```

## OT 15k 10mz

``` r
OT_15k_10mz_1 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_40w_10mz_15k_auto_HeLa_1ng_directedDIA_Report.csv"))

OT_15k_10mz_5 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_40w_10mz_15k_auto_HeLa_5ng_directedDIA_Report.csv"))

OT_15k_10mz_10 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_40w_10mz_15k_auto_HeLa_10ng_directedDIA_Report.csv"))

OT_15k_10mz_100 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_40w_10mz_15k_auto_HeLa_100ng_directedDIA_Report.csv"))
```

## OT 30k 20mz

``` r
OT_30k_20mz_1 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_20w_20mz_30k_auto_1ng_HeLa_directedDIA_Report.csv"))

OT_30k_20mz_5 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_20w_20mz_30k_auto_5ng_HeLa_directedDIA_Report.csv"))

OT_30k_20mz_10 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_20w_20mz_30k_auto_10ng_HeLa_directedDIA_Report.csv"))

OT_30k_20mz_100 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_20w_20mz_30k_auto_100ng_HeLa_directedDIA_Report.csv"))
```

## OT 120k varmz

``` r
OT_120k_varmz_1 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_4w_varmz_120k_auto_HeLa_1ng_directedDIA_Report.csv"))

OT_120k_varmz_5 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_4w_varmz_120k_auto_HeLa_5ng_directedDIA_Report.csv"))

OT_120k_varmz_10 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_4w_varmz_120k_auto_HeLa_10ng_directedDIA_Report.csv"))

OT_120k_varmz_100 <- 
  read.csv2(paste0(ot_path, "20220924_HB_Evo_Whisper100_40SPD_OTDIA_4w_varmz_120k_auto_HeLa_100ng_directedDIA_Report.csv"))
```

# Clean data tables

Imported tables contain a lot of information, which is not useful for
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
Highlighted ", sum(x_filter_highlight$filtered) , " peptide(s) found in less than ", ncol(x_clean_with_mean) -1, " replicates")
  
  return(x_filter_highlight)
}
```

The `filter_pep()` function is used for every data-set.

## LIT Normal 10mz

``` r
LIT_Normal_10mz_1 <- filter_pep(LIT_Normal_10mz_1)
```

    ## 
    ## Highlighted 51 peptide(s) found in less than 4 replicates

``` r
LIT_Normal_10mz_5 <- filter_pep(LIT_Normal_10mz_5)
```

    ## 
    ## Highlighted 62 peptide(s) found in less than 4 replicates

``` r
LIT_Normal_10mz_10 <- filter_pep(LIT_Normal_10mz_10)
```

    ## 
    ## Highlighted 88 peptide(s) found in less than 4 replicates

``` r
LIT_Normal_10mz_100 <- filter_pep(LIT_Normal_10mz_100)
```

    ## 
    ## Highlighted 393 peptide(s) found in less than 4 replicates

## OT 15k 10mz

``` r
OT_15k_10mz_1 <- filter_pep(OT_15k_10mz_1)
```

    ## 
    ## Highlighted 0 peptide(s) found in less than 3 replicates

``` r
OT_15k_10mz_5 <- filter_pep(OT_15k_10mz_5)
```

    ## 
    ## Highlighted 2 peptide(s) found in less than 3 replicates

``` r
OT_15k_10mz_10 <- filter_pep(OT_15k_10mz_10)
```

    ## 
    ## Highlighted 5 peptide(s) found in less than 3 replicates

``` r
OT_15k_10mz_100 <- filter_pep(OT_15k_10mz_100)
```

    ## 
    ## Highlighted 43 peptide(s) found in less than 3 replicates

## OT 30k 20mz

``` r
OT_30k_20mz_1 <- filter_pep(OT_30k_20mz_1)
```

    ## 
    ## Highlighted 3 peptide(s) found in less than 3 replicates

``` r
OT_30k_20mz_5 <- filter_pep(OT_30k_20mz_5)
```

    ## 
    ## Highlighted 6 peptide(s) found in less than 3 replicates

``` r
OT_30k_20mz_10 <- filter_pep(OT_30k_20mz_10)
```

    ## 
    ## Highlighted 16 peptide(s) found in less than 3 replicates

``` r
OT_30k_20mz_100 <- filter_pep(OT_30k_20mz_100)
```

    ## 
    ## Highlighted 36 peptide(s) found in less than 3 replicates

## OT 120k varmz

``` r
OT_120k_varmz_1 <- filter_pep(OT_120k_varmz_1)
```

    ## 
    ## Highlighted 9 peptide(s) found in less than 3 replicates

``` r
OT_120k_varmz_5 <- filter_pep(OT_120k_varmz_5)
```

    ## 
    ## Highlighted 8 peptide(s) found in less than 3 replicates

``` r
OT_120k_varmz_10 <- filter_pep(OT_120k_varmz_10)
```

    ## 
    ## Highlighted 8 peptide(s) found in less than 3 replicates

``` r
OT_120k_varmz_100 <- filter_pep(OT_120k_varmz_100)
```

    ## 
    ## Highlighted 8 peptide(s) found in less than 3 replicates

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

Information about the scanning mode/resolution, and window size are
included afterwards.

## LIT Normal 10mz

``` r
LIT_Normal_10mz_1_cv <- compute_CV(LIT_Normal_10mz_1)
LIT_Normal_10mz_1_cv$method <- "LIT Normal 10 m/z"
LIT_Normal_10mz_1_cv$input<- 1

LIT_Normal_10mz_5_cv <- compute_CV(LIT_Normal_10mz_5)
LIT_Normal_10mz_5_cv$method <- "LIT Normal 10 m/z"
LIT_Normal_10mz_5_cv$input<- 5

LIT_Normal_10mz_10_cv <- compute_CV(LIT_Normal_10mz_10)
LIT_Normal_10mz_10_cv$method <- "LIT Normal 10 m/z"
LIT_Normal_10mz_10_cv$input<- 10

LIT_Normal_10mz_100_cv <- compute_CV(LIT_Normal_10mz_100)
LIT_Normal_10mz_100_cv$method <- "LIT Normal 10 m/z"
LIT_Normal_10mz_100_cv$input<- 100
```

## OT 15k 10mz

``` r
OT_15k_10mz_1_cv <- compute_CV(OT_15k_10mz_1)
OT_15k_10mz_1_cv$method <- "OT 15k 10 m/z"
OT_15k_10mz_1_cv$input<- 1

OT_15k_10mz_5_cv <- compute_CV(OT_15k_10mz_5)
OT_15k_10mz_5_cv$method <- "OT 15k 10 m/z"
OT_15k_10mz_5_cv$input<- 5

OT_15k_10mz_10_cv <- compute_CV(OT_15k_10mz_10)
OT_15k_10mz_10_cv$method <- "OT 15k 10 m/z"
OT_15k_10mz_10_cv$input<- 10

OT_15k_10mz_100_cv <- compute_CV(OT_15k_10mz_100)
OT_15k_10mz_100_cv$method <- "OT 15k 10 m/z"
OT_15k_10mz_100_cv$input<- 100
```

## OT 30k 20mz

``` r
OT_30k_20mz_1_cv <- compute_CV(OT_30k_20mz_1)
OT_30k_20mz_1_cv$method <- "OT 30k 20 m/z"
OT_30k_20mz_1_cv$input<- 1

OT_30k_20mz_5_cv <- compute_CV(OT_30k_20mz_5)
OT_30k_20mz_5_cv$method <- "OT 30k 20 m/z"
OT_30k_20mz_5_cv$input<- 5

OT_30k_20mz_10_cv <- compute_CV(OT_30k_20mz_10)
OT_30k_20mz_10_cv$method <- "OT 30k 20 m/z"
OT_30k_20mz_10_cv$input<- 10

OT_30k_20mz_100_cv <- compute_CV(OT_30k_20mz_100)
OT_30k_20mz_100_cv$method <- "OT 30k 20 m/z"
OT_30k_20mz_100_cv$input<- 100
```

## OT 120k varmz

``` r
OT_120k_varmz_1_cv <- compute_CV(OT_120k_varmz_1)
OT_120k_varmz_1_cv$method <- "OT 120k var m/z"
OT_120k_varmz_1_cv$input<- 1

OT_120k_varmz_5_cv <- compute_CV(OT_120k_varmz_5)
OT_120k_varmz_5_cv$method <- "OT 120k var m/z"
OT_120k_varmz_5_cv$input<- 5

OT_120k_varmz_10_cv <- compute_CV(OT_120k_varmz_10)
OT_120k_varmz_10_cv$method <- "OT 120k var m/z"
OT_120k_varmz_10_cv$input<- 10

OT_120k_varmz_100_cv <- compute_CV(OT_120k_varmz_100)
OT_120k_varmz_100_cv$method <- "OT 120k var m/z"
OT_120k_varmz_100_cv$input<- 100
```

# Join

Data frames containing CV are combined with the function `rbind()` and
leveled with the function `factor()`.

``` r
id_figure3 <- 
  rbind(LIT_Normal_10mz_1_cv, LIT_Normal_10mz_5_cv, LIT_Normal_10mz_10_cv, LIT_Normal_10mz_100_cv,
        OT_15k_10mz_1_cv, OT_15k_10mz_5_cv, OT_15k_10mz_10_cv, OT_15k_10mz_100_cv,
        OT_30k_20mz_1_cv, OT_30k_20mz_5_cv, OT_30k_20mz_10_cv, OT_30k_20mz_100_cv, 
        OT_120k_varmz_1_cv, OT_120k_varmz_5_cv, OT_120k_varmz_10_cv, OT_120k_varmz_100_cv)

id_figure3$method %>% table()
```

    ## .
    ## LIT Normal 10 m/z   OT 120k var m/z     OT 15k 10 m/z     OT 30k 20 m/z 
    ##             47660             47269             42695             56726

``` r
id_figure3$input %>% table()
```

    ## .
    ##     1     5    10   100 
    ## 22471 45938 53564 72377

``` r
id_figure3$method <- factor(id_figure3$method, 
                          levels = unique(id_figure3$method))
```

A new column called “CV_filter” is added to identify peptides with a CV
over 0.15 (1), CV between 0.1 and 0.15 (2) and CV under 0.1 (3).

``` r
id_figure3$CV_filter <- NA
id_figure3$CV_filter[id_figure3$CV < 0.1] <- "3"
id_figure3$CV_filter[id_figure3$CV >= 0.1 & id_figure3$CV <= 0.15] <- "2"
id_figure3$CV_filter[id_figure3$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out and
the number of identified peptides for each replicate and each CV_filter
is calculated by summing the number of non-missing peptides. The mean
and the standard deviation of the replicates are also calculated for
each method.

``` r
id_figure3 <- id_figure3 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(method, input, CV_filter) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
```

    ## `summarise()` has grouped output by 'method', 'input'. You can override using
    ## the `.groups` argument.

``` r
id_figure3$id_S4[id_figure3$id_S4 == 0] <- NA

id_figure3$mean_id <- rowMeans(id_figure3[4:7], na.rm = TRUE)
id_figure3$sd_id <- apply(id_figure3[4:7], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

# Plot

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

## Turbo

``` r
id_figure3$y_errorbar <- id_figure3$mean_id
id_figure3$y_errorbar[id_figure3$CV_filter == "1"] <-
  id_figure3$y_errorbar[id_figure3$CV_filter == "1"] +
  id_figure3$y_errorbar[id_figure3$CV_filter == "2"]+
  id_figure3$y_errorbar[id_figure3$CV_filter == "3"]
id_figure3$y_errorbar[id_figure3$CV_filter == "2"] <-
  id_figure3$y_errorbar[id_figure3$CV_filter == "2"]+
  id_figure3$y_errorbar[id_figure3$CV_filter == "3"]

figure3 <- id_figure3 %>% 
  ggplot(aes(x = as.factor(input), y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
 facet_grid(~method)+
  labs(x = "Input (ng)",
       y = "Identified peptides",
       fill = "CV") +
  theme_bw() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 27),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 27),
        legend.key.size = unit(0.9, "cm"),
        legend.position = "top")+
  scale_fill_manual(values = c("grey82", "#ffaa8b", "#ff713e"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 24000))+
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure3
```

![](fig3_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
#save as pdf
#cairo_pdf(filename = "fig3.pdf", width = 20.83, height = 9.375)
#figure3
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
