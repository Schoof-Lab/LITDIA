Sup. Fig. 1 Optimizing a balance between Windowing Schemes and ITs
================

Comparison of the number of identified peptides on DIA-LIT-based method
on Normal scanning mode for WhisperTM 20 SPD and Rapid scanning mode for
WhisperTM 40 SPD from 1 ng of tryptic HeLa lysate with different numbers
of windows and IT at fixed cycle time. Identified peptides with a
coefficient of variation (CV) between 10% and 15% are colored with light
red and those with a CV below 10% with dark red.

# Load packages

``` r
library("tidyverse")
library("psych")
library("pheatmap")
library("ggpubr")
library("patchwork")
library("ghibli")
library("reprex")
library("ggh4x")
library("ggplot2")
```

# Load data

The path to every folder used is stored to keep everything compact. LIT
data is stored in the folder located in lit_path. To reproduce analysis
the path must be changed to corresponding path on the local computer.

``` r
lit_path <- "C://data/lit"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in RStudio with the function
`read.csv2()`.

40 SPD

``` r
#38ms 34w
LITDIA_Rapid_40SPD_23ms_40w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_auto_34w.csv"))
#30ms 30w
LITDIA_Rapid_40SPD_30ms_30w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_30ms_30w.csv"))
```

20 SPD

``` r
#76ms 20w
LITDIA_Normal_20SPD_76ms_20w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV_76ms_20w.csv"))
#50ms 30w
LITDIA_Normal_20SPD_50ms_30w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV_50ms_30w.csv"))
#38ms 40w
LITDIA_Normal_20SPD_38ms_40w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV_38ms_40w.csv"))
```

# Clean data tables

Imported tables contain a lot of information, which are not useful for
our analysis. The function `filter_pep()` only keeps the file names
(replicates), quantitative values for each file name, and the peptide
sequence. It is also calculating the mean between the four replicates.
Spectronaut does some imputation by default , but we want to keep
missing values. Imputed values are replaced by *NA*. Peptides missing in
more than one replicate are highlighted as *TRUE*.

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

The `filter_pep()` function is used for every data-set.

40 SPD

``` r
LITDIA_Rapid_40SPD_23ms_40w_1 <- filter_pep(LITDIA_Rapid_40SPD_23ms_40w_1)
```

    ## 
    ## Highlighted 16 peptide(s) found in less than 2 replicates

``` r
LITDIA_Rapid_40SPD_30ms_30w_1 <- filter_pep(LITDIA_Rapid_40SPD_30ms_30w_1)
```

    ## 
    ## Highlighted 29 peptide(s) found in less than 2 replicates

20 SPD

``` r
LITDIA_Normal_20SPD_76ms_20w_1 <- filter_pep(LITDIA_Normal_20SPD_76ms_20w_1)
```

    ## 
    ## Highlighted 5 peptide(s) found in less than 2 replicates

``` r
LITDIA_Normal_20SPD_50ms_30w_1 <- filter_pep(LITDIA_Normal_20SPD_50ms_30w_1)
```

    ## 
    ## Highlighted 12 peptide(s) found in less than 2 replicates

``` r
LITDIA_Normal_20SPD_38ms_40w_1 <- filter_pep(LITDIA_Normal_20SPD_38ms_40w_1)
```

    ## 
    ## Highlighted 8 peptide(s) found in less than 2 replicates

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

Information about the scanning mode, the windowing scheme and the
numbers of samples per day are included afterwards.

40 SPD

``` r
LITDIA_Rapid_40SPD_23ms_40w_1_CV <- compute_CV(LITDIA_Rapid_40SPD_23ms_40w_1)
LITDIA_Rapid_40SPD_23ms_40w_1_CV$method <- "23ms_40w"
LITDIA_Rapid_40SPD_23ms_40w_1_CV$SPD <- "40 SPD"
LITDIA_Rapid_40SPD_30ms_30w_1_CV <- compute_CV(LITDIA_Rapid_40SPD_30ms_30w_1)
LITDIA_Rapid_40SPD_30ms_30w_1_CV$method <- "30ms_30w"
LITDIA_Rapid_40SPD_30ms_30w_1_CV$SPD <- "40 SPD"
```

20 SPD

``` r
LITDIA_Normal_20SPD_76ms_20w_1_CV <- compute_CV(LITDIA_Normal_20SPD_76ms_20w_1)
LITDIA_Normal_20SPD_76ms_20w_1_CV$method <- "76ms_20w"
LITDIA_Normal_20SPD_76ms_20w_1_CV$SPD <- "20 SPD"
LITDIA_Normal_20SPD_50ms_30w_1_CV <- compute_CV(LITDIA_Normal_20SPD_50ms_30w_1)
LITDIA_Normal_20SPD_50ms_30w_1_CV$method <- "50ms_30w"
LITDIA_Normal_20SPD_50ms_30w_1_CV$SPD <- "20 SPD"
LITDIA_Normal_20SPD_38ms_40w_1_CV <- compute_CV(LITDIA_Normal_20SPD_38ms_40w_1)
LITDIA_Normal_20SPD_38ms_40w_1_CV$method <- "38ms_40w"
LITDIA_Normal_20SPD_38ms_40w_1_CV$SPD <- "20 SPD"
```

# Join

Data frames containing CV are combined with the function `rbind()` and
arranged with the function `factor()`.

``` r
CV_figure_sub1 <- 
  rbind(LITDIA_Rapid_40SPD_23ms_40w_1_CV, LITDIA_Rapid_40SPD_30ms_30w_1_CV, 
        LITDIA_Normal_20SPD_38ms_40w_1_CV, LITDIA_Normal_20SPD_50ms_30w_1_CV,
        LITDIA_Normal_20SPD_76ms_20w_1_CV)
CV_figure_sub1$method %>% table()
```

    ## .
    ## 23ms_40w 30ms_30w 38ms_40w 50ms_30w 76ms_20w 
    ##     4594     3997     7023     6365     5406

``` r
CV_figure_sub1$CV_filter <- NA
```

A new column called “CV_filter” is added to identify peptides with a CV
over 0.15 (1), CV between 0.1 and 0.15 (2) and CV under 0.1 (3).

``` r
CV_figure_sub1$CV_filter[CV_figure_sub1$CV < 0.1] <- "3"
CV_figure_sub1$CV_filter[CV_figure_sub1$CV >= 0.1 & CV_figure_sub1$CV <= 0.15] <- "2"
CV_figure_sub1$CV_filter[CV_figure_sub1$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out and
the number of identified peptides for each replicate and each CV_filter
is calculated by summing the number of non-missing peptides. The mean
and the standard deviation of the replicates are also calculated for
each method.

``` r
CV_figure_sub1 <- CV_figure_sub1 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(method, CV_filter, SPD) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
CV_figure_sub1$id_S4[CV_figure_sub1$id_S4 == 0] <- NA
CV_figure_sub1$mean_id <- rowMeans(CV_figure_sub1[4:7], na.rm = TRUE)
CV_figure_sub1$sd_id <- apply(CV_figure_sub1[4:7], FUN =  sd, MARGIN = 1, na.rm = TRUE)
CV_figure_sub1
```

    ## # A tibble: 15 × 9
    ## # Groups:   method, CV_filter [15]
    ##    method   CV_filter SPD    id_S1 id_S2 id_S3 id_S4 mean_id sd_id
    ##    <chr>    <chr>     <chr>  <int> <int> <int> <int>   <dbl> <dbl>
    ##  1 23ms_40w 1         40 SPD  1624  1613  1630    NA   1622.  8.62
    ##  2 23ms_40w 2         40 SPD   908   903   910    NA    907   3.61
    ##  3 23ms_40w 3         40 SPD  2026  2014  2025    NA   2022.  6.66
    ##  4 30ms_30w 1         40 SPD  1890  1886  1819    NA   1865  39.9 
    ##  5 30ms_30w 2         40 SPD   799   794   780    NA    791   9.85
    ##  6 30ms_30w 3         40 SPD  1269  1267  1233    NA   1256. 20.2 
    ##  7 38ms_40w 1         20 SPD  1818  1809  1825    NA   1817.  8.02
    ##  8 38ms_40w 2         20 SPD  1229  1231  1230    NA   1230   1   
    ##  9 38ms_40w 3         20 SPD  3944  3940  3946    NA   3943.  3.06
    ## 10 50ms_30w 1         20 SPD  1681  1688  1704    NA   1691  11.8 
    ## 11 50ms_30w 2         20 SPD  1140  1142  1143    NA   1142.  1.53
    ## 12 50ms_30w 3         20 SPD  3489  3492  3496    NA   3492.  3.51
    ## 13 76ms_20w 1         20 SPD   829   825   821    NA    825   4   
    ## 14 76ms_20w 2         20 SPD   809   808   810    NA    809   1   
    ## 15 76ms_20w 3         20 SPD  3753  3757  3752    NA   3754   2.65

# Plot

For the violin plot the data from 20 SPD and 40 SPD were combined with
patchwork. So two different plots were created and the combined to
`figure_sup1`.

## Plot 40 SPD

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

``` r
CV_figure_sub1$y_errorbar <- CV_figure_sub1$mean_id
CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "1"] <-
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "1"] +
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "2"]+
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "3"]
CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "2"] <-
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "2"]+
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "3"]

figure_sub1_40SPD <- CV_figure_sub1 %>% 
  filter(SPD == "40 SPD") %>% 
  ggplot(aes(x = method, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7)+
  scale_x_discrete(labels = c("23ms_40w" = "  23ms
40w",
                              "30ms_30w" = "  30ms
30w"))+
  facet_grid(~SPD)+
  labs(x = "Method",
       y = "Identified peptides",
       fill = "CV < 20%")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.text.x = element_text(hjust = 0.62),
        axis.title.y = element_text(vjust = 2.2),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 22),
        legend.key.size = unit(0.9, "cm"),
        axis.title.x = element_text(vjust = -1),
        plot.margin = margin(b = 20,
                             l = 20))+
  scale_fill_manual(values = c("grey82", "#ff9090", "#ff3030"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 8000)) +
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure_sub1_40SPD
```

![](Sup.-Fig.-1-Optimizing-a-balance-between-Windowing-Schemes-and-ITs_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Plot 20 SPD

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

``` r
CV_figure_sub1$y_errorbar <- CV_figure_sub1$mean_id
CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "1"] <-
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "1"] +
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "2"]+
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "3"]
CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "2"] <-
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "2"]+
  CV_figure_sub1$y_errorbar[CV_figure_sub1$CV_filter == "3"]

figure_sub1_20SPD <- CV_figure_sub1 %>% 
  filter(SPD == "20 SPD") %>% 
  ggplot(aes(x = method, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7)+
  scale_x_discrete(labels = c("38ms_40w" = "  38ms
40w",
                              "50ms_30w" = "  50ms
30w",
                              "76ms_20w" = "  76ms
20w"))+
  facet_grid(~SPD)+
  labs(x = "Injection Time & Number of Windows",
       y = "Identidied peptides",
       fill = "CV")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.text.x = element_text(hjust = 0.62),
        axis.title.y = element_text(vjust = 2.2),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 22),
        legend.key.size = unit(0.9, "cm"),
        axis.title.x = element_text(vjust = -1),
        plot.margin = margin(b = 20,
                             l = 20))+
  scale_fill_manual(values = c("grey82", "#ff9090", "#ff3030"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 8000)) +
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure_sub1_20SPD
```

![](Sup.-Fig.-1-Optimizing-a-balance-between-Windowing-Schemes-and-ITs_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

# Join Plot 40 SPD and Plot 20 SPD

Figure_sub1_20SPD and Figure_sub1_40SPD are combined using patchwork to
create figure-sub1.

``` r
limit <- c(min(c(layer_scales(figure_sub1_40SPD)$y$get_limits(),
                 layer_scales(figure_sub1_20SPD)$y$get_limits())),
           max(c(layer_scales(figure_sub1_40SPD)$y$get_limits(),
                 layer_scales(figure_sub1_20SPD)$y$get_limits()))+500)
figure_sub1 <- 
  (((figure_sub1_40SPD + 
        theme(axis.title.x = element_blank(),
        legend.position = 'none')
    + plot_layout(widths = c(2,3))) +
  (figure_sub1_20SPD +
    theme(axis.title.x = element_text(hjust = -0.5, vjust = -1),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
          axis.title.y = element_blank()))))
figure_sub1
```

![](Sup.-Fig.-1-Optimizing-a-balance-between-Windowing-Schemes-and-ITs_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Sup. Fig. is saved as PDF

``` r
cairo_pdf(filename = "Sup. Fig. 1 Optimizing a balance between Windowing Schemes and ITs.pdf", width = 20.83, height = 25)
figure_sub1
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
    ##  [1] ggh4x_0.2.1     reprex_2.0.1    ghibli_0.3.2    patchwork_1.1.1
    ##  [5] ggpubr_0.4.0    pheatmap_1.0.12 psych_2.2.5     forcats_0.5.1  
    ##  [9] stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4     readr_2.1.2    
    ## [13] tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lubridate_1.8.0    lattice_0.20-45    assertthat_0.2.1   digest_0.6.29     
    ##  [5] utf8_1.2.2         R6_2.5.1           cellranger_1.1.0   backports_1.4.1   
    ##  [9] evaluate_0.15      highr_0.9          httr_1.4.3         pillar_1.7.0      
    ## [13] rlang_1.0.2        readxl_1.4.0       rstudioapi_0.13    car_3.0-13        
    ## [17] rmarkdown_2.14     labeling_0.4.2     munsell_0.5.0      broom_0.8.0       
    ## [21] compiler_4.2.0     modelr_0.1.8       xfun_0.31          pkgconfig_2.0.3   
    ## [25] mnormt_2.0.2       tmvnsim_1.0-2      htmltools_0.5.2    tidyselect_1.1.2  
    ## [29] fansi_1.0.3        crayon_1.5.1       tzdb_0.3.0         dbplyr_2.1.1      
    ## [33] withr_2.5.0        prismatic_1.1.0    grid_4.2.0         nlme_3.1-157      
    ## [37] jsonlite_1.8.0     gtable_0.3.0       lifecycle_1.0.1    DBI_1.1.2         
    ## [41] magrittr_2.0.3     scales_1.2.0       carData_3.0-5      cli_3.3.0         
    ## [45] stringi_1.7.6      farver_2.1.0       ggsignif_0.6.3     fs_1.5.2          
    ## [49] xml2_1.3.3         ellipsis_0.3.2     generics_0.1.2     vctrs_0.4.1       
    ## [53] RColorBrewer_1.1-3 tools_4.2.0        glue_1.6.2         hms_1.1.1         
    ## [57] abind_1.4-5        parallel_4.2.0     fastmap_1.1.0      yaml_2.3.5        
    ## [61] colorspace_2.0-3   rstatix_0.7.0      rvest_1.0.2        knitr_1.39        
    ## [65] haven_2.5.0
