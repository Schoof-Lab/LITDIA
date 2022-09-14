Lukas R. Woltereck & Samuel Grégoire

Sup. Fig. 1 Optimizing a balance between Windowing Schemes and ITs

Comparison of the number of identified peptides on DIA-LIT-based method
on Normal scanning mode for WhisperTM 20 SPD and Rapid scanning mode for
WhisperTM 40 SPD from 1 ng of tryptic HeLa lysate with different numbers
of windows and IT at fixed cycle time.

# Load packages

``` r
library("tidyverse")
library("patchwork")
```

# Load data

The path to every folder used is stored to keep everything compact. LIT
data is stored in the folder located in lit\_path. To reproduce analysis
the path must be changed to corresponding path on the local computer.

``` r
#lit_path <- "C:/Data/LIT/"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in RStudio with the function
`read.csv2()`.

## 40 SPD

``` r
#38ms 34w
LITDIA_Rapid_40SPD_23ms_40w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_auto_34w.csv"))
#30ms 30w
LITDIA_Rapid_40SPD_30ms_30w_1 <- read.csv2(paste0(lit_path, "20220509_LIT_Rapid_DIA_1ng_40SPD_whisper100_1CV_30ms_30w.csv"))
```

## 20 SPD

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
Spectronaut does imputation by default, but we want to keep missing
values. Imputed values are replaced by *NA*. Peptides missing in more
than one replicate are highlighted as *TRUE*.

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

## 40 SPD

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

## 20 SPD

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
throughput (samples per day) are included afterwards.

## 40 SPD

``` r
LITDIA_Rapid_40SPD_23ms_40w_1_CV <- compute_CV(LITDIA_Rapid_40SPD_23ms_40w_1)
LITDIA_Rapid_40SPD_23ms_40w_1_CV$method <- "23ms_40w"
LITDIA_Rapid_40SPD_23ms_40w_1_CV$SPD <- "40 SPD"
LITDIA_Rapid_40SPD_30ms_30w_1_CV <- compute_CV(LITDIA_Rapid_40SPD_30ms_30w_1)
LITDIA_Rapid_40SPD_30ms_30w_1_CV$method <- "30ms_30w"
LITDIA_Rapid_40SPD_30ms_30w_1_CV$SPD <- "40 SPD"
```

## 20 SPD

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
CV_figure_sup1 <- 
  rbind(LITDIA_Rapid_40SPD_23ms_40w_1_CV, LITDIA_Rapid_40SPD_30ms_30w_1_CV, 
        LITDIA_Normal_20SPD_38ms_40w_1_CV, LITDIA_Normal_20SPD_50ms_30w_1_CV,
        LITDIA_Normal_20SPD_76ms_20w_1_CV)

CV_figure_sup1$method %>% table()
```

    ## .
    ## 23ms_40w 30ms_30w 38ms_40w 50ms_30w 76ms_20w 
    ##     4594     3997     7023     6365     5406

``` r
CV_figure_sup1$CV_filter <- NA
```

A new column called “CV\_filter” is added to identify peptides with a CV
over 0.15 (1), CV between 0.1 and 0.15 (2) and CV under 0.1 (3).

``` r
CV_figure_sup1$CV_filter[CV_figure_sup1$CV < 0.1] <- "3"
CV_figure_sup1$CV_filter[CV_figure_sup1$CV >= 0.1 & CV_figure_sup1$CV <= 0.15] <- "2"
CV_figure_sup1$CV_filter[CV_figure_sup1$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out and
the number of identified peptides for each replicate and each CV\_filter
is calculated by summing the number of non-missing peptides. The mean
and the standard deviation of the replicates are also calculated for
each method.

``` r
CV_figure_sup1 <- CV_figure_sup1 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(method, CV_filter, SPD) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
```

    ## `summarise()` has grouped output by 'method', 'CV_filter'. You can override
    ## using the `.groups` argument.

``` r
CV_figure_sup1$id_S4[CV_figure_sup1$id_S4 == 0] <- NA
CV_figure_sup1$mean_id <- rowMeans(CV_figure_sup1[4:7], na.rm = TRUE)
CV_figure_sup1$sd_id <- apply(CV_figure_sup1[4:7], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

# Plot

## Plot 40 SPD

The mean and the standard deviation of identified peptides for each
CV\_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the “stacked” aspect of the col plot, error
bars y position is adjusted by using the y\_errorbar variable instead of
the mean\_id.

``` r
CV_figure_sup1$y_errorbar <- CV_figure_sup1$mean_id
CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "1"] <-
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "1"] +
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "2"]+
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "3"]
CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "2"] <-
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "2"]+
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "3"]

figure_sup1_40SPD <- CV_figure_sup1 %>% 
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
```

## Plot 20 SPD

``` r
CV_figure_sup1$y_errorbar <- CV_figure_sup1$mean_id
CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "1"] <-
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "1"] +
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "2"]+
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "3"]
CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "2"] <-
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "2"]+
  CV_figure_sup1$y_errorbar[CV_figure_sup1$CV_filter == "3"]

figure_sup1_20SPD <- CV_figure_sup1 %>% 
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
figure_sup1_20SPD
```

![](fig_sup1_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

# Join Plot 40 SPD and Plot 20 SPD

figure\_sup1\_20SPD and figure\_sup1\_40SPD are combined using patchwork
to create figure\_sup1.

``` r
limit <- c(min(c(layer_scales(figure_sup1_40SPD)$y$get_limits(),
                 layer_scales(figure_sup1_20SPD)$y$get_limits())),
           max(c(layer_scales(figure_sup1_40SPD)$y$get_limits(),
                 layer_scales(figure_sup1_20SPD)$y$get_limits()))+500)
figure_sup1 <- 
  (((figure_sup1_40SPD + 
        theme(axis.title.x = element_blank(),
        legend.position = 'none')
    + plot_layout(widths = c(2,3))) +
  (figure_sup1_20SPD +
    theme(axis.title.x = element_text(hjust = -0.5, vjust = -1),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
          axis.title.y = element_blank()))))
figure_sup1
```

![](fig_sup1_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
#save as pdf
#cairo_pdf(filename = "Sup. Fig. 1 Optimizing a balance between Windowing Schemes and ITs.pdf", width = 20.83, height = 25)
#figure_sup1
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
