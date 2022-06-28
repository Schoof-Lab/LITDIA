Fig. 5 Improvement of the Number of Points Across Chromatographic Peaks
================

Distribution of numbers of points across a peak by using Whisper100TM 20
SPD with DIA-LIT-based methods on Rapid, Normal, and Enhanced scanning
modes and Whisper100TM 40 SPD with DIA-LIT-based methods on Turbo,
Rapid, and Normal scanning modes on 1, 5, and 10 ng input material.

# Load packages

``` r
library("tidyverse")
library("patchwork")
#set bw theme as default theme
theme_set(theme_bw())
```

The path to every folder used is stored to keep everything compact. PPP
data is stored in the folder located in ppp_path. To reproduce analysis
the path must be changed to corresponding path on the local computer.

``` r
ppp_path <- "C://data/ppp"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in RStudio with the function
`read.csv2()`.

# Figure A

Distribution of numbers of points across a peak by using Whisper100TM 20
SPD with DIA-LIT-based methods on Rapid, Normal, and Enhanced scanning
modes and Whisper100TM 40 SPD with DIA-LIT-based methods on Turbo,
Rapid, and Normal scanning modes on 1, 5, and 10 ng input material.

## Load data A

The data containing the number of points across the peak are exported
separately from Spectronaut as tsv files. They are loaded in RStudio
with the function `read.csv2()`. Annotations regarding the origin file
are included.

1 ng

``` r
pp_Turbo_40SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Turbo_DIA_1ng_40SPD_1CV.tsv"), sep = "/")
pp_Turbo_40SPD_1$indiv <- "Turbo_40SPD_1_"
pp_Rapid_40SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_1ng_40SPD_1CV.tsv"), sep = "/")
pp_Rapid_40SPD_1$indiv <- "Rapid_40SPD_1_"
pp_Normal_40SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_40SPD_1CV.tsv"), sep = "/")
pp_Normal_40SPD_1$indiv <- "Normal_40SPD_1_"
pp_Rapid_20SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_1ng_20SPD_1CV.tsv"), sep = "/")
pp_Rapid_20SPD_1$indiv <- "Rapid_20SPD_1_"
pp_Normal_20SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_20SPD_1CV.tsv"), sep = "/")
pp_Normal_20SPD_1$indiv <- "Normal_20SPD_1_"
pp_Enhanced_20SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Enhanced_DIA_1ng_20SPD_1CV.tsv"), sep = "/")
pp_Enhanced_20SPD_1$indiv <- "Enhanced_20SPD_1_"
```

5 ng

``` r
pp_Turbo_40SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Turbo_DIA_5ng_40SPD_1CV.tsv"), sep = "/")
pp_Turbo_40SPD_5$indiv <- "Turbo_40SPD_5_"
pp_Rapid_40SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_5ng_40SPD_1CV.tsv"), sep = "/")
pp_Rapid_40SPD_5$indiv <- "Rapid_40SPD_5_"
pp_Normal_40SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_5ng_40SPD_1CV.tsv"), sep = "/")
pp_Normal_40SPD_5$indiv <- "Normal_40SPD_5_"
pp_Rapid_20SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_5ng_20SPD_1CV.tsv"), sep = "/")
pp_Rapid_20SPD_5$indiv <- "Rapid_20SPD_5_"
pp_Normal_20SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_5ng_20SPD_1CV.tsv"), sep = "/")
pp_Normal_20SPD_5$indiv <- "Normal_20SPD_5_"
pp_Enhanced_20SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Enhanced_DIA_5ng_20SPD_1CV.tsv"), sep = "/")
pp_Enhanced_20SPD_5$indiv <- "Enhanced_20SPD_5_"
```

10 ng

``` r
pp_Turbo_40SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Turbo_DIA_10ng_40SPD_1CV.tsv"), sep = "/")
pp_Turbo_40SPD_10$indiv <- "Turbo_40SPD_10_"
pp_Rapid_40SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_10ng_40SPD_1CV.tsv"), sep = "/")
pp_Rapid_40SPD_10$indiv <- "Rapid_40SPD_10_"
pp_Normal_40SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_10ng_40SPD_1CV.tsv"), sep = "/")
pp_Normal_40SPD_10$indiv <- "Normal_40SPD_10_"
pp_Rapid_20SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_10ng_20SPD_1CV.tsv"), sep = "/")
pp_Rapid_20SPD_10$indiv <- "Rapid_20SPD_10_"
pp_Normal_20SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_10ng_20SPD_1CV.tsv"), sep = "/")
pp_Normal_20SPD_10$indiv <- "Normal_20SPD_10_"
pp_Enhanced_20SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Enhanced_DIA_10ng_20SPD_1CV.tsv"), sep = "/")
pp_Enhanced_20SPD_10$indiv <- "Enhanced_20SPD_10_"
```

## Colnames A

Columns are renamed.

``` r
colnames(pp_Turbo_40SPD_1) <- colnames(pp_Turbo_40SPD_5) <- 
colnames(pp_Turbo_40SPD_10) <-  colnames(pp_Rapid_40SPD_1) <- 
colnames(pp_Rapid_40SPD_5) <- colnames(pp_Rapid_40SPD_10) <-  
colnames(pp_Normal_40SPD_1) <- colnames(pp_Normal_40SPD_5) <- 
colnames(pp_Normal_40SPD_10) <- colnames(pp_Rapid_20SPD_1) <- 
colnames(pp_Rapid_20SPD_5) <- colnames(pp_Rapid_20SPD_10) <-  
colnames(pp_Normal_20SPD_1) <- colnames(pp_Normal_20SPD_5) <-
colnames(pp_Normal_20SPD_10) <-  colnames(pp_Enhanced_20SPD_1) <- colnames(pp_Enhanced_20SPD_5) <- colnames(pp_Enhanced_20SPD_10) <- 
c("Points_per_Peaks", "indiv")
```

## Join A

Data frames containing PPP are combined with the function `rbind()`.

``` r
ppp <- rbind(pp_Turbo_40SPD_1, pp_Turbo_40SPD_5, pp_Turbo_40SPD_10, 
            pp_Rapid_40SPD_1, pp_Rapid_40SPD_5, pp_Rapid_40SPD_10, 
            pp_Normal_40SPD_1, pp_Normal_40SPD_5, pp_Normal_40SPD_10,
            pp_Rapid_20SPD_1, pp_Rapid_20SPD_5, pp_Rapid_20SPD_10,
            pp_Normal_20SPD_1, pp_Normal_20SPD_5, pp_Normal_20SPD_10, 
            pp_Enhanced_20SPD_1, pp_Enhanced_20SPD_5, pp_Enhanced_20SPD_10)
```

The number of samples per day, the injection time and the scanning mode
are added to the data frame and the function `factor()` is used to set
the order.

``` r
ppp$input <- NA
ppp$input[grepl("1_", ppp$indiv)] <- "1 ng"
ppp$input[grepl("5_", ppp$indiv)] <- "5 ng"
ppp$input[grepl("10_", ppp$indiv)] <- "10 ng"
ppp$input <- factor(ppp$input, levels = unique(ppp$input))
ppp$IT <- NA
ppp$IT[grepl("20SPD", ppp$indiv)] <- "20 SPD"
ppp$IT[grepl("40SPD", ppp$indiv)] <- "40 SPD"
ppp$IT <- factor(ppp$IT, levels = c("40 SPD", "20 SPD"))
ppp$method <- NA
ppp$method[grepl("Turbo", ppp$indiv)] <- "Turbo"
ppp$method[grepl("Rapid", ppp$indiv)] <- "Rapid"
ppp$method[grepl("Normal", ppp$indiv)] <- "Normal"
ppp$method[grepl("Enhanced", ppp$indiv)] <- "Enhanced"
ppp$method <- factor(ppp$method, levels = c("Turbo", "Rapid", "Normal", "Enhanced"))
```

## Plot A

PPP distribution from 20 SPD and 40 SPD are plotted separately using
violin plots in `ggplot2` and then combined together with patchwork. The
min and max value for each plot are stored in “limit” so both plot can
have the same scale while still including every values.

``` r
# 40 SPD
palette1 <-c("#C8F5FF", "#75CFFF", "#B3E4FF") 
plot_40SPD <- ppp %>%
  filter(IT == "40 SPD") %>% 
  ggplot(aes(x = input, y = Points_per_Peaks,
             fill = as.factor(indiv)))+
  scale_y_log10()+
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values = rep(palette1,3)) +
  facet_grid(IT ~ method)+
  geom_hline(yintercept = 6,
             linetype = "dashed") +
  labs(title = "A")+
  theme(legend.position = "none",
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.margin = margin(t = 20,  
                             r = 50,
                             b = 20,
                             l = 10),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.1))
# 20 SPD
palette2 <- c("#00E5E6", "#009999", "#00CCCC")
plot_20SPD <- ppp %>%
  filter(IT == "20 SPD") %>% 
  ggplot(aes(x = input, y = Points_per_Peaks,
             fill = as.factor(indiv)))+
  geom_violin(draw_quantiles = 0.5)+
  scale_y_log10() +
  scale_fill_manual(values = rep(palette2,3)) +
  facet_grid(IT ~ method) +
  geom_hline(yintercept = 6,
             linetype = "dashed") +
  labs(x = "Input",
       y = "Points across a peak") +
  theme(legend.position = "none",
        axis.title = element_text(size = 27),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.margin = margin(t = 15,  
                             r = 50,
                             b = 10,
                             l = 10))
limit <- 10^c(min(c(layer_scales(plot_40SPD)$y$get_limits(), 
                   layer_scales(plot_20SPD)$y$get_limits())),
              max(c(layer_scales(plot_40SPD)$y$get_limits(),
                   layer_scales(plot_20SPD)$y$get_limits())))
```

Plot_40SPD and Plot_20SPD are combined using patchwork to create figure
5a.

``` r
figure5a <- (plot_40SPD + 
    scale_y_log10(limits = limit,
                   breaks = c(3, 6, 10, 30))+
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))/
  (plot_20SPD + 
     scale_y_log10(limits = limit,
                   breaks = c(3, 6, 10, 30))+
     theme(axis.title.y = element_text(hjust = -1.8, vjust = 1)))
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.
    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

``` r
figure5a
```

![](Fig.-5-Improvement-of-the-Number-of-Points-Across-Chromatographic-Peaks_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Figure B

Comparison of the number of precursors with a number of points across a
peak equal to or greater than 6 on DIA-LIT-based method on Turbo, Rapid,
and Normal scanning mode with different input material (1 ng, 5 ng, and
10 ng of tryptic HeLa lysate) on WhisperTM 100 40 SPD and on Rapid,
Normal, and Enhanced scanning mode with different input material (1 ng,
5 ng, and 10 ng of tryptic HeLa lysate) on WhisperTM 100 20 SPD.

To create a stacked barplot, the data frame is separated in two
different groups. `temp2` includes all peptides with at least six points
across the peak. `temp1` all the others. Both are combined in `temp`
with `rbind()`.

``` r
ppp$acc <- ppp$Points_per_Peaks >= 6
temp1 <- ppp %>% 
  group_by(indiv, IT, input, method) %>% 
  summarise(n_prec = sum(!is.na(Points_per_Peaks)))
```

    ## `summarise()` has grouped output by 'indiv', 'IT', 'input'. You can override
    ## using the `.groups` argument.

``` r
temp2 <- ppp %>% 
  group_by(indiv, IT, input, method) %>% 
  summarise(n_acc = sum(acc))
```

    ## `summarise()` has grouped output by 'indiv', 'IT', 'input'. You can override
    ## using the `.groups` argument.

``` r
temp1$Precursors <- "< 6"
temp2$Precursors <- "\U2265 6"
colnames(temp1)[5] <- "n"
colnames(temp2)[5] <- "n"
temp1$n <- temp1$n - temp2$n
temp <- rbind(temp2, temp1)
temp$Precursors <- factor(temp$Precursors, 
                         levels = c("< 6", "\U2265 6"))
```

PPP numbers from 20 SPD and 40 SPD are plotted separately using stacked
col plot in `ggplot2` and then combined together with patchwork. n
column contains the number of precursors for all the replicates, this
number is divided by 4 obtain the mean.

``` r
id_pr_40SPD <- temp %>% 
  filter(IT == "40 SPD" ) %>% 
  ggplot(aes(x = method, y = n/4, fill = Precursors, group = method)) +
  geom_col() +
  theme_bw() +
  facet_grid(IT ~ input) +
  labs(x = "Scanning mode",
       y = "Precursors number",
       title = "B",
       fill = "Points across a peak")+
  scale_fill_manual(values = c("grey82", "#ff3030"))+
  theme(plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.095),
        axis.title = element_blank(),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.title = element_text(size = 27, vjust = 1.3),
        legend.text = element_text(size = 22),
        legend.position = "top")
id_pr_20SPD <- temp %>% 
  filter(IT == "20 SPD" ) %>% 
  ggplot(aes(x = method , y = n/4, fill = Precursors, group = method)) +
  geom_col() +
  theme_bw() +
  facet_grid(IT ~ input) +
  labs(x = "Scanning mode",
       y = "Precursors number")+
  scale_fill_manual(values = c("grey82", "#ff3030"))+
  theme(axis.title = element_text(size = 27),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 27),
        strip.background.x = element_blank(),
        strip.background = element_rect(fill = "#E5E4F7"),
        strip.text.x = element_blank(),
        axis.title.x = element_text(vjust = -1),
        legend.position = "none")
limit <- c(min(c(layer_scales(id_pr_40SPD)$y$get_limits(),
                 layer_scales(id_pr_20SPD)$y$get_limits())),
           max(c(layer_scales(id_pr_40SPD)$y$get_limits(),
                 layer_scales(id_pr_20SPD)$y$get_limits()))+500)
figure5b <-
 (id_pr_40SPD +
    theme(axis.title = element_blank())+
    scale_y_continuous(limits = limit))/
 (id_pr_20SPD +
    theme(axis.title.y = element_text(hjust = -2.15,
                                      vjust = 3.5))+
    scale_y_continuous(limits = limit))
figure5b
```

![](Fig.-5-Improvement-of-the-Number-of-Points-Across-Chromatographic-Peaks_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# Plot A and B

Plot 5A and 5B are combined using patchwork to create figure 5 and save
it as pdf.

``` r
figure5 <- wrap_plots(figure5a,figure5b, nrow = 2)
figure5 
```

![](Fig.-5-Improvement-of-the-Number-of-Points-Across-Chromatographic-Peaks_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
cairo_pdf(filename = "Fig. 5 Improvement of the Number of Points Across Chromatographic Peaks.pdf", width = 16.67, height = 18.75,)
figure5 
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
    ## [37] digest_0.6.29    stringi_1.7.6    grid_4.2.0       cli_3.3.0       
    ## [41] tools_4.2.0      magrittr_2.0.3   crayon_1.5.1     pkgconfig_2.0.3 
    ## [45] ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0 
    ## [49] assertthat_0.2.1 rmarkdown_2.14   httr_1.4.3       rstudioapi_0.13 
    ## [53] R6_2.5.1         compiler_4.2.0
