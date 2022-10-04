Fig. 7 Improvement of the Number of Points Across Chromatographic Peaks.
================
Samuel Grégoire & Lukas R. Woltereck

Fig. 7 Improvement of the Number of Points Across Chromatographic Peaks.

Distribution of numbers of points across a peak by using Whisper100TM 20
SPD with LIT-DIA-based methods on Rapid, Normal, and Enhanced scanning
modes and Whisper100TM 40 SPD with LIT-DIA-based methods on Turbo,
Rapid, and Normal scanning modes on 1, 5, and 10 ng input material. B.
Comparison of the number of precursors with a number of points across a
peak equal to or greater than 6 on the LIT-DIA-based method on Turbo,
Rapid, and Normal scanning modes with different input material (1 ng, 5
ng, and 10 ng of tryptic HeLa digest) on WhisperTM 100 40 SPD and on
Rapid, Normal, and Enhanced scanning mode with different input material
(1 ng, 5 ng, and 10 ng of tryptic HeLa digest) on WhisperTM 100 20 SPD.

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
#ppp_path <- "C:/Data/PPP/"
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

``` r
#1 ng
ppp_Turbo_40SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Turbo_DIA_1ng_40SPD_1CV.tsv"), sep = "/")
ppp_Turbo_40SPD_1$indiv <- "Turbo_40SPD_1_"
ppp_Rapid_40SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_1ng_40SPD_1CV.tsv"), sep = "/")
ppp_Rapid_40SPD_1$indiv <- "Rapid_40SPD_1_"
ppp_Normal_40SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_40SPD_1CV.tsv"), sep = "/")
ppp_Normal_40SPD_1$indiv <- "Normal_40SPD_1_"
ppp_Rapid_20SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_1ng_20SPD_1CV.tsv"), sep = "/")
ppp_Rapid_20SPD_1$indiv <- "Rapid_20SPD_1_"
ppp_Normal_20SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_20SPD_1CV.tsv"), sep = "/")
ppp_Normal_20SPD_1$indiv <- "Normal_20SPD_1_"
ppp_Enhanced_20SPD_1 <- read.csv2(paste0(ppp_path, "PPP_LIT_Enhanced_DIA_1ng_20SPD_1CV.tsv"), sep = "/")
ppp_Enhanced_20SPD_1$indiv <- "Enhanced_20SPD_1_"

#5 ng
ppp_Turbo_40SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Turbo_DIA_5ng_40SPD_1CV.tsv"), sep = "/")
ppp_Turbo_40SPD_5$indiv <- "Turbo_40SPD_5_"
ppp_Rapid_40SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_5ng_40SPD_1CV.tsv"), sep = "/")
ppp_Rapid_40SPD_5$indiv <- "Rapid_40SPD_5_"
ppp_Normal_40SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_5ng_40SPD_1CV.tsv"), sep = "/")
ppp_Normal_40SPD_5$indiv <- "Normal_40SPD_5_"
ppp_Rapid_20SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_5ng_20SPD_1CV.tsv"), sep = "/")
ppp_Rapid_20SPD_5$indiv <- "Rapid_20SPD_5_"
ppp_Normal_20SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_5ng_20SPD_1CV.tsv"), sep = "/")
ppp_Normal_20SPD_5$indiv <- "Normal_20SPD_5_"
ppp_Enhanced_20SPD_5 <- read.csv2(paste0(ppp_path, "PPP_LIT_Enhanced_DIA_5ng_20SPD_1CV.tsv"), sep = "/")
ppp_Enhanced_20SPD_5$indiv <- "Enhanced_20SPD_5_"

#10 ng
ppp_Turbo_40SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Turbo_DIA_10ng_40SPD_1CV.tsv"), sep = "/")
ppp_Turbo_40SPD_10$indiv <- "Turbo_40SPD_10_"
ppp_Rapid_40SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_10ng_40SPD_1CV.tsv"), sep = "/")
ppp_Rapid_40SPD_10$indiv <- "Rapid_40SPD_10_"
ppp_Normal_40SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_10ng_40SPD_1CV.tsv"), sep = "/")
ppp_Normal_40SPD_10$indiv <- "Normal_40SPD_10_"
ppp_Rapid_20SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Rapid_DIA_10ng_20SPD_1CV.tsv"), sep = "/")
ppp_Rapid_20SPD_10$indiv <- "Rapid_20SPD_10_"
ppp_Normal_20SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_10ng_20SPD_1CV.tsv"), sep = "/")
ppp_Normal_20SPD_10$indiv <- "Normal_20SPD_10_"
ppp_Enhanced_20SPD_10 <- read.csv2(paste0(ppp_path, "PPP_LIT_Enhanced_DIA_10ng_20SPD_1CV.tsv"), sep = "/")
ppp_Enhanced_20SPD_10$indiv <- "Enhanced_20SPD_10_"
```

## Colnames A

Columns are renamed.

``` r
colnames(ppp_Turbo_40SPD_1) <- colnames(ppp_Turbo_40SPD_5) <- 
colnames(ppp_Turbo_40SPD_10) <-  colnames(ppp_Rapid_40SPD_1) <- 
colnames(ppp_Rapid_40SPD_5) <- colnames(ppp_Rapid_40SPD_10) <-  
colnames(ppp_Normal_40SPD_1) <- colnames(ppp_Normal_40SPD_5) <- 
colnames(ppp_Normal_40SPD_10) <- colnames(ppp_Rapid_20SPD_1) <- 
colnames(ppp_Rapid_20SPD_5) <- colnames(ppp_Rapid_20SPD_10) <-  
colnames(ppp_Normal_20SPD_1) <- colnames(ppp_Normal_20SPD_5) <-
colnames(ppp_Normal_20SPD_10) <-  colnames(ppp_Enhanced_20SPD_1) <- 
  colnames(ppp_Enhanced_20SPD_5) <- colnames(ppp_Enhanced_20SPD_10) <- 
  c("Points_per_Peaks", "indiv")
```

## Join A

Data frames containing PPP are combined with the function `rbind()`.

``` r
ppp <- rbind(ppp_Turbo_40SPD_1, ppp_Turbo_40SPD_5, ppp_Turbo_40SPD_10, 
            ppp_Rapid_40SPD_1, ppp_Rapid_40SPD_5, ppp_Rapid_40SPD_10, 
            ppp_Normal_40SPD_1, ppp_Normal_40SPD_5, ppp_Normal_40SPD_10,
            ppp_Rapid_20SPD_1, ppp_Rapid_20SPD_5, ppp_Rapid_20SPD_10,
            ppp_Normal_20SPD_1, ppp_Normal_20SPD_5, ppp_Normal_20SPD_10, 
            ppp_Enhanced_20SPD_1, ppp_Enhanced_20SPD_5, ppp_Enhanced_20SPD_10)
```

The throughput (samples per day), the injection time and the scanning
mode are added to the data frame and the function `factor()` is used to
set the order.

``` r
ppp$input <- NA
ppp$input[grepl("1_", ppp$indiv)] <- 1
ppp$input[grepl("5_", ppp$indiv)] <- 5
ppp$input[grepl("10_", ppp$indiv)] <- 10
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

PPP distributions from 20 SPD and 40 SPD are plotted separately using
violin plots in `ggplot2` and then combined together with patchwork. The
min and max value for each plot are stored in “limit” so both plot can
have the same scale while still including every values.

``` r
# 40 SPD
plot_40SPD <- ppp %>%
  filter(IT == "40 SPD") %>% 
  ggplot(aes(x = input, y = Points_per_Peaks))+
  scale_y_log10()+
  geom_violin(draw_quantiles = 0.5,
              fill = "#D1E5EB") +
  facet_grid(IT ~ method)+
  geom_hline(yintercept = 6,
             linetype = "dashed") +
  labs(title = "A")+
  theme(legend.position = "none",
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.1))
# 20 SPD
plot_20SPD <- ppp %>%
  filter(IT == "20 SPD") %>% 
  ggplot(aes(x = input, y = Points_per_Peaks,
             fill = "#D1E5EB"))+
  geom_violin(draw_quantiles = 0.5,
              fill = "#D1E5EB")+
  scale_y_log10() +
  facet_grid(IT ~ method) +
  geom_hline(yintercept = 6,
             linetype = "dashed") +
  labs(x = "Input (ng)",
       y = "Points across a chromatographic peak") +
  theme(legend.position = "none",
        axis.title = element_text(size = 27),
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"))
limit <- 10^c(min(c(layer_scales(plot_40SPD)$y$get_limits(), 
                   layer_scales(plot_20SPD)$y$get_limits())),
              max(c(layer_scales(plot_40SPD)$y$get_limits(),
                   layer_scales(plot_20SPD)$y$get_limits())))
```

Plot_40SPD and Plot_20SPD are combined using patchwork to create figure
7a.

``` r
figure7a <- (plot_40SPD + 
    scale_y_log10(limits = limit,
                   breaks = c(3, 6, 10, 30))+
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))/
  (plot_20SPD + 
     scale_y_log10(limits = limit,
                   breaks = c(3, 6, 10, 30))+
     theme(axis.title.y = element_text(hjust = -0.05)))
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.
    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

``` r
figure7a
```

![](Figure7_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Figure B

Comparison of the number of precursors with a number of points across a
peak equal to or greater than 6 on DIA-LIT-based method on Turbo, Rapid,
and Normal scanning mode with different input material (1 ng, 5 ng, and
10 ng of tryptic HeLa lysate) on WhisperTM 100 40 SPD and on Rapid,
Normal, and Enhanced scanning mode with different input material (1 ng,
5 ng, and 10 ng of tryptic HeLa lysate) on WhisperTM 100 20 SPD.

To create a stacked barplot, the data frame is separated in two
different groups. `ppp_acc` includes all precursors with at least six
points across the peak and `ppp_tot` all precursors with ppp
information. Both are combined in `ppp_join` with `rbind()`.

``` r
ppp$acc <- ppp$Points_per_Peaks >= 6
ppp_tot <- ppp %>% 
  group_by(indiv, IT, input, method) %>% 
  summarise(n_prec = sum(!is.na(Points_per_Peaks)))
```

    ## `summarise()` has grouped output by 'indiv', 'IT', 'input'. You can override
    ## using the `.groups` argument.

``` r
ppp_acc <- ppp %>% 
  group_by(indiv, IT, input, method) %>% 
  summarise(n_acc = sum(acc))
```

    ## `summarise()` has grouped output by 'indiv', 'IT', 'input'. You can override
    ## using the `.groups` argument.

``` r
ppp_tot$Precursors <- "< 6"
ppp_acc$Precursors <- "\U2265 6"
colnames(ppp_tot)[5] <- "n"
colnames(ppp_acc)[5] <- "n"
ppp_tot$n <- ppp_tot$n - ppp_acc$n
ppp_join <- rbind(ppp_acc, ppp_tot)
ppp_join$Precursors <- factor(ppp_join$Precursors, 
                         levels = c("< 6", "\U2265 6"))
```

To visualize `ppp_join` a stacked col plot was created, using `ggplot2`.

PPP numbers from 20 SPD and 40 SPD are plotted separately using stacked
col plot in `ggplot2` and then combined together with patchwork. n
column contains the number of precursors for all the replicates, this
number is divided by 4 to obtain the mean.

``` r
id_pr_40SPD <- ppp_join %>% 
  filter(IT == "40 SPD" ) %>% 
  ggplot(aes(x = method, y = n/4, fill = Precursors, group = method)) +
  geom_col() +
  theme_bw() +
  facet_grid(IT ~ input,
             labeller = as_labeller(c("40 SPD" = "40 SPD",
                                      "1" = "1 ng",
                                      "5" = "5 ng",
                                      "10" = "10 ng"))) +
  labs(x = "Scanning mode",
       y = "Precursors number", 
       title = "B",
       fill = "Points across a chromatographic peak")+
  scale_fill_manual(values = c("grey82", "#ff713e"))+
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

id_pr_20SPD <- ppp_join %>% 
  filter(IT == "20 SPD" ) %>% 
  ggplot(aes(x = method , y = n/4, fill = Precursors, group = method)) +
  geom_col() +
  theme_bw() +
  facet_grid(IT ~ input) +
  labs(x = "Scanning mode",
       y = "Precursors number")+
  scale_fill_manual(values = c("grey82", "#ff713e"))+
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
figure7b <-
 (id_pr_40SPD +
    theme(axis.title = element_blank())+
    scale_y_continuous(limits = limit))/
 (id_pr_20SPD +
    theme(axis.title.y = element_text(hjust = -3.5,
                                      vjust = 3))+
    scale_y_continuous(limits = limit))
figure7b
```

![](Figure7_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Plot A and B

Plot 7A and 7B are combined using patchwork to create figure 7 and save
it as pdf.

``` r
figure7 <- wrap_plots(figure7a, figure7b, nrow = 2)
figure7 
```

![](Figure7_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
#save as pdf
#cairo_pdf(filename = "fig7.pdf", width = 16.67, height = 18.75)
#figure7 
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
