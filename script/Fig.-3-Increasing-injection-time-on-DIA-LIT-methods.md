Fig. 3 Increasing injection time on DIA LIT methods
================

# Load packages

``` r
library("tidyverse")
library("patchwork")
#set bw theme as default theme
theme_set(theme_bw())
```

The path to every folder used is stored to keep everything compact. LIT
data is stored in the folder located in lit_path. To reproduce analysis
the path must be changed to corresponding path on the local computer.

``` r
lit_path <- "C://data/lit"
ppp_path <- "C://data/ppp"
```

Data tables were exported from Spectronaut (version 15.5.211111.50606)
in csv files and need to be loaded in RStudio with the function
`read.csv2()`.

# Figure A

Comparison of the number of identified peptides on DIA-LIT-based method
on Normal scanning mode for WhisperTM 100 20 SPD from 1 ng of tryptic
HeLa lysate with different ITs at fixed 40 isolation windows. Identified
peptides with a coefficient of variation (CV) between 10% and 15% are
colored with light red and those with a CV below 10% with dark red. The
cycle times for the methods are indicated by their injection times: 38
ms 2.4 s, for 60 ms 3.28 s, for 80 ms 4.09 s, and for 100 ms 4.91 s.

# Load data A

``` r
LITDIA_Normal_20SPD_60ms <- read.csv2(paste0(lit_path, "20220215_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV45_60ms.csv"))
LITDIA_Normal_20SPD_80ms <- read.csv2(paste0(lit_path, "20220215_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV45_80ms.csv"))
LITDIA_Normal_20SPD_100ms <- read.csv2(paste0(lit_path, "20220215_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV45_100ms.csv"))
LITDIA_Normal_20SPD_38ms <- read.csv2(paste0(lit_path, "20220525_LIT_Normal_DIA_1ng_20SPD_whisper100_1CV45_auto.csv"))
```

## Clean data tables A

Imported tables contain a lot of information, which are not useful for
our analysis. The function `filter_pep()` only keeps the file names
(replicates), quantitative values for each file name, and the peptide
sequence. It also calculates the mean between the four replicates.
Spectronaut does some imputation by default, but we want to keep missing
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

``` r
LITDIA_Normal_20SPD_60ms <- filter_pep(LITDIA_Normal_20SPD_60ms)
```

    ## 
    ## Highlighted 94 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_20SPD_80ms <- filter_pep(LITDIA_Normal_20SPD_80ms)
```

    ## 
    ## Highlighted 73 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_20SPD_100ms <- filter_pep(LITDIA_Normal_20SPD_100ms)
```

    ## 
    ## Highlighted 98 peptide(s) found in less than 3 replicates

``` r
LITDIA_Normal_20SPD_38ms <- filter_pep(LITDIA_Normal_20SPD_38ms)
```

    ## 
    ## Highlighted 39 peptide(s) found in less than 3 replicates

## Number of identified peptides A

Empty vectors are created to store the number of identified peptides in
each replicates for each data-set.

``` r
id_LITDIA_Normal_20SPD_60ms <- id_LITDIA_Normal_20SPD_80ms <- 
  id_LITDIA_Normal_20SPD_100ms <- id_LITDIA_Normal_20SPD_38ms <-
  rep(NA, 4)
```

The numbers of identified peptides are counted by summing non-missing
values for each replicate and stored in the vectors created above.

``` r
for(i in 2:(ncol(LITDIA_Normal_20SPD_60ms)-2)){
  id_LITDIA_Normal_20SPD_60ms[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_60ms[, i]))}
for(i in 2:(ncol(LITDIA_Normal_20SPD_80ms)-2)){
  id_LITDIA_Normal_20SPD_80ms[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_80ms[, i]))}
for(i in 2:(ncol(LITDIA_Normal_20SPD_100ms)-2)){
  id_LITDIA_Normal_20SPD_100ms[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_100ms[, i]))}
for(i in 2:(ncol(LITDIA_Normal_20SPD_38ms)-2)){
  id_LITDIA_Normal_20SPD_38ms[i-1] <- sum(!is.na(LITDIA_Normal_20SPD_38ms[, i]))}
```

## Compute CV A

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

Information about the injection time are included afterwards.

``` r
LITDIA_Normal_20SPD_38ms_CV <- compute_CV(LITDIA_Normal_20SPD_38ms)
LITDIA_Normal_20SPD_38ms_CV$IT <- "38ms"
LITDIA_Normal_20SPD_60ms_CV <- compute_CV(LITDIA_Normal_20SPD_60ms)
LITDIA_Normal_20SPD_60ms_CV$IT <- "60ms"
LITDIA_Normal_20SPD_80ms_CV <- compute_CV(LITDIA_Normal_20SPD_80ms)
LITDIA_Normal_20SPD_80ms_CV$IT <- "80ms"
LITDIA_Normal_20SPD_100ms_CV <- compute_CV(LITDIA_Normal_20SPD_100ms)
LITDIA_Normal_20SPD_100ms_CV$IT <- "100ms"
```

## Join A

Data frames containing CV are combined with the function `rbind()`.

``` r
CV_figure3 <- 
  rbind(LITDIA_Normal_20SPD_38ms_CV, LITDIA_Normal_20SPD_60ms_CV, 
        LITDIA_Normal_20SPD_80ms_CV, LITDIA_Normal_20SPD_100ms_CV)
# Factorise to set order
CV_figure3$IT <- factor(CV_figure3$IT, 
                          levels = c("38ms", "60ms", "80ms", "100ms"))
```

A new column called ???CV_filter??? is added to identify peptides with a CV
over 0.15 (1), between 0.1 and 0.15 (2) and under 0.1 (3).

``` r
CV_figure3$CV_filter <- NA
CV_figure3$CV_filter[CV_figure3$CV < 0.1] <- "3"
CV_figure3$CV_filter[CV_figure3$CV >= 0.1 & CV_figure3$CV <= 0.15] <- "2"
CV_figure3$CV_filter[CV_figure3$CV > 0.15] <- "1"
```

Peptides without CV (found in just one replicate) are filtered out. The
number of identified peptides for each replicate and each CV_filter is
calculated by summing the number of non-missing peptides. The mean and
the standard deviation of the replicates are also calculated for each
method.

``` r
CV_figure3$SPD <- "20 SPD"
CV_figure3 <- CV_figure3 %>%
  filter(!is.na(CV_filter)) %>% 
  group_by(SPD, IT, CV_filter) %>% 
  summarise(id_S1 = sum(!is.na(S1)),
            id_S2 = sum(!is.na(S2)),
            id_S3 = sum(!is.na(S3)),
            id_S4 = sum(!is.na(S4)))
CV_figure3$mean_id <- rowMeans(CV_figure3[4:7], na.rm = TRUE)
CV_figure3$sd_id <- apply(CV_figure3[4:7], FUN =  sd, MARGIN = 1, na.rm = TRUE)
```

## Plot A

The mean and the standard deviation of identified peptides for each
CV_filter and method are plotted in a stacked col plot using `ggplot2`.
Error bars representing the standard deviation are displayed using
`geom_errorbar()`. Due to the ???stacked??? aspect of the col plot, error
bars y position is adjusted by using the y_errorbar variable instead of
the mean_id.

``` r
CV_figure3$y_errorbar <- CV_figure3$mean_id
CV_figure3$y_errorbar[CV_figure3$CV_filter == "1"] <-
  CV_figure3$y_errorbar[CV_figure3$CV_filter == "1"] +
  CV_figure3$y_errorbar[CV_figure3$CV_filter == "2"]+
  CV_figure3$y_errorbar[CV_figure3$CV_filter == "3"]
CV_figure3$y_errorbar[CV_figure3$CV_filter == "2"] <-
  CV_figure3$y_errorbar[CV_figure3$CV_filter == "2"]+
  CV_figure3$y_errorbar[CV_figure3$CV_filter == "3"]
figure3a <- CV_figure3 %>% 
  ggplot(aes(x = IT, y = mean_id, fill = CV_filter))+
  geom_col(width = 0.7) +
  facet_grid(~SPD) +
  theme_bw() +
  labs(x = "Injection time",
       y = "Identified peptides",
       fill = "CV",
       title = "A") +
  theme(plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.095),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 27),
        strip.background = element_rect(fill = "#E5E4F7"),
        axis.text.x = element_text(hjust = 0.62),
        axis.title.y = element_text(vjust = 2.2),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 22),
        legend.key.size = unit(0.9, "cm"),
        axis.title.x = element_text(vjust = -1),
        legend.position = "top",
        plot.margin = margin(b = 20,
                             l = 20))+
  scale_fill_manual(values = c("grey82", "#ff9090", "#ff3030"),
                    labels = c("\U2265 15%", "10 - 15%", "< 10%" ))+
  scale_y_continuous(limits = c(0, 8000)) +
  geom_errorbar(aes(ymin = y_errorbar - sd_id, 
                         ymax = y_errorbar + sd_id),
                    width = 0.2,
                    size = 1)
figure3a
```

![](Fig.-3-Increasing-injection-time-on-DIA-LIT-methods_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

# Figure B

Comparison of the number of precursors with a number of points across a
peak equal to or greater than 6 on DIA-LIT-based method on *Normal*
scanning mode from 1 ng input of tryptic HeLa lysate on
Whisper<sup>TM</sup> 100 20 SPD with different ITs (38 ms, 60 ms, 80 ms
and 100 ms) at fixed 40 isolation windows.

## Load B

The data containing the number of points across the peak are exported
separately from Spectronaut (version 15.5.211111.50606) as tsv file.
They are loaded in RStudio with the function `read.csv2()`.

``` r
ppp_Turbo_20SPD_1ng_38ms <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_20SPD_1CV_38ms.tsv"), sep = "/")
ppp_Turbo_20SPD_1ng_38ms$indiv <- "20SPD_1ng_38ms"
ppp_Turbo_20SPD_1ng_60ms <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_20SPD_1CV_60ms.tsv"), sep = "/")
ppp_Turbo_20SPD_1ng_60ms$indiv <- "20SPD_1ng_60ms"
ppp_Turbo_20SPD_1ng_80ms <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_20SPD_1CV_80ms.tsv"), sep = "/")
ppp_Turbo_20SPD_1ng_80ms$indiv <- "20SPD_1ng_80ms"
ppp_Turbo_20SPD_1ng_100ms <- read.csv2(paste0(ppp_path, "PPP_LIT_Normal_DIA_1ng_20SPD_1CV_100ms.tsv"), sep = "/")
ppp_Turbo_20SPD_1ng_100ms$indiv <- "20SPD_1ng_100ms"
```

## Join B

Data frames containing CV are combined with the function `rbind()`.

``` r
colnames(ppp_Turbo_20SPD_1ng_38ms) <- colnames(ppp_Turbo_20SPD_1ng_60ms) <- colnames(ppp_Turbo_20SPD_1ng_80ms) <- colnames(ppp_Turbo_20SPD_1ng_100ms) <-
  c("Points_per_Peaks", "indiv")
ppfig3 <- rbind(ppp_Turbo_20SPD_1ng_38ms, ppp_Turbo_20SPD_1ng_60ms, ppp_Turbo_20SPD_1ng_80ms,
            ppp_Turbo_20SPD_1ng_100ms)
```

The number of samples per day and the injection time are added to the
data frame and the function `factor()` is used to set the order.

``` r
ppfig3$spd <- NA
ppfig3$spd[grepl("20SPD", ppfig3$indiv)] <- "20 SPD"
ppfig3$spd[grepl("40SPD", ppfig3$indiv)] <- "40 SPD"
ppfig3$spd <- factor(ppfig3$spd, levels = c("40 SPD", "20 SPD"))
ppfig3$IT <- NA
ppfig3$IT[grepl("38ms", ppfig3$indiv)] <- "38 ms"
ppfig3$IT[grepl("60ms", ppfig3$indiv)] <- "60 ms"
ppfig3$IT[grepl("80ms", ppfig3$indiv)] <- "80 ms"
ppfig3$IT[grepl("100ms", ppfig3$indiv)] <- "100 ms"
ppfig3$IT <- factor(ppfig3$IT, levels = unique(ppfig3$IT))
```

## Plot B

To create a stacked barplot, the data frame is separated in two
different groups. `temp2` includes all peptides with at least six points
across the peak. `temp1` all the others. Both are combined in `temp`
with `rbind()`.

``` r
ppfig3$acc <- ppfig3$Points_per_Peaks >= 6
temp1 <- ppfig3 %>% 
  group_by(indiv, IT, spd) %>% 
  summarise(n_prec = sum(!is.na(Points_per_Peaks)))
temp2 <- ppfig3 %>% 
  group_by(indiv, IT, spd) %>% 
  summarise(n_acc = sum(acc))
temp1$Precursors <- "< 6"
temp2$Precursors <- "\U2265 6"
colnames(temp1)[4] <- "n"
colnames(temp2)[4] <- "n"
temp1$n <- temp1$n - temp2$n
temp <- rbind(temp2, temp1)
temp$Precursors <- factor(temp$Precursors, 
                         levels = c("< 6", "\U2265 6"))
```

To visualize `temp` a stacked col plot was created, using `ggplot2`.

``` r
figure3b <- temp %>% 
  ggplot(aes(x = IT, y = n/4, fill = Precursors, group = IT)) +
  geom_col(width = 0.7) +
  theme_bw() +
  facet_grid(~ spd)+
  labs(x = "Injection time",
       y = "Precursors number",
       title = "B",
       fill = "Points across a peak")+
  scale_fill_manual(values = c("grey82", "#ff3030"))+
  theme(plot.title = element_text(size = 32, face = "bold",
                                  hjust = -0.095),
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 2.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 22),
        legend.position = "top")
figure3b
```

![](Fig.-3-Increasing-injection-time-on-DIA-LIT-methods_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

# Join A and B

Plot 3A and 3B are combined using patchwork to create figure 3 and save
it as pdf.

``` r
figure3 <- 
(figure3a + theme(plot.margin = unit(c(0,18,0,0), "pt"))) +
(figure3b + theme(plot.margin = unit(c(0,0,0,18), "pt")))
figure3
```

![](Fig.-3-Increasing-injection-time-on-DIA-LIT-methods_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
cairo_pdf(filename = "Fig. 3 Increasing injection time on DIA LIT methods.pdf", width = 20.83, height = 7.81)
figure3
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
