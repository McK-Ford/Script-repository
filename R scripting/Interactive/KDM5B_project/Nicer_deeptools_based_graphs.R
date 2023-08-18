library(tidyverse) #can replace much of the above!!!

###########################
## Stranded profile data ##
###########################
#can handle groseq this way too...

format_deeptools_profile_tabs <- function(tab_file){
  tab <- read_table(tab_file, col_names = FALSE, col_types = cols(X1=col_skip()), skip=2) #imports a deeptools outFileNameData file
  reformat_tab <- tab %>%
    pivot_longer(cols=2:length(tab),
                 names_to = "Bins",
                 names_prefix="X",
                 values_to = "score",
                 values_drop_na = TRUE) %>% #pivots table into longform with bins as one column and values as another
    mutate(Bins=(as.numeric(Bins)-2), #because extra columns interfere with the bin naming (bin 1 is imported as X3 for example)
           .keep="unused") #need to fix the numbering by taking X off and subtracting 2
  return(reformat_tab)
}

H3K79me2 <- format_deeptools_profile_tabs("McKayla Ford-selected/G13.H3K79me2.WT.scaled.tab.txt")
H3K4me2 <- format_deeptools_profile_tabs("McKayla Ford-selected/G15.H3K4me2.WT.scaled.tab.txt")
H3K4me3 <- format_deeptools_profile_tabs("McKayla Ford-selected/G15.H3K4me3.WT.scaled.tab.txt")
Groseq_m <- format_deeptools_profile_tabs("McKayla Ford-selected/groseq_rep3_minus.scaled.tab.txt")
Groseq_p <- format_deeptools_profile_tabs("McKayla Ford-selected/groseq_rep3_plus.scaled.tab.txt")
R3_m <- format_deeptools_profile_tabs("McKayla Ford-selected/R3_rpm_groseq_minus.scaled.tab.txt")
R3_p <- format_deeptools_profile_tabs("McKayla Ford-selected/R3_rpm_groseq_plus.scaled.tab.txt")



m_long <- format_deeptools_profile_tabs("skew_10_minus.scaled.tab.txt")
p_long <- format_deeptools_profile_tabs("skew_10_plus.scaled.tab.txt")

#confirms I didn't mess up the +- sitch.
skew_bound <- rbind(p_long, m_long)
ggplot(skew_bound, aes(x=Bins, y=score, color=X2)) +
  theme_bw() + geom_line()
#################################################################################3
Groseq <- rbind(Groseq_m, Groseq_p)
ggplot(Groseq, aes(x=Bins, y=score, color=X2)) +
  theme_bw() + geom_line()
R3 <- rbind(R3_p, R3_m)
ggplot(R3, aes(x=Bins, y=score, color=X2)) +
  theme_bw() + geom_line() #okay, R3 danko mapped is clearly highly messed up, is it inverted because of the minus bed? Maybe.
#so ignore R3 for now.
##########################################################################################
#okay now for a proper graph
TSS_m <- m_long %>% filter(X2=="qnd_TSS_skew_minus.bed") %>% select(!X2)
TSS_p <- p_long %>% filter(X2=="qnd_TSS_skew_plus.bed") %>% select(!X2) #split into TSS / Dist
TSS <- TSS_m #bind plus and minus of a given class, here TSS, together.
TSS$score <- (TSS_m$score+TSS_p$score)/2
TSS$class <-rep_len(c("TSS"), dim(TSS)[1]) #bind a column just saying the class to the df.

Dist_m <- m_long %>% filter(X2=="qnd_Dist_skew_minus.bed") %>% select(!X2)
Dist_p <- p_long %>% filter(X2=="qnd_Dist_skew_plus.bed") %>% select(!X2)
Dist <- Dist_m
Dist$score <- (Dist_m$score+Dist_p$score)/2
Dist$class <-rep_len(c("Dist"), dim(Dist)[1])
classed_skew <- rbind(TSS, Dist) #bind the dfs back together

ggplot(classed_skew, aes(x=Bins, y=score, color=class)) + #plot bins by score separated by class.
  theme_bw() + geom_line(size=1.5) +
  scale_x_continuous(breaks=seq(
    from=0, to=150, by=25)) +
  scale_y_continuous(breaks=seq(
    from=-1, to=1, by=0.05)) +
  annotate("text", x=25, y=-0.15, label="TSS", size=5) +
  annotate("text", x=125, y=-0.15, label="3' CpG", size=5)+
  annotate("text", x=150, y=-0.15, label="+ 0.25 kb", size=5)+
  annotate("text", x=0, y=-0.15, label="- 0.25 kb", size=5) + 
  ylab("Skew")

TSS_m <- Groseq_m %>% filter(X2=="qnd_TSS_skew_minus.bed") %>% select(!X2)
TSS_p <- Groseq_p %>% filter(X2=="qnd_TSS_skew_plus.bed") %>% select(!X2) #split into TSS / Dist
TSS <- TSS_m #bind plus and minus of a given class, here TSS, together.
TSS$score <- (TSS_m$score+TSS_p$score)/2
TSS$class <-rep_len(c("TSS"), dim(TSS)[1]) #bind a column just saying the class to the df.

Dist_m <- Groseq_m %>% filter(X2=="qnd_Dist_skew_minus.bed") %>% select(!X2)
Dist_p <- Groseq_p %>% filter(X2=="qnd_Dist_skew_plus.bed") %>% select(!X2)
Dist <- Dist_m
Dist$score <- (Dist_m$score+Dist_p$score)/2
Dist$class <-rep_len(c("Dist"), dim(Dist)[1])
classed_groseq <- rbind(TSS, Dist) #bind the dfs back together

ggplot(classed_groseq, aes(x=Bins, y=score, color=class)) + #plot bins by score separated by class.
  theme_bw() + geom_line(size=1.5) +
  scale_x_continuous(breaks=seq(
    from=0, to=150, by=25)) +
  annotate("text", x=25, y=-0.15, label="TSS", size=5) +
  annotate("text", x=125, y=-0.15, label="3' CpG", size=5)+
  annotate("text", x=150, y=-0.15, label="+ 0.25 kb", size=5)+
  annotate("text", x=0, y=-0.15, label="- 0.25 kb", size=5) + 
  ylab("Groseq")

########################################################################################################################
TSS <- H3K4me3 %>% filter(X2=="qnd_Dist_skew.bed") %>% select (!X2)
TSS$class <-rep_len(c("TSS"), dim(TSS)[1])
Dist <- H3K4me3 %>% filter(X2=="qnd_TSS_skew.bed") %>% select (!X2)
Dist$class <-rep_len(c("Dist"), dim(Dist)[1])
H3K4me3 <- rbind(TSS, Dist) #bind the dfs back together

ggplot(H3K4me3, aes(x=Bins, y=score, color=class)) + #plot bins by score separated by class.
  theme_bw() + geom_line(size=1.5) +
  scale_x_continuous(breaks=seq(
    from=0, to=150, by=25)) +
  annotate("text", x=25, y=-0.15, label="TSS", size=5) +
  annotate("text", x=125, y=-0.15, label="3' CpG", size=5)+
  annotate("text", x=150, y=-0.15, label="+ 0.25 kb", size=5)+
  annotate("text", x=0, y=-0.15, label="- 0.25 kb", size=5) + 
  ylab("H3K4me3")
#I would like to figure out how to put error bars into this... I know it's taking the mean but is it taking counts or what?
#right now deeptools is too much of a black box to me, I need to fix that. I know it's mean score but can I get more?
TSS <- H3K4me2 %>% filter(X2=="qnd_Dist_skew.bed") %>% select (!X2)
TSS$class <-rep_len(c("TSS"), dim(TSS)[1])
Dist <- H3K4me2 %>% filter(X2=="qnd_TSS_skew.bed") %>% select (!X2)
Dist$class <-rep_len(c("Dist"), dim(Dist)[1])
H3K4me2 <- rbind(TSS, Dist) #bind the dfs back together

ggplot(H3K4me2, aes(x=Bins, y=score, color=class)) + #plot bins by score separated by class.
  theme_bw() + geom_line(size=1.5) +
  scale_x_continuous(breaks=seq(
    from=0, to=150, by=25)) +
  annotate("text", x=25, y=-0.15, label="TSS", size=5) +
  annotate("text", x=125, y=-0.15, label="3' CpG", size=5)+
  annotate("text", x=150, y=-0.15, label="+ 0.25 kb", size=5)+
  annotate("text", x=0, y=-0.15, label="- 0.25 kb", size=5) + 
  ylab("H3K4me2")

TSS <- H3K79me2 %>% filter(X2=="qnd_Dist_skew.bed") %>% select (!X2)
TSS$class <-rep_len(c("TSS"), dim(TSS)[1])
Dist <- H3K79me2 %>% filter(X2=="qnd_TSS_skew.bed") %>% select (!X2)
Dist$class <-rep_len(c("Dist"), dim(Dist)[1])
H3K79me2 <- rbind(TSS, Dist) #bind the dfs back together

ggplot(H3K79me2, aes(x=Bins, y=score, color=class)) + #plot bins by score separated by class.
  theme_bw() + geom_line(size=1.5) +
  scale_x_continuous(breaks=seq(
    from=0, to=150, by=25)) +
  annotate("text", x=25, y=-0.15, label="TSS", size=5) +
  annotate("text", x=125, y=-0.15, label="3' CpG", size=5)+
  annotate("text", x=150, y=-0.15, label="+ 0.25 kb", size=5)+
  annotate("text", x=0, y=-0.15, label="- 0.25 kb", size=5) + 
  ylab("H3K79me2")


##########################
## Stranded matrix data ##
##########################
#information on which group is which is held in the introductory string, which I need to figure out how to access.
#I can view it by eye obviously but that's inefficient so we kinda want to avoid that.
# long_tab <- function (tab, start) {
#   reformat_tab <- tab %>%
#     pivot_longer(cols=start:length(tab),
#                  names_to = "Bins",
#                  names_prefix="X",
#                  values_to = "score",
#                  values_drop_na = TRUE) %>% #pivots table into longform with bins as one column and values as another
#     mutate(Bins=(as.numeric(Bins)-2), #because extra columns interfere with the bin naming (bin 1 is imported as X3 for example)
#            .keep="unused") #need to fix the numbering by taking X off and subtracting 2
#   return(reformat_tab)
# }
# 
# format_deeptools_matrix_tabs <- function (tab_file) {
#   tab <- read_table(tab_file, col_names = FALSE, skip = 1)
#   #get metadata string from first line of file
#   first_line <- readLines(tab_file, n=1)
#   extracted <- unlist(str_extract_all(first_line, "\\w*")) #this makes a vector with punctuation removed
#   joined <- str_squish(str_flatten(extracted, collapse = " ")) #now we have a cleaner, easier to parse string
#   #get group labels
#   labels <- str_trim(str_extract(joined, "(?<=group_labels).*(?=group_boundaries)"))
#   labels_vec <- as.list(str_subset(str_split_fixed(labels, " ", n=Inf), "[^(bed)]"))
#   #get group bounds
#   group_bounds <- str_trim(str_extract(joined, "(?<=group_boundaries).*(?=sample_labels)"))
#   nums <- as.numeric(str_split_fixed(group_bounds, " ", n=Inf))
#   split_by_group <- split(tab, cumsum(1:nrow(tab) %in% (nums+1)))
#   long_tab_list <- lapply(split_by_group, long_tab, start=7)
#   for (i in 1:length(long_tab_list)) {
#     d1 = dim(long_tab_list[[i]])
#     class = (rep(labels_vec[[i]], d1[1]))
#     long_tab_list[[i]] = cbind.data.frame(long_tab_list[[i]], class) #yay this works to get a class column
#   }
#   return(long_tab_list)
# }
# 
# skew_plus <- format_deeptools_matrix_tabs("skew_10_plus.matrix")
# skew_minus <- format_deeptools_matrix_tabs("skew_10_minus.matrix")
# 
# TSS_m <- skew_minus[[1]]
# TSS_p <- skew_plus[[1]]
# TSS <- rbind(TSS_m, TSS_p) #bind plus and minus of a given class, here TSS, together.
# TSS$classes <-rep_len(c("TSS"), dim(TSS)[1]) #bind a column just saying the class to the df. B/C as much as I tried
# #the current one is formatted horribly. It does function as a dblchck to make sure I didn't screw anything up, at least.
# 
# Dist_m <- skew_minus[[2]]
# Dist_p <- skew_plus[[2]]
# Dist <- rbind(Dist_m, Dist_p)
# Dist$classes <-rep_len(c("Dist"), dim(Dist)[1])
# Dist$end_bin <- ((Dist$X3-Dist$X2)/10)+100
# Dist <- arrange(Dist, (Dist$X3-Dist$X2))
# 
# #okay testing method, do this instead of full matrix.
# #wide_for_sampling <- pivot_wider(classed_matrix, names_from = Bins, values_from = score, names_prefix = "H")
# #samp <- sample(nrow(wide_for_sampling), size=10, replace = FALSE) 
# #subs <- wide_for_sampling[samp,]
# #long_subs <- subs %>% pivot_longer(cols=9:408, names_to = "Bins", values_to = "score", names_prefix = "H")
# 
# ggplot(Dist, aes(y=X4)) + #plot bins by score separated by class.
#   geom_raster(interpolate=TRUE, aes(fill=score, x=as.numeric(Bins)))
# #this just doesn't work, it won't sort and the size is too big.