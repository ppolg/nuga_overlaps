###########################################################################

# Checking the overlaps of Mtb genes:
#   • Read in FASTA, GFF
#   • Search ORFs for which overlap
#   • Filter out 4-nt
#   • Check for RBSs as 5 consecutive purine 5-15 nt upstream of TSS
#   • Make necessary tables/figures
#     - Export them

###########################################################################

########
# Libs #
########

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","Biostrings")
invisible(lapply(packages, require, character.only = TRUE))

########
# Func #
########

# Read gff and optionally filter for type
get_gff_input <- function(gff_file,filter_type="none") {
  input_gff <- paste(here::here("NUGA_data/"), gff_file, ".gff", sep = "") 
  gff <- read.gff(input_gff)
  
  # Filter if argument passed, to get only rRNA, tRNA, CDS etc...
  if(filter_type != "none"){ 
    gff %>% dplyr::filter(type == filter_type)
  }
  else{invisible(gff)}
}

# nin - because it just looks clean and neat!
`%nin%` = Negate(`%in%`)

########
# Init #
########


# For classifying nucleotides as pyrimidine/purine
replacements <- c("A"="R","G"="R","T"="Y","U"="Y","C"="Y" )

# For checking RBS:
RBS_seq <- "RRRRR"   #  The strand of purines
RBS_range <- c(5,19) #  The range in which we look for purines

# Range of plot colours
colours <- hash::hash(
  red = "#A3280A",
  orange = "#E3812B",
  brown = "#8A5122",
  yellow = "#E0D253",
  grey = "#858482",
  green = "#195928",
  blue = "#4A749E",
  purple = "#612882"
)

# ggplot2 theme for figures
theme_alex <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_blank(),
           axis.title = element_text(face = "bold",size = rel(1.4)),
           axis.title.y = element_text(angle=90,vjust = 3, size = rel(1)),
           axis.title.x = element_text(vjust = -0.1, size = rel(1)),
           axis.text = element_text(face="bold",size = rel(1.3)), 
           axis.line = element_line(colour="black",size=1),
           axis.ticks = element_line(colour="black",size = 1),
           axis.ticks.length = unit(.35, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "top",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.x = unit(0.3, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="italic", color = "black"),
           plot.margin=unit(c(12,6,6,6),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}

########
# Main #
########

# Read Mtb H37Rv FASTA (from NCBI)
mtb_fasta <- readDNAStringSet(here("NUGA_data/Mtb.fasta"))

mtb_fasta_string_plus <- toString(mtb_fasta)
mtb_fasta_string_minus <- toString(reverseComplement(mtb_fasta))

# Read Mtb H37Rv GFF (From Mycobrowser release 4)
mtb_gff <- get_gff_input("Mtb") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";Function=", 2)[,1], "Name=", 2)[,2],) %>%
  mutate(Rv_name = str_split_fixed(str_split_fixed(attributes, ";Name=", 2)[,1], "Locus=", 2)[,2],) %>%
  mutate(product_name = str_split_fixed(str_split_fixed(attributes, ";Comments=", 2)[,1], "Product=", 2)[,2],) %>%
  mutate(function_group = str_split_fixed(str_split_fixed(attributes, ";Protein Data Bank", 2)[,1], "Functional_Category=", 2)[,2],) %>%
  dplyr::filter(start != 1960667) %>% #Two ncRNAs with same TSS messing up my code
  dplyr::filter(locus_name != "oxyR'") %>% #Misannotated pseudogene
  filter(type == "CDS") %>%
  select(type,start,end,strand,attributes,locus_name, Rv_name, product_name, function_group)

#-------------#
# PLUS STRAND #
#-------------#

gff_plus <- mtb_gff %>%
  dplyr::filter(strand == "+") %>%
  dplyr::mutate(
    prev_start = dplyr::lag(start, order_by=start),
    prev_end = dplyr::lag(end, order_by=start),
    next_start = dplyr::lead(start, order_by=start),
    next_end = dplyr::lead(end, order_by=start),
    overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
    overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(start, .by_group = TRUE)

# Grab upstream for purines

gff_plus <- gff_plus %>%
  dplyr::mutate(seq_upstream = toString(subseq(mtb_fasta,start-RBS_range[2],start-RBS_range[1])))

gff_plus$seq_upstream <- str_replace_all(gff_plus$seq_upstream,replacements)

gff_plus <- gff_plus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

# Get overlap

plus_overlaps <- gff_plus %>%
  dplyr::filter(overlap_next == TRUE | overlap_prev == TRUE ) %>%
  dplyr::mutate(
    n_overlap_prev = ifelse(overlap_prev == TRUE,prev_end - start +1,NA),
    n_overlap_next = ifelse(overlap_next == TRUE,end - next_start +1,NA),
    seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(subseq(mtb_fasta,start,prev_end)), NA),
    seq_overlap_next = ifelse(overlap_next == TRUE, toString(subseq(mtb_fasta,next_start,end)), NA)
  )

# Group into overlapping operons
#   i.e. if multiple overlaps, still same group

plus_overlaps$group = NA
plus_overlaps$group[1] = as.numeric(1)
for (row in 2:nrow(plus_overlaps)) {
  if (plus_overlaps$overlap_prev[row] == TRUE) {
    plus_overlaps$group[row] <- plus_overlaps$group[row-1]
  }
  else {
    plus_overlaps$group[row] <- plus_overlaps$group[row-1] + 1
  }
}

#--------------#
# MINUS STRAND #
#--------------#

gff_minus <- mtb_gff %>%
  dplyr::filter(strand == "-") %>%
  dplyr::mutate(
    prev_start = dplyr::lead(start, order_by=start),
    prev_end = dplyr::lead(end, order_by=start),
    next_start = dplyr::lag(start, order_by=start),
    next_end = dplyr::lag(end, order_by=start),
    overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
    overlap_next = ifelse(start <= next_end & start > next_start,TRUE,FALSE )) %>%
  group_by(start) %>%
  dplyr::arrange(desc(start))

# Grab range upstream

gff_minus <- gff_minus %>%
  dplyr::mutate(seq_upstream = toString(reverseComplement(subseq(mtb_fasta,end+RBS_range[1],end+RBS_range[2]))))

gff_minus$seq_upstream <- str_replace_all(gff_minus$seq_upstream,replacements)

gff_minus <- gff_minus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

# Get overlap

minus_overlaps <- gff_minus %>%
  dplyr::filter(overlap_next == TRUE | overlap_prev == TRUE ) %>%
  dplyr::mutate(
    n_overlap_prev = ifelse(overlap_prev == TRUE,end - prev_start +1,NA),
    n_overlap_next = ifelse(overlap_next == TRUE,next_end - start +1,NA),
    seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(reverseComplement(subseq(mtb_fasta,prev_start,end))), NA),
    seq_overlap_next = ifelse(overlap_next == TRUE, toString(reverseComplement(subseq(mtb_fasta,start,next_end))), NA)
  )

# Group into overlapping operons
#   i.e. if multiple overlaps, still same group

minus_overlaps$group = NA
minus_overlaps$group[1] = as.numeric(plus_overlaps$group[nrow(plus_overlaps)]+1)
for (row in 2:nrow(minus_overlaps)) {
  if (minus_overlaps$overlap_prev[row] == TRUE) {
    minus_overlaps$group[row] <- minus_overlaps$group[row-1]
  }
  else {
    minus_overlaps$group[row] <- minus_overlaps$group[row-1] + 1
  }
}

#----------#
# COMBINED #
#----------#

#Get all gff together
gff_annot <- rbind(gff_plus,gff_minus) %>%
  mutate(start_codon = ifelse(strand == "+",toString(subseq(mtb_fasta,start,start+2)),toString(reverseComplement(subseq(mtb_fasta,end-2,end)))))

# Get all overlaps together
overlaps_all <- rbind(plus_overlaps,minus_overlaps)

# overlaps of 4 nt
overlaps_4 <- overlaps_all %>%
  dplyr::filter(n_overlap_prev==4)

# For the pie chart of 4nt overlaps
overlaps_4_count <- overlaps_4 %>%
  dplyr::select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(seq_overlap_prev) %>%
  group_by(seq_overlap_prev) %>%
  summarise(n = length(seq_overlap_prev))

# For the pie chart of 4nt overlaps, no RBS
overlaps_4_norbs_count <- overlaps_4 %>%
  dplyr::filter(has_RBS == FALSE) %>%
  select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(seq_overlap_prev) %>%
  group_by(seq_overlap_prev) %>%
  summarise(n = length(seq_overlap_prev))

# For the pie chart of all start codons

gff_weirds <- gff_annot %>% # Select all that are non-NTG start
  filter(start_codon %nin% c("ATG","GTG","TTG","CTG"))

gff_noweirds <- gff_annot %>% # Select everyhting that IS NTG start
  filter(start %nin% gff_weirds$start)

start_codon_count <- gff_noweirds %>%
  dplyr::select(type,start,end,strand,locus_name,start_codon) %>%
  ungroup() %>%
  arrange(start_codon) %>%
  group_by(start_codon) %>%
  summarise(n = length(start_codon)) %>%
  add_row(start_codon = "other", n=nrow(gff_weirds)) %>%
  mutate(start_codon = factor(start_codon, levels=c("ATG","GTG","TTG","other")))


# Categorise overlap function - all overlaps
overlaps_func <- overlaps_all %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset= "overlap")

# Categorise overlap function - 4 nt
overlaps_func_4 <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_4")

# Categorise overlap function - 4 nt, no RBS
overlaps_func_4_noRBS <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "" & has_RBS == FALSE) %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_4_noRBS")

# Categorise overlap function - every ORF, ever
total_func <- mtb_gff %>%
  dplyr::filter( function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset="all")

# Combine
func_all <- rbind(overlaps_func,overlaps_func_4, overlaps_func_4_noRBS, total_func)

########
# Plot #
########

# Overlap length distribution
ggplot(overlaps_all,aes(n_overlap_next)) +
  geom_histogram(binwidth = 1, color="black", size=1, fill="grey80") +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Count (log10)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_log10(expand = c(0,0), breaks = c(1,10,100,1000), limits=c(1,1e3)) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# Function categories distribution
ggplot(func_all, aes(x=function_group, y=n_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_alex() +
  coord_flip() +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits=c(0,32),expand=c(0,0),breaks=c(0,10,20,30), labels = c("0%","10%","20%","30%")) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))

############
# Stat/Out #
############

# Hypergeo for functional groups
func_annot <- as.data.frame(func_all) %>% filter(dataset == "all") %>% select(function_group,n)
func_overlap <- as.data.frame(func_all) %>% filter(dataset == "overlap") %>% select(function_group,n)
func_4nt <- func_all %>% filter(dataset == "overlap_4") %>% select(function_group,n)
func_norbs <- func_all %>% filter(dataset == "overlap_4_noRBS") %>% select(function_group,n)

func_4nt <- func_4nt %>% add_row (function_group = "unknown", n=0, .before = 9)
func_norbs <- func_norbs %>% add_row (function_group = "unknown", n=0, .before = 9)

func_pvals <- func_annot %>% select(function_group)

p1 = c()
for(row in 1:nrow(func_overlap)){
  n = 1.0-phyper(func_overlap$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_overlap$n))
  p1 <- base::append(p1,n)
}

p2 = c()
for(row in 1:nrow(func_4nt)){
  n = 1.0-phyper(func_4nt$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_4nt$n))
  p2 <- base::append(p2,n)
}

p3 = c()
for(row in 1:nrow(func_norbs)){
  n = 1.0-phyper(func_norbs$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_norbs$n))
  p3 <- base::append(p3,n)
}

p1_fdr <- p.adjust(p1,method = "BH")
p2_fdr <- p.adjust(p2,method = "BH")
p3_fdr <- p.adjust(p3,method = "BH")

p1_by <- p.adjust(p1,method = "BY")
p2_by <- p.adjust(p2,method = "BY")
p3_by <- p.adjust(p3,method = "BY")

func_pvals <- func_annot %>% select(function_group) %>%
  mutate(p_overlap := p1,
         p_overlap_4nt := p2,
         p_overlap_4bt_norbs := p3,
         p_overlap_BH := p1_fdr,
         p_overlap_4nt_BH := p2_fdr,
         p_overlap_4bt_norbs_BH := p3_fdr,
         p_overlap_BY := p1_by,
         p_overlap_4nt_BY := p2_by,
         p_overlap_4bt_norbs_BY := p3_by)

#------------#
# Print figs #
#------------#



#-------------#
# Print lists #
#-------------#

# Write lists of overlaps - all
print_overlaps <- overlaps_all %>%
  select(locus_name,product_name,start,end,strand,overlap_prev,n_overlap_prev,seq_overlap_prev,
         overlap_next,n_overlap_next,seq_overlap_next,group) %>%
  dplyr::rename("Locus name" = locus_name,
                "Product" = product_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand,
                "Overlaps ORF downstream" = overlap_next,
                "Downstream overlap length" = n_overlap_next,
                "Sequence overlapping downstream" = seq_overlap_next,
                "Overlaps ORF upstream" = overlap_prev,
                "Upstream overlap length" = n_overlap_prev,
                "Sequence overlapping upstream" = seq_overlap_prev,
                "Group" = group)

fwrite(print_overlaps, file = here::here("NUGA_Out/overlaps_all.csv"))

# Write lists of overlaps - 4nt
print_overlaps_4 <- overlaps_all %>%
  filter(n_overlap_prev == 4 | n_overlap_next == 4) %>%
  select(locus_name,product_name,start,end,strand,overlap_prev,seq_overlap_prev,
         overlap_next,seq_overlap_next,group) %>%
  dplyr::rename("Locus name" = locus_name,
                "Product" = product_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand,
                "Overlaps ORF downstream" = overlap_next,
                "Sequence overlapping downstream" = seq_overlap_next,
                "Overlaps ORF upstream" = overlap_prev,
                "Sequence overlapping upstream" = seq_overlap_prev,
                "Group" = group)

fwrite(print_overlaps_4, file = here::here("NUGA_Out/overlaps_4.csv"))

# Write tables for pie charts to make later

fwrite(start_codon_count, file = here::here("NUGA_Out/pie_start_codons.csv"))
fwrite(overlaps_4_count, file = here::here("NUGA_Out/pie_4nt.csv"))
fwrite(overlaps_4_norbs_count, file = here::here("NUGA_Out/pie_noRBS.csv"))

# Write table for p values
fwrite(func_pvals, file = here::here("NUGA_Out/p_values.csv"))
