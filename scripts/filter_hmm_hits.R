#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# Set default values for optional parameters
default_window_size <- 100000
default_ratio_cutoff <- 50
default_min_score <- 200
default_min_size <- 200

# Check for required arguments and optional parameters
if (length(args) < 5) {
  stop(paste("Usage: Rscript UCE_hmm_hits_filter.R <hits.txt> <assembly.fasta> <species> <out_dir> <figure_dir>",
             "[window_size=", default_window_size, "]",
             "[ratio_cutoff=", default_ratio_cutoff, "]",
             "[min_score=", default_min_score, "]",
             "[min_size=", default_min_size, "]"), call.=FALSE)
}

## Set input and output files
## Inputs
hits_file <- args[1]
assembly_file <- args[2]
species <- args[3]

## Outputs
out_dir <- args[4]
figure_dir <- args[5]

## Optional parameters with defaults
win_size <- if(length(args) >= 6) as.numeric(args[6]) else default_window_size
ratio_cutoff <- if(length(args) >= 7) as.numeric(args[7]) else default_ratio_cutoff
min_score <- if(length(args) >= 8) as.numeric(args[8]) else default_min_score
min_size <- if(length(args) >= 9) as.numeric(args[9]) else default_min_size

library(Biostrings)
library(GenomicRanges)
library(tidyverse)

## Read hmm hits
hits_txt <- read_table(hits_file,
                       comment = "#",
                       col_names = ("chr v1 locus v2 UCE_from UCE_to hit_from hit_to envfrom envto seqlength strand eval score bias des" %>%
                                      str_split(" "))[[1]])

## Remove ambiguous hits
multiple_hits_loci <- hits_txt %>%
  group_by(locus) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>%
  pull(locus)

p1 <- hits_txt %>%
  filter(locus %in% multiple_hits_loci) %>%
  group_by(locus) %>%
  arrange(eval) %>%
  mutate(idx = 1:n()) %>%
  filter(idx <=2) %>%
  mutate(idx = ifelse(idx == 1, "first", "second")) %>%
  select(locus, idx, eval) %>%
  pivot_wider(names_from = idx, values_from = eval) %>%
  mutate(ratio = first/second) %>%
  mutate(ratio = -log10(ratio)) %>%
  ggplot(aes(x = ratio)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = quantile(ratio, 0.05)),linetype = "dashed")

R.devices::suppressGraphics(ggsave(file.path(figure_dir, paste0("first_to_second_eval_ratio_", species, ".png")), p1))

loci_remove_abiguity <- hits_txt %>%
  filter(locus %in% multiple_hits_loci) %>%
  group_by(locus) %>%
  arrange(eval) %>%
  mutate(idx = 1:n()) %>%
  filter(idx <=2) %>%
  mutate(idx = ifelse(idx == 1, "first", "second")) %>%
  select(locus, idx, eval) %>%
  pivot_wider(names_from = idx, values_from = eval) %>%
  mutate(ratio = first/second) %>%
  mutate(ratio = -log10(ratio)) %>%
  filter(ratio <= ratio_cutoff) %>%
  pull(locus)

## Remove low score hits
p2 <- hits_txt %>%
  filter(!locus %in% loci_remove_abiguity) %>%
  group_by(locus) %>%
  filter(eval == min(eval)) %>%
  ggplot(aes(x = score)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = quantile(score, 0.05)), linetype = "dashed")

R.devices::suppressGraphics(ggsave(file.path(figure_dir, paste0("hits_score_", species, ".png")), p2))

loci_remove_lowScore <- hits_txt %>%
  filter(!locus %in% loci_remove_abiguity) %>%
  group_by(locus) %>%
  filter(eval == min(eval)) %>%
  filter(score < min_score) %>%
  pull(locus)

## Remove short hits
p3 <- hits_txt %>%
  filter(!locus %in% loci_remove_abiguity) %>%
  filter(!locus %in% loci_remove_lowScore) %>%
  group_by(locus) %>%
  filter(eval == min(eval)) %>%
  mutate(size = abs(hit_to - hit_from)) %>%
  ggplot(aes(x = size)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = quantile(size, 0.05)), linetype = "dashed")
R.devices::suppressGraphics(ggsave(file.path(figure_dir, paste0("hits_size_", species, ".png")), p3))

loci_remove_tooSmall <- hits_txt %>%
  filter(!locus %in% loci_remove_abiguity) %>%
  filter(!locus %in% loci_remove_lowScore) %>%
  group_by(locus) %>%
  filter(eval == min(eval)) %>%
  mutate(size = abs(hit_to - hit_from)) %>%
  filter(size < min_size) %>%
  pull(locus)

## Get good hits

accurate_hits <- hits_txt %>%
  filter(!locus %in% loci_remove_abiguity) %>%
  filter(!locus %in% loci_remove_lowScore) %>%
  filter(!locus %in% loci_remove_tooSmall) %>%
  group_by(locus) %>%
  filter(eval == min(eval))

## Plot distribution on genome
p4 <- accurate_hits %>%
  select(locus, chr, hit_from, hit_to) %>%
  mutate(pos = (hit_to+hit_from)/2) %>%
  mutate(bin = round(pos/win_size, digits = 0)) %>%
  group_by(chr, bin) %>%
  summarize(n = n()) %>%
  filter(n != 0) %>%
  mutate(pos = bin*win_size) %>%
  ggplot(aes(x = pos/10^6, y = n)) +
  geom_point() +
  facet_wrap(~chr, ncol = 1, scales = "free") +
  xlab("Position(Mbp)") +
  ylab(paste0("Number of UCE loci in ", win_size/1000, "Kbp windows")) +
  ylim(0, 10) +
  ggtitle(species)
R.devices::suppressGraphics(ggsave(file.path(figure_dir, paste0("uce_distribution_", win_size/1000, "Kbp_", species, ".png")), p4))

## Write hits to file
accurate_hits %>%
  select(locus, chr, hit_from, hit_to) %>%
  write_tsv(file.path(out_dir, paste0("uce_positions_", species, ".tsv")))

## Extract sequences
assembly_bss <- readDNAStringSet(assembly_file)
names(assembly_bss) <- names(assembly_bss) %>% str_remove_all(" .+") %>% str_remove_all(" $")

## Write each uce to a separate file
fasta_dir <- file.path(out_dir, species)
dir.create(fasta_dir, showWarnings = FALSE)
for (i in 1:nrow(accurate_hits)){
  seq_hitted <- subseq(assembly_bss[[accurate_hits[[i,"chr"]]]],
                       start = min(accurate_hits[[i,"hit_from"]], accurate_hits[[i,"hit_to"]]),
                       end = max(accurate_hits[[i,"hit_from"]], accurate_hits[[i,"hit_to"]]))
  if(accurate_hits[[i,"strand"]] == "-"){
    seq_hitted <- reverseComplement(seq_hitted)
  }
  locus <- accurate_hits[[i,"locus"]]
  locus_bss <- DNAStringSet(seq_hitted)
  names(locus_bss) <- species
  locus_bss %>%
    writeXStringSet(file.path(fasta_dir, paste0(locus, ".fasta")))
}
