library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)

fasta_file <- "C:\\Users\\16421\\Desktop\\123\\CDS.fasta"
output_dir <- "C:\\Users\\16421\\Desktop\\123\\"

cds_seqs <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)
gene_names <- names(cds_seqs)
sequences <- sapply(cds_seqs, function(x) x[1])

processed_seqs <- c()
for (i in 1:length(sequences)) {
  seq <- toupper(sequences[i])
  seq_len <- nchar(seq)
  
  if (seq_len >= 3) {
    last3 <- substr(seq, seq_len-2, seq_len)
    if (last3 %in% c("TAA", "TAG")) {
      seq <- substr(seq, 1, seq_len-3)
    }
  }
  
  new_len <- nchar(seq)
  rem <- new_len %% 3
  if (rem != 0) {
    seq <- substr(seq, 1, new_len - rem)
  }
  
  processed_seqs[i] <- seq
}
names(processed_seqs) <- gene_names

codon_df <- data.frame()
for (gene in gene_names) {
  seq <- processed_seqs[[gene]]
  s <- strsplit(seq, "")[[1]]
  n <- length(s)
  
  pos1 <- s[seq(1, n, 3)]
  pos2 <- s[seq(2, n, 3)]
  pos3 <- s[seq(3, n, 3)]
  
  calc <- function(p) {
    n <- length(p)
    c(A = sum(p == "A")/n*100,
      T = sum(p == "T")/n*100,
      G = sum(p == "G")/n*100,
      C = sum(p == "C")/n*100)
  }
  
  p1 <- calc(pos1)
  p2 <- calc(pos2)
  p3 <- calc(pos3)
  
  codon_df <- rbind(codon_df, data.frame(
    Gene = gene,
    Pos1A = p1["A"], Pos1T = p1["T"], Pos1G = p1["G"], Pos1C = p1["C"],
    Pos2A = p2["A"], Pos2T = p2["T"], Pos2G = p2["G"], Pos2C = p2["C"],
    Pos3A = p3["A"], Pos3T = p3["T"], Pos3G = p3["G"], Pos3C = p3["C"]
  ))
}

plot_data <- codon_df %>%
  pivot_longer(cols = -Gene, names_to = "Feature", values_to = "Content") %>%
  mutate(Feature = factor(Feature, levels = c(
    "Pos1A","Pos1T","Pos1G","Pos1C",
    "Pos2A","Pos2T","Pos2G","Pos2C",
    "Pos3A","Pos3T","Pos3G","Pos3C"
  )))

p <- ggplot(plot_data, aes(x = Gene, y = Feature)) +
  geom_tile(aes(fill = Content), color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Content, 1)), size = 3, color = "black") +
  scale_fill_gradient(low = "#87CEEB", high = "#FA8072") + 
  labs(x = "Mitochondrial PCGs", y = "Codon Position", fill = "Content (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 10, face = "italic"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    
    panel.grid = element_blank()
  )

ggsave(
  filename = paste0(output_dir, "Codon_Usage_Heatmap.png"),
  plot = p,
  width = 10,
  height = 6,
  dpi = 900
)
