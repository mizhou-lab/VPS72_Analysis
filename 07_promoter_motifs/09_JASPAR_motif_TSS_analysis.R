# Analyzes JASPAR motif hits within ±50 bp of the VPS72 TSS and tests TF class enrichment relative to the full promoter background.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(writexl)
})

# Paths
base_dir <- "data/jaspar"
in_csv   <- file.path(base_dir, "jaspar_data.csv")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# Load and prep
# NOTE: if start/end are 1..1000 for a -1000..-1 window, use: mid_pos = (start + end)/2 - 1001
jas <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(mid_pos = (start + end) / 2)

rel_thr <- 0.85

bg     <- jas %>% filter(relScore >= rel_thr)
tss50  <- bg  %>% filter(mid_pos >= -50, mid_pos <= 50)

# Export TSS ±50bp hits
write_xlsx(list("All_hits_±50bp_TSS" = tss50 %>% arrange(desc(relScore))),
           file.path(base_dir, "SuppTable_VPS72_JASPAR_TSS_50bp.xlsx"))

# Correlation scatter
p_scatter <- ggplot(tss50, aes(x = mid_pos, y = relScore)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Distance from VPS72 TSS (bp)", y = "JASPAR motif relScore",
       title = "Predicted TF motif strength near the VPS72 TSS (±50 bp)") +
  theme_classic(base_size = 12)

# Label top TFs
label_tf <- tss50 %>%
  group_by(TF) %>%
  summarise(n_sites = n(), max_relScore = max(relScore),
            mid_pos = mid_pos[which.max(relScore)][1], .groups = "drop") %>%
  filter(n_sites >= 3, max_relScore >= 0.9) %>%
  slice_head(n = 10)

p_labeled <- p_scatter +
  geom_text_repel(data = label_tf, aes(x = mid_pos, y = max_relScore, label = TF),
                  size = 3, max.overlaps = 10, box.padding = 0.4,
                  point.padding = 0.3, segment.color = "grey50")

ggsave(file.path(base_dir, "SuppFig_VPS72_JASPAR_TSS_50bp.png"),         p_scatter,  width = 6.5, height = 4.5, dpi = 300)
ggsave(file.path(base_dir, "SuppFig_VPS72_JASPAR_TSS_50bp_labeled.png"), p_labeled,  width = 6.5, height = 4.5, dpi = 300)

# TF class enrichment (Fisher exact: ±50 bp vs promoter background)
classes <- sort(unique(na.omit(bg$class)))

enrich <- lapply(classes, function(cl) {
  a <- sum(tss50$class == cl, na.rm = TRUE)
  b <- sum(tss50$class != cl, na.rm = TRUE)
  c <- sum(bg$class == cl,    na.rm = TRUE) - a
  d <- sum(bg$class != cl,    na.rm = TRUE) - b
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  data.frame(class = cl, n_tss50 = a, n_bg = sum(bg$class == cl, na.rm = TRUE),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}) %>%
  bind_rows() %>%
  mutate(padj = p.adjust(p, method = "BH"),
         log2OR = log2(odds_ratio),
         class = ifelse(is.na(class) | class == "", "Unclassified", class)) %>%
  arrange(padj)

write_csv(enrich, file.path(base_dir, "SuppTable_VPS72_TFclass_enrichment_TSS50.csv"))

p_enrich <- ggplot(enrich %>% slice_head(n = 12), aes(x = reorder(class, log2OR), y = log2OR)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(size = n_tss50), alpha = 0.9) +
  coord_flip() +
  labs(x = NULL, y = "log2 Odds ratio (±50 bp vs promoter background)",
       title = "TF class enrichment near VPS72 TSS",
       subtitle = paste0("relScore >= ", rel_thr, "; point size = # motifs within ±50 bp")) +
  theme_classic(base_size = 12)

ggsave(file.path(base_dir, "SuppFig_VPS72_TFclass_enrichment_TSS50.png"), p_enrich, width = 10, height = 4.6, dpi = 300)
ggsave(file.path(base_dir, "SuppFig_VPS72_TFclass_enrichment_TSS50.pdf"), p_enrich, width = 10, height = 4.6)

message("JASPAR outputs saved to: ", base_dir)
