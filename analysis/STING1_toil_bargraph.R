###############################################################################
# STING1 (TMEM173) expression: tumor vs normal, bar graphs from REAL data
#
# Data: UCSC Xena "TCGA TARGET GTEx" (Toil recompute) -- TCGA and GTEx samples
#       processed through the same RSEM pipeline, so tumor and normal are
#       directly comparable.
#   expression: TcgaTargetGtex_rsem_gene_tpm   (values are log2(TPM + 0.001))
#   phenotype:  TcgaTargetGTEX_phenotype.txt
#   GENCODE v23 is used by Toil, so STING1 is still called TMEM173 there.
#
# Output:
#   Fig_STING1_pancancer_bar.pdf/.png  -- 9 tumor types, Tumor vs Normal
#   Fig_STING1_OV_FT_bar.pdf/.png      -- TCGA-OV vs GTEx fallopian tube / ovary
#   STING1_toil_tidy.csv               -- one row per sample (for Prism, etc.)
#   STING1_toil_stats.csv              -- n, mean, SD, median, log2FC, p, FDR
#
# Bars = mean +/- SD, every sample shown as a dot.
# Tests: two-sided Wilcoxon rank-sum, BH-adjusted across tumor types.
#
# install.packages(c("ggplot2", "dplyr", "tidyr", "UCSCXenaTools"))
###############################################################################
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(UCSCXenaTools)
})

outdir <- "."
cache  <- "xena_cache"

# ---- 1. Download ------------------------------------------------------------
host     <- "https://toil.xenahubs.net"
expr_ds  <- "TcgaTargetGtex_rsem_gene_tpm"
pheno_ds <- "TcgaTargetGTEX_phenotype.txt"

# Only the one gene is fetched, not the whole ~60k-gene matrix
expr <- fetch_dense_values(host, expr_ds, identifiers = "TMEM173",
                           use_probeMap = TRUE)
expr <- data.frame(sample = colnames(expr), log2tpm001 = as.numeric(expr[1, ]))

pheno <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterDatasets = pheno_ds) %>%
  XenaQuery() %>%
  XenaDownload(destdir = cache) %>%
  XenaPrepare()
pheno <- pheno %>%
  rename(category = detailed_category, site = `_primary_site`,
         sample_type = `_sample_type`, study = `_study`) %>%
  select(sample, category, site, sample_type, study)

# Toil values are log2(TPM + 0.001); convert to log2(TPM + 1) like GEPIA
dat <- expr %>%
  inner_join(pheno, by = "sample") %>%
  mutate(expression = log2(pmax(2^log2tpm001 - 0.001, 0) + 1))

# ---- 2. Tumor / normal definitions ------------------------------------------
# Normal = GTEx tissue + TCGA adjacent "Solid Tissue Normal" (as in GEPIA).
# Edit gtex_site if you want a different normal comparator.
# AML: GTEx has no bone marrow; whole blood is used here -- consider dropping
# AML, since the argument is about solid tumors.
map <- tribble(
  ~cancer, ~tcga_category,                        ~gtex_site,
  "AML",   "Acute Myeloid Leukemia",              "Blood",
  "BRCA",  "Breast Invasive Carcinoma",           "Breast",
  "COAD",  "Colon Adenocarcinoma",                "Colon",
  "GBM",   "Glioblastoma Multiforme",             "Brain",
  "KIRP",  "Kidney Papillary Cell Carcinoma",     "Kidney",
  "LUSC",  "Lung Squamous Cell Carcinoma",        "Lung",
  "OV",    "Ovarian Serous Cystadenocarcinoma",   "Ovary",
  "PRAD",  "Prostate Adenocarcinoma",             "Prostate",
  "UCS",   "Uterine Carcinosarcoma",              "Uterus"
)
tumor_types <- c("Primary Tumor",
                 "Primary Blood Derived Cancer - Peripheral Blood")

pan <- bind_rows(lapply(seq_len(nrow(map)), function(i) {
  m <- map[i, ]
  tum <- dat %>% filter(study == "TCGA", category == m$tcga_category,
                        sample_type %in% tumor_types)
  nor <- dat %>% filter((study == "GTEX" & site == m$gtex_site) |
                        (study == "TCGA" & category == m$tcga_category &
                         sample_type == "Solid Tissue Normal"))
  if (nrow(tum) == 0 || nrow(nor) == 0) {
    print(sort(unique(dat$category))); print(sort(unique(dat$site)))
    stop("No samples matched for ", m$cancer, " -- check names above")
  }
  bind_rows(mutate(tum, cancer = m$cancer, group = "Tumor"),
            mutate(nor, cancer = m$cancer, group = "Normal"))
}))
pan$cancer <- factor(pan$cancer, levels = map$cancer)
pan$group  <- factor(pan$group,  levels = c("Tumor", "Normal"))

# OV panel: TCGA-OV vs GTEx fallopian tube (cell of origin) and GTEx ovary
ov <- bind_rows(
  dat %>% filter(study == "GTEX", site == "Fallopian Tube") %>%
    mutate(group = "Normal FT\n(GTEx)"),
  dat %>% filter(study == "GTEX", site == "Ovary") %>%
    mutate(group = "Normal ovary\n(GTEx)"),
  dat %>% filter(study == "TCGA",
                 category == "Ovarian Serous Cystadenocarcinoma",
                 sample_type == "Primary Tumor") %>%
    mutate(group = "HGSC\n(TCGA-OV)")
)
if (n_distinct(ov$group) < 3) stop("OV panel: a group has no samples")
ov$group <- factor(ov$group, levels = unique(ov$group))
print(table(ov$group))   # report the (small) fallopian tube n in the legend

# ---- 3. Statistics ----------------------------------------------------------
stars <- function(p) as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, .05, Inf),
                                      labels = c("****", "***", "**", "*", "ns")))

summ <- pan %>% group_by(cancer, group) %>%
  summarise(n = n(), mean = mean(expression), sd = sd(expression),
            median = median(expression), .groups = "drop")

pan_stats <- pan %>% group_by(cancer) %>%
  summarise(p = wilcox.test(expression[group == "Tumor"],
                            expression[group == "Normal"])$p.value,
            log2FC = mean(expression[group == "Tumor"]) -
                     mean(expression[group == "Normal"]),
            .groups = "drop") %>%
  mutate(FDR = p.adjust(p, "BH"), lab = stars(FDR))

stats_out <- summ %>%
  pivot_wider(names_from = group, values_from = c(n, mean, sd, median)) %>%
  left_join(pan_stats %>% select(-lab), by = "cancer")
write.csv(stats_out, file.path(outdir, "STING1_toil_stats.csv"), row.names = FALSE)
write.csv(bind_rows(pan %>% mutate(panel = "pancancer"),
                    ov %>% mutate(panel = "OV_FT", cancer = "OV")) %>%
            select(panel, cancer, group, sample, study, site, category,
                   sample_type, expression),
          file.path(outdir, "STING1_toil_tidy.csv"), row.names = FALSE)
print(as.data.frame(stats_out))

# ---- 4. Pan-cancer bar graph ------------------------------------------------
cols  <- c(Tumor = "#EABF00", Normal = "#2874C5")
dodge <- position_dodge(width = .8)
ntab  <- summ %>% select(cancer, group, n) %>%
  pivot_wider(names_from = group, values_from = n)
xlabs <- setNames(sprintf("%s\nT=%d\nN=%d", ntab$cancer, ntab$Tumor, ntab$Normal),
                  ntab$cancer)
ymax  <- max(pan$expression)
pan_stats$y <- ymax * 1.05

p1 <- ggplot() +
  geom_col(data = summ, aes(cancer, mean, fill = group),
           position = dodge, width = .75, color = "black", linewidth = .4,
           alpha = .85) +
  geom_point(data = pan, aes(cancer, expression, fill = group),
             shape = 21, size = .6, stroke = .1, alpha = .35, color = "black",
             position = position_jitterdodge(jitter.width = .3, dodge.width = .8)) +
  geom_errorbar(data = summ, aes(cancer, ymin = mean - sd, ymax = mean + sd,
                                 group = group),
                position = dodge, width = .25, linewidth = .5) +
  geom_text(data = pan_stats, aes(cancer, y, label = lab), size = 4) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = xlabs) +
  scale_y_continuous(expand = expansion(mult = c(0, .08))) +
  labs(x = NULL, y = expression("STING1 expression  " * log[2] * "(TPM+1)")) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "top",
        legend.justification = "right",
        axis.text = element_text(size = 9, color = "black"),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 12))

ggsave(file.path(outdir, "Fig_STING1_pancancer_bar.pdf"), p1, width = 11, height = 4.8)
ggsave(file.path(outdir, "Fig_STING1_pancancer_bar.png"), p1, width = 11, height = 4.8, dpi = 300)

# ---- 5. OV vs fallopian tube bar graph --------------------------------------
ov_summ <- ov %>% group_by(group) %>%
  summarise(n = n(), mean = mean(expression), sd = sd(expression), .groups = "drop")
tum_lab <- levels(ov$group)[3]
ov_p <- sapply(levels(ov$group)[1:2], function(g)
  wilcox.test(ov$expression[ov$group == g],
              ov$expression[ov$group == tum_lab])$p.value)
top <- max(ov$expression)
brk <- data.frame(x = c(1, 2), xend = 3, y = top * c(1.14, 1.06),
                  lab = sprintf("%s (p = %.2g)", stars(ov_p), ov_p))
ov_xlabs <- setNames(sprintf("%s\nn=%d", ov_summ$group, ov_summ$n), ov_summ$group)

p2 <- ggplot() +
  geom_col(data = ov_summ, aes(group, mean, fill = group), width = .65,
           color = "black", linewidth = .4, alpha = .85) +
  geom_jitter(data = ov, aes(group, expression, fill = group), shape = 21,
              size = 1.2, stroke = .15, alpha = .5, width = .18, height = 0) +
  geom_errorbar(data = ov_summ, aes(group, ymin = mean - sd, ymax = mean + sd),
                width = .2, linewidth = .5) +
  geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y), linewidth = .4) +
  geom_text(data = brk, aes(x = (x + xend) / 2, y = y, label = lab),
            vjust = -.4, size = 3.5) +
  scale_fill_manual(values = c("#2874C5", "#7FA7D9", "#EABF00")) +
  scale_x_discrete(labels = ov_xlabs) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(x = NULL, y = expression("STING1 expression  " * log[2] * "(TPM+1)")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12))

ggsave(file.path(outdir, "Fig_STING1_OV_FT_bar.pdf"), p2, width = 4.2, height = 4.8)
ggsave(file.path(outdir, "Fig_STING1_OV_FT_bar.png"), p2, width = 4.2, height = 4.8, dpi = 300)
message("done")
