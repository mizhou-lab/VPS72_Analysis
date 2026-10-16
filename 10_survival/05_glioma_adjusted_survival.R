# Glioma (GBM and LGG): IDH, subtype, MGMT, grade, 1p/19q and proliferation Cox models.
# Adapted from VPS72_Glioma_Survival_Supplement.R; figure & workbook code removed.
# TCGA survival only; CPTAC can contribute descriptive IDH comparisons, not OS fits.
source("00_config.R"); need(c("survival","dplyr","tidyr","readr","tibble"))
suppressPackageStartupMessages({library(survival);library(dplyr);library(tidyr);library(readr);library(tibble)})
base_dir <- INPUT_DIR
f_expr <- list.files(base_dir,pattern="^VPS72__mRNA_Expression.*\\.txt$",full.names=TRUE)
if(!length(f_expr)) stop("No glioma VPS72__mRNA_Expression*.txt files under INPUT_DIR")
exclude_studies <- c(brain_cptac_2020="pediatric cohort")
study_priority <- c("gbm_tcga_gdc","lgg_tcga_gdc","gbm_tcga_pan_can_atlas_2018","lgg_tcga_pan_can_atlas_2018","brain_cptac_gdc")
f_clin <- required(GLIOMA_CLIN_FILE)
f_p9 <- input(GLIOMA_MUT_FILES[1]); f_p10 <- input(GLIOMA_MUT_FILES[2]); f_sig <- input(GLIOMA_PROLIF_FILE)
# ---- 1. Load -----------------------------------------------------------------
read_expr <- function(f) {
  x <- read_tsv(f, col_types = cols(.default = "c"))
  unit <- sub("^VPS72: mRNA Expression,? ?", "", names(x)[4])
  unit <- sub("^\\((.*)\\)$", "\\1", unit)            # "(FPKM)" -> "FPKM"
  x %>% setNames(c("study", "patient", "sample", "expr")) %>%
    mutate(expr = suppressWarnings(as.numeric(expr)), expr_unit = unit,
           source_file = basename(f))
}
expr <- bind_rows(lapply(f_expr, read_expr)) %>%
  filter(!is.na(expr), !study %in% names(exclude_studies))
message("Expression studies used: ",
        paste(sprintf("%s (n=%d)", names(table(expr$study)), table(expr$study)), collapse = ", "))

clin <- read_tsv(f_clin, col_types = cols(.default = "c"))

read_plot <- function(f) read_tsv(f, col_types = cols(.default = "c")) %>%
  setNames(c("sample", "IDH1_mutation", "VPS72_CNA")) %>% distinct(sample, .keep_all = TRUE)
mut_files <- c(f_p9,f_p10); mut_files <- mut_files[file.exists(mut_files)]
mut_cna <- if(length(mut_files)) bind_rows(lapply(mut_files,read_plot)) %>% distinct(sample,.keep_all=TRUE) else
  tibble(sample=character(),IDH1_mutation=character(),VPS72_CNA=character())

# TCGA pan-glioma reference (Ceccarelli et al. 2016): IDH, MGMT, grade, subtype
pub <- clin %>%
  filter(`Study ID` == "lgggbm_tcga_pub") %>%
  transmute(patient = `Patient ID`,
            IDH_pub = case_when(grepl("^(WT|IDHwt|Wildtype|Wild-type)$", `IDH status`, ignore.case=TRUE) ~ "Wildtype",
                                 grepl("mut", `IDH status`, ignore.case=TRUE) ~ "Mutant", TRUE ~ NA_character_),
            MGMT    = `MGMT promoter status`,
            grade   = `Neoplasm Histologic Grade`,
            subtype = `Original Subtype`,          # GBM: Verhaak subtypes incl. G-CIMP
            idh_codel = `IDH/codel subtype`,       # IDHwt / IDHmut-non-codel / IDHmut-codel
            codel   = case_when(`IDH/codel subtype` == "IDHmut-codel" ~ "Codeleted",
                                !is.na(`IDH/codel subtype`) ~ "Non-codeleted")) %>%
  distinct(patient, .keep_all = TRUE)

# Proliferation: Module11_Prolif_score, matched on 15-character TCGA barcode
prolif <- if (file.exists(f_sig)) {
  read_tsv(f_sig, show_col_types = FALSE) %>%
    filter(SetName == "Module11_Prolif_score") %>%
    select(-Source, -SetName) %>%
    pivot_longer(everything(), names_to = "barcode", values_to = "prolif") %>%
    mutate(sample15 = substr(barcode, 1, 15), prolif = as.numeric(prolif)) %>%
    filter(!is.na(prolif)) %>% distinct(sample15, .keep_all = TRUE) %>%
    select(sample15, prolif)
} else {
  message("Scores_160_Signatures.tsv not found: proliferation analyses skipped")
  tibble(sample15 = character(), prolif = numeric())
}

# ---- 2. One row per patient --------------------------------------------------
base <- expr %>%
  inner_join(clin %>% select(`Study ID`, `Sample ID`, OS_months = `OS Months`,
                             OS_status = `OS Status`, age = AGE, sex = Sex,
                             sample_type = `Sample Type`),
             by = c("study" = "Study ID", "sample" = "Sample ID")) %>%
  left_join(mut_cna, by = "sample") %>%
  mutate(OS_months = as.numeric(OS_months),
         event     = parse_event(OS_status),
         age       = as.numeric(age),
         log2expr  = log2(expr + 1),
         cohort    = ifelse(grepl("cptac", study), "CPTAC", "TCGA"),
         cancer    = ifelse(grepl("lgg", study, ignore.case = TRUE), "LGG", "GBM"),
         IDH1_call = case_when(IDH1_mutation == "No mutation" ~ "Wildtype",
                               !is.na(IDH1_mutation) & IDH1_mutation != "Not profiled" ~ "Mutant"))

tcga <- base %>%
  filter(cohort == "TCGA", is.na(sample_type) | sample_type %in% c("Primary Tumor", "Primary")) %>%
  arrange(match(study, study_priority), patient, sample) %>%
  distinct(patient, .keep_all = TRUE) %>%            # one study per patient
  left_join(pub, by = "patient") %>%
  mutate(IDH = coalesce(IDH_pub, IDH1_call),
         IDH_source = case_when(!is.na(IDH_pub) ~ "lgggbm_tcga_pub (IDH1/2)",
                                !is.na(IDH1_call) ~ "IDH1 mutation call"),
         sample15 = substr(sample, 1, 15), n_samples = 1L) %>%
  left_join(prolif, by = "sample15") %>% select(-sample15)

cptac <- base %>%
  filter(cohort == "CPTAC") %>%
  group_by(study, patient, cohort, cancer, expr_unit, source_file) %>%
  summarise(sample = paste(sort(unique(sample)), collapse = ";"), n_samples = n(),
            expr = mean(expr), log2expr = mean(log2expr),
            OS_months = first_nonmissing(OS_months), event = first_nonmissing(event),
            age = first_nonmissing(age), sex = first_nonmissing(sex),
            VPS72_CNA = first_nonmissing(VPS72_CNA),
            IDH1_mutation = paste(sort(unique(na.omit(IDH1_mutation))), collapse = ";"),
            IDH = case_when(any(IDH1_call == "Mutant", na.rm = TRUE) ~ "Mutant",
                            any(IDH1_call == "Wildtype", na.rm = TRUE) ~ "Wildtype"),
            .groups = "drop") %>%
  mutate(IDH_source = ifelse(is.na(IDH), NA, "IDH1 mutation call"))

# Cohorts usable for survival: censored patients must have follow-up time
surv_ok <- base %>% filter(cohort == "TCGA") %>% group_by(study) %>%
  summarise(ok = any(event == 0 & !is.na(OS_months)), .groups="drop") %>% filter(ok) %>% pull(study)

dat <- bind_rows(tcga, cptac) %>%
  group_by(study) %>%                      # scale within study (FPKM vs RSEM)
  mutate(vps72_z  = as.numeric(scale(log2expr)),
         high     = as.numeric(log2expr >= median(log2expr)),
         prolif_z = if (all(is.na(prolif))) NA_real_ else as.numeric(scale(prolif))) %>%
  ungroup() %>%
  mutate(VPS72_group = factor(ifelse(high == 1, "High", "Low"), levels = c("Low", "High")),
         IDH  = factor(IDH, levels = c("Wildtype", "Mutant")),
         sex  = factor(sex),
         MGMT = factor(MGMT, levels = c("Unmethylated", "Methylated")),
         in_survival = study %in% surv_ok & !is.na(OS_months) & !is.na(event) & OS_months > 0)

dsurv <- filter(dat, in_survival)
if(!nrow(dsurv)) stop("No eligible adult TCGA glioma survival records after data matching")
cancers <- intersect(c("LGG", "GBM"), unique(dsurv$cancer))   # LGG = panels A-F, GBM = G-L

message("Patients: ", paste(capture.output(print(count(dat, cancer, cohort, in_survival))), collapse = "\n"))

# ---- 3. Descriptive tests ----------------------------------------------------
desc_tab <- list()
for (ca in unique(dat$cancer)) {
  d <- filter(dat, cancer == ca)
  cos <- if (n_distinct(d$cohort) > 1) c("All cohorts", unique(d$cohort)) else "All cohorts"
  for (co in cos) {                                             # VPS72 by IDH
    x <- if (co == "All cohorts") d else filter(d, cohort == co)
    x <- filter(x, !is.na(IDH))
    if (n_distinct(x$IDH) == 2)
      desc_tab[[length(desc_tab) + 1]] <- tibble(
        cancer = ca, test = "VPS72 by IDH status (Wilcoxon)", subset = co,
        n = nrow(x), statistic = NA_real_,
        detail = sprintf("median z: wildtype %.2f (n=%d), mutant %.2f (n=%d)",
                         median(x$vps72_z[x$IDH == "Wildtype"]), sum(x$IDH == "Wildtype"),
                         median(x$vps72_z[x$IDH == "Mutant"]), sum(x$IDH == "Mutant")),
        p = wilcox.test(vps72_z ~ IDH, data = x)$p.value)
  }
  s <- filter(d, cohort == "TCGA", !is.na(subtype))              # VPS72 by subtype
  if (ca == "GBM" && n_distinct(s$subtype) > 1) {
    kw <- kruskal.test(vps72_z ~ factor(subtype), data = s)
    med <- s %>% group_by(subtype) %>% summarise(m = median(vps72_z), n = n(), .groups = "drop")
    desc_tab[[length(desc_tab) + 1]] <- tibble(
      cancer = ca, test = "VPS72 by transcriptional subtype (Kruskal-Wallis)", subset = "TCGA",
      n = nrow(s), statistic = unname(kw$statistic),
      detail = paste(sprintf("%s %.2f (n=%d)", med$subtype, med$m, med$n), collapse = "; "),
      p = kw$p.value)
  }
  s2 <- filter(d, cohort == "TCGA", !is.na(idh_codel))           # LGG: IDH/codel subtype
  if (ca == "LGG" && n_distinct(s2$idh_codel) > 1) {
    kw <- kruskal.test(vps72_z ~ factor(idh_codel), data = s2)
    med <- s2 %>% group_by(idh_codel) %>% summarise(m = median(vps72_z), n = n(), .groups = "drop")
    desc_tab[[length(desc_tab) + 1]] <- tibble(
      cancer = ca, test = "VPS72 by IDH/1p19q subtype (Kruskal-Wallis)", subset = "TCGA",
      n = nrow(s2), statistic = unname(kw$statistic),
      detail = paste(sprintf("%s %.2f (n=%d)", med$idh_codel, med$m, med$n), collapse = "; "),
      p = kw$p.value)
  }
  pr <- filter(d, !is.na(prolif))                                # VPS72 vs proliferation
  if (nrow(pr) > 10) {
    ct <- suppressWarnings(cor.test(pr$log2expr, pr$prolif, method = "spearman", exact = FALSE))
    desc_tab[[length(desc_tab) + 1]] <- tibble(
      cancer = ca, test = "VPS72 vs proliferation (Spearman)", subset = "TCGA",
      n = nrow(pr), statistic = unname(ct$estimate),
      detail = "statistic = Spearman rho; proliferation = Module11_Prolif_score", p = ct$p.value)
  }
}
desc_tab <- bind_rows(desc_tab)

# ---- 4. Cox models -----------------------------------------------------------
# Each model: covariates + patient subset. "Comparator" models refit IDH + age
# on the subset used by the model that adds a covariate.
model_specs <- function(ca) {
  cc  <- function(d) filter(d, !is.na(IDH), !is.na(age))
  sp <- list(
    list(id = "M1", label = "Unadjusted",                 cov = NULL,                  sub = function(d) d),
    list(id = "M2", label = "Unadjusted (complete cases)", cov = NULL,                  sub = cc),
    list(id = "M3", label = "+ IDH + age",                cov = c("IDH", "age"),       sub = cc),
    list(id = "M4", label = "+ IDH + age + sex",          cov = c("IDH", "age", "sex"), sub = cc),
    list(id = "M5", label = "IDH-wildtype only (+ age)",  cov = "age",
         sub = function(d) filter(cc(d), IDH == "Wildtype")),
    list(id = "M6c", label = "+ IDH + age [MGMT subset]", cov = c("IDH", "age"),
         sub = function(d) filter(cc(d), !is.na(MGMT))),
    list(id = "M6", label = "+ IDH + age + MGMT",         cov = c("IDH", "age", "MGMT"),
         sub = function(d) filter(cc(d), !is.na(MGMT))),
    list(id = "M8c", label = "+ IDH + age [proliferation subset]", cov = c("IDH", "age"),
         sub = function(d) filter(cc(d), !is.na(prolif_z))),
    list(id = "M8", label = "+ IDH + age + proliferation", cov = c("IDH", "age", "prolif_z"),
         sub = function(d) filter(cc(d), !is.na(prolif_z))))
  if (ca == "GBM") sp <- c(sp, list(
    list(id = "M7c", label = "+ IDH + age [subtype subset]", cov = c("IDH", "age"),
         sub = function(d) filter(cc(d), !is.na(subtype))),
    list(id = "M7", label = "+ IDH + age + subtype",      cov = c("IDH", "age", "subtype"),
         sub = function(d) filter(cc(d), !is.na(subtype)))))
  if (ca == "LGG") sp <- c(sp, list(
    list(id = "M9c", label = "+ IDH + age [grade/codel subset]", cov = c("IDH", "age"),
         sub = function(d) filter(cc(d), !is.na(grade), !is.na(codel))),
    list(id = "M9", label = "+ IDH + age + grade + 1p/19q", cov = c("IDH", "age", "grade", "codel"),
         sub = function(d) filter(cc(d), !is.na(grade), !is.na(codel)))))
  sp
}

fit_one <- function(d, spec, term) {
  d <- spec$sub(d)
  if (nrow(d) < 20 || sum(d$event, na.rm=TRUE) < 10) return(NULL)
  # drop covariates with a single level in this subset (e.g. IDH in IDH-wildtype)
  cov <- Filter(function(v) length(unique(na.omit(d[[v]]))) > 1, spec$cov)
  st  <- if (n_distinct(d$study) > 1) " + strata(study)" else ""
  rhs <- paste(c(term, cov), collapse = " + ")
  f   <- as.formula(paste("Surv(OS_months, event) ~", rhs, st))
  fit <- tryCatch(coxph(f, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  lr_p <- if (length(cov)) {                      # likelihood-ratio test for VPS72
    f0 <- as.formula(paste("Surv(OS_months, event) ~", paste(cov, collapse = " + "), st))
    tryCatch(anova(coxph(f0, data = d), fit)[2, "Pr(>|Chi|)"], error=function(e) NA_real_)
  } else s$logtest["pvalue"]
  tibble(model_id = spec$id, model = spec$label, term = term,
         n = fit$n, events = fit$nevent,
         HR = s$conf.int[term, 1], lower95 = s$conf.int[term, 3], upper95 = s$conf.int[term, 4],
         p_wald = s$coefficients[term, "Pr(>|z|)"], p_LRT_VPS72 = unname(lr_p),
         covariates = paste(cov, collapse = ", "),
         formula = paste(deparse(f, width.cutoff = 500), collapse = ""))
}

cox_tab <- bind_rows(lapply(cancers, function(ca) {
  d <- filter(dsurv, cancer == ca)
  bind_rows(lapply(model_specs(ca), function(sp)
    bind_rows(fit_one(d, sp, "vps72_z"), fit_one(d, sp, "high")))) %>%
    mutate(cancer = ca, .before = 1)
})) %>%
  mutate(VPS72_coding = ifelse(term == "vps72_z", "Per SD (primary)", "High vs low (median)"),
         .after = term)

# Proportional hazards (Schoenfeld) for the main adjusted model
ph_tab <- bind_rows(lapply(cancers, function(ca) {
  d <- filter(dsurv, cancer == ca, !is.na(IDH), !is.na(age))
  z <- tryCatch(cox.zph(coxph(Surv(OS_months, event) ~ vps72_z + IDH + age, data = d))$table,
                error=function(e) NULL)
  if(is.null(z)) return(tibble(cancer=character(),model=character(),variable=character(),chisq=numeric(),df=numeric(),p=numeric()))
  tibble(cancer = ca, model = "vps72_z + IDH + age", variable = rownames(z),
         chisq = z[, "chisq"], df = z[, "df"], p = z[, "p"])
}))

print(cox_tab %>% filter(term == "vps72_z") %>%
        transmute(cancer, model, n, events, HR = round(HR, 2),
                  CI = sprintf("%.2f-%.2f", lower95, upper95), p = signif(p_wald, 2)), n = 50)


# Numerical source tables only (no plots / Excel workbooks).
write_table(dat, "05_glioma_analysis_samples.tsv")
write_table(cox_tab, "05_glioma_adjusted_Cox_models.tsv")
write_table(desc_tab, "05_glioma_descriptive_tests.tsv")
write_table(ph_tab, "05_glioma_PH_diagnostics.tsv")
cohort_tab <- dat %>% group_by(cancer,cohort,study,in_survival) %>%
  summarise(patients=n(),events=sum(event==1,na.rm=TRUE),IDH_mutant=sum(IDH=="Mutant",na.rm=TRUE),
            IDH_wildtype=sum(IDH=="Wildtype",na.rm=TRUE),.groups="drop")
write_table(cohort_tab,"05_glioma_cohort_counts.tsv")
