# Numeric Cox, copy-number, and survival summaries for the manuscript's Table 3.
source("00_config.R")
for (fp in c(F_EXPR, F_CLINICAL)) if (!file.exists(fp)) stop("Required input absent: ", fp)

pat12 <- function(x) substr(gsub("\\.", "-", as.character(x)), 1L, 12L)
is_primary <- function(x) grepl("^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-01", toupper(gsub("\\.", "-", as.character(x))))
zc <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x - mean(x, na.rm = TRUE)) / s }
number <- function(x) suppressWarnings(as.numeric(as.character(x)))
parse_os <- function(x) {
  s <- toupper(trimws(as.character(x)))
  fifelse(grepl("^1(:|$)|DECEASED|DEAD", s), 1L,
    fifelse(grepl("^0(:|$)|LIVING|ALIVE", s), 0L, NA_integer_))
}
require_columns <- function(d, cols, label) {
  missing <- setdiff(cols, names(d))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

# Sample-level cBioPortal expression, filtered to primary tumors before choosing one sample per patient.
expr_cols <- names(fread(F_EXPR, nrows = 0L))
if (!all(c("SAMPLE_ID", "VPS72", "MCL1", "SETDB1") %in% expr_cols))
  stop("Expression file missing SAMPLE_ID/VPS72/MCL1/SETDB1 columns")
expr <- fread(F_EXPR, select = c("SAMPLE_ID", "VPS72", "MCL1", "SETDB1"))
require_columns(expr, c("SAMPLE_ID", "VPS72", "MCL1", "SETDB1"), "Expression")
expr <- expr[is_primary(SAMPLE_ID)]
if (!nrow(expr)) stop("Expression export contains no TCGA primary-tumor (-01) sample barcodes")
expr[, `:=`(patient = pat12(SAMPLE_ID), sample_id = as.character(SAMPLE_ID),
  VPS72 = number(VPS72), MCL1 = number(MCL1), SETDB1 = number(SETDB1))]
expr <- expr[is.finite(VPS72)]
setorder(expr, patient, sample_id)
expr <- expr[!duplicated(patient), .(patient, expr_sample_id = sample_id, VPS72, MCL1, SETDB1)]

cl <- fread(F_CLINICAL, check.names = FALSE)
C_PAT <- "Patient ID"; C_SAMP <- "Sample ID"; C_AGE <- "Diagnosis Age"; C_SEX <- "Sex"
C_STAGE <- "Neoplasm Disease Stage American Joint Committee on Cancer Code"
C_ANEU <- "Aneuploidy Score"; C_FGA <- "Fraction Genome Altered"
C_TYPE <- "TCGA PanCanAtlas Cancer Type Acronym"; C_STYPE <- "Sample Type"
C_OSM <- "Overall Survival (Months)"; C_OSS <- "Overall Survival Status"
require_columns(cl, c(C_PAT, C_AGE, C_SEX, C_STAGE, C_ANEU, C_TYPE, C_STYPE, C_OSM, C_OSS), "Clinical")
if (!C_SAMP %in% names(cl)) cl[, (C_SAMP) := NA_character_]
if (!C_FGA %in% names(cl)) cl[, (C_FGA) := NA_real_]
cl <- cl[grepl("Primary", as.character(get(C_STYPE)), ignore.case = TRUE) |
  (is.na(get(C_STYPE)) & is_primary(get(C_SAMP)))]
if (!nrow(cl)) stop("No clinical primary-tumor records")
cov <- cl[, .(patient = pat12(get(C_PAT)), clinical_sample_id = as.character(get(C_SAMP)),
  cancer = as.character(get(C_TYPE)), age = number(get(C_AGE)), sex = as.character(get(C_SEX)),
  stage_raw = as.character(get(C_STAGE)), aneu = number(get(C_ANEU)), fga = number(get(C_FGA)),
  OS_months = number(get(C_OSM)), OS = parse_os(get(C_OSS)))]
cov <- cov[grepl("^TCGA-", patient) & !is.na(cancer) & !cancer %in% c("LAML", "DLBC")]
for (nm in c("cancer", "OS", "OS_months")) {
  clash <- cov[, .(nvals = uniqueN(get(nm)[!is.na(get(nm))])), by = patient][nvals > 1L]
  if (nrow(clash)) stop("Conflicting ", nm, " for patients: ", paste(head(clash$patient, 10), collapse = ", "))
}
setorder(cov, patient, clinical_sample_id)
cov <- cov[!duplicated(patient)]
cov[, sex := factor(fifelse(grepl("^m", sex, ignore.case = TRUE), "Male",
  fifelse(grepl("^f", sex, ignore.case = TRUE), "Female", NA_character_)))]
cov[, stage4 := factor(fifelse(grepl("^STAGE IV", toupper(stage_raw)), "IV",
  fifelse(grepl("^STAGE III", toupper(stage_raw)), "III",
  fifelse(grepl("^STAGE II", toupper(stage_raw)), "II",
  fifelse(grepl("^STAGE I", toupper(stage_raw)), "I", "Unknown")))),
  levels = c("I", "II", "III", "IV", "Unknown"))]

dat <- merge(expr, cov, by = "patient")
dat <- dat[is.finite(OS_months) & OS_months > 0 & !is.na(OS)]
if (!nrow(dat)) stop("No matched primary-tumor patients with valid survival records")

# Regional copy-number proxy: mean continuous log2 CN across the 20 neighbors (not arm-wide 1q gain).
dat[, `:=`(VPS72_cn = NA_real_, regional_1q21 = NA_real_)]
if (file.exists(F_CNV)) {
  header <- names(fread(F_CNV, nrows = 0L))
  missing_cn <- setdiff(c("SAMPLE_ID", "VPS72", NEIGHBORS), header)
  if (length(missing_cn)) stop("Copy-number file missing: ", paste(missing_cn, collapse = ", "))
  cn <- fread(F_CNV, select = c("SAMPLE_ID", "VPS72", NEIGHBORS))
  cn <- cn[is_primary(SAMPLE_ID)]
  cn[, patient := pat12(SAMPLE_ID)]
  cn[, (c("VPS72", NEIGHBORS)) := lapply(.SD, number), .SDcols = c("VPS72", NEIGHBORS)]
  cn[, regional_1q21 := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = NEIGHBORS]
  cn[!is.finite(regional_1q21), regional_1q21 := NA_real_]
  setorder(cn, patient, SAMPLE_ID)
  cn <- cn[!duplicated(patient), .(patient, VPS72_cn = VPS72, regional_1q21)]
  dat[, c("VPS72_cn", "regional_1q21") := NULL]
  dat <- merge(dat, cn, by = "patient", all.x = TRUE)
} else {
  message("Optional copy-number file not present; CN-adjusted models will be skipped: ", F_CNV)
}

# Optional independent proliferation signature.
dat[, prolif := NA_real_]
if (file.exists(F_PROLIF)) {
  sg <- fread(F_PROLIF)
  if ("SetName" %in% names(sg) && any(sg$SetName == "Module11_Prolif_score")) {
    r <- sg[SetName == "Module11_Prolif_score"]
    gene_cols <- setdiff(names(r), c("Source", "SetName"))
    if (nrow(r) != 1L) stop("Expected one Module11_Prolif_score row")
    p <- data.table(patient = pat12(gene_cols), prolif = number(unlist(r[, ..gene_cols], use.names = FALSE)))
    setorder(p, patient)
    p <- p[!duplicated(patient)]
    dat[, prolif := NULL]
    dat <- merge(dat, p, by = "patient", all.x = TRUE)
  }
}

dat[, `:=`(vps72_z = zc(VPS72), mcl1_z = zc(MCL1), setdb1_z = zc(SETDB1),
  reg_z = zc(regional_1q21), prolif_z = zc(prolif)), by = cancer]
dat[, VPS72_group := fifelse(VPS72 >= median(VPS72, na.rm = TRUE), "High", "Low"), by = cancer]
dat[, cn_group := fifelse(is.na(VPS72_cn), NA_character_,
  fifelse(VPS72_cn > CN_GAIN_CUT, "Gain", fifelse(VPS72_cn < CN_LOSS_CUT, "Loss", "Neutral")))]
dat <- dat[is.finite(vps72_z)]
if (!nrow(dat)) stop("No patients with nonconstant within-cancer VPS72 expression")
write_result(dat, "matched_patient_survival_data.tsv")

skips <- list()
record_skip <- function(label, reason) {
  skips[[length(skips) + 1L]] <<- data.table(Model = label, Reason = reason)
  message("Skipped ", label, ": ", reason)
  NULL
}
fit_one <- function(input, covariates, label, term = "vps72_z") {
  required <- c("OS_months", "OS", "cancer", term, covariates)
  if (!all(required %in% names(input))) return(record_skip(label, "Missing model columns"))
  d <- input[complete.cases(input[, ..required])]
  if (nrow(d) < MIN_MODEL_N || sum(d$OS) < MIN_MODEL_EVENTS)
    return(record_skip(label, "Insufficient complete cases or deaths"))
  predictors <- c(term, covariates)
  constant <- predictors[vapply(predictors, function(nm) uniqueN(d[[nm]]) < 2L, logical(1))]
  if (length(constant)) return(record_skip(label, paste("Constant covariates:", paste(constant, collapse = ", "))))
  strata_term <- if (uniqueN(d$cancer) > 1L) " + strata(cancer)" else ""
  form <- as.formula(paste0("Surv(OS_months, OS) ~ ", paste(predictors, collapse = " + "), strata_term))
  fit <- tryCatch(coxph(form, data = d), error = function(e) NULL)
  if (is.null(fit)) return(record_skip(label, "Cox model failed"))
  co <- coef(fit)
  if (!term %in% names(co) || !is.finite(co[[term]]))
    return(record_skip(label, "Exposure coefficient not estimable"))
  se <- sqrt(diag(vcov(fit)))[[term]]
  if (!is.finite(se) || se <= 0) return(record_skip(label, "Exposure standard error invalid"))
  p_wald <- 2 * pnorm(-abs(co[[term]] / se))
  base_terms <- if (length(covariates)) paste(covariates, collapse = " + ") else "1"
  f0 <- as.formula(paste0("Surv(OS_months, OS) ~ ", base_terms, strata_term))
  lrt <- tryCatch({
    m0 <- coxph(f0, data = d)
    df_diff <- attr(logLik(fit), "df") - attr(logLik(m0), "df")
    if (df_diff <= 0) NA_real_ else pchisq(2 * (as.numeric(logLik(fit)) - as.numeric(logLik(m0))),
      df = df_diff, lower.tail = FALSE)
  }, error = function(e) NA_real_)
  vif <- NA_real_
  if (length(covariates)) {
    lmfit <- tryCatch(lm(as.formula(paste(term, "~", paste(covariates, collapse = " + "))), data = d),
      error = function(e) NULL)
    if (!is.null(lmfit)) {
      rsq <- summary(lmfit)$r.squared
      if (is.finite(rsq)) vif <- if (rsq >= 1) Inf else 1 / (1 - rsq)
    }
  }
  ans <- data.table(Model = label, Exposure = term, Unit = if (term == "gain") "Gain vs neutral" else "Per 1 SD",
    N = nrow(d), Events = sum(d$OS), HR = exp(co[[term]]), CI_lower = exp(co[[term]] - qnorm(0.975) * se),
    CI_upper = exp(co[[term]] + qnorm(0.975) * se), P_Wald = p_wald, P_LRT = lrt,
    VIF = vif, Covariates = if (length(covariates)) paste(covariates, collapse = ", ") else "None")
  attr(ans, "patient_ids") <- d$patient
  ans
}
MODELS <- list(
  M1 = character(),
  M2 = c("age", "sex"),
  M3 = c("age", "sex", "stage4"),
  M4 = c("age", "sex", "stage4", "aneu", "fga"),
  M5 = c("mcl1_z", "setdb1_z"),
  M6 = c("reg_z"),
  M7 = c("mcl1_z", "setdb1_z", "reg_z"),
  M8 = c("prolif_z"),
  M9 = c("age", "sex", "stage4", "aneu", "fga", "mcl1_z", "setdb1_z", "reg_z")
)
main <- list(); same_sample <- list()
for (model in names(MODELS)) {
  fitted <- fit_one(dat, MODELS[[model]], model)
  if (is.null(fitted)) next
  main[[model]] <- fitted
  if (model != "M1") {
    sub <- dat[patient %in% attr(fitted, "patient_ids")]
    cmp <- fit_one(sub, character(), paste0(model, "_VPS72_unadjusted_same_patients"))
    if (!is.null(cmp)) {
      if (cmp$N != fitted$N) stop("Comparator sample mismatch in ", model)
      same_sample[[model]] <- cmp
    }
  }
}
strata <- list()
for (label in c("Gain", "Neutral")) {
  x <- fit_one(dat[cn_group == label], c("age", "sex"), paste0("S_", label, "_VPS72"))
  if (!is.null(x)) strata[[label]] <- x
}
if (any(dat$cn_group %in% c("Gain", "Neutral"), na.rm = TRUE)) {
  dcn <- copy(dat[cn_group %in% c("Gain", "Neutral")])
  dcn[, gain := as.integer(cn_group == "Gain")]
  x <- fit_one(dcn, c("age", "sex"), "S_Gain_vs_Neutral", term = "gain")
  if (!is.null(x)) strata[["CN_compare"]] <- x
}
per_cancer <- list()
for (ct in sort(unique(dat$cancer))) {
  x <- fit_one(dat[cancer == ct], c("age", "sex", "mcl1_z", "setdb1_z"), paste0("Cancer_", ct))
  if (!is.null(x)) { x[, Cancer := ct]; per_cancer[[ct]] <- x }
}
empty_cox <- data.table(Model = character(), Exposure = character(), Unit = character(),
  N = integer(), Events = integer(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric(),
  P_Wald = numeric(), P_LRT = numeric(), VIF = numeric(), Covariates = character())
combine_cox <- function(xs) if (length(xs)) rbindlist(xs, fill = TRUE) else copy(empty_cox)
main_dt <- combine_cox(main)
strata_dt <- combine_cox(strata)
per_dt <- combine_cox(per_cancer)
if (nrow(per_dt)) per_dt[, FDR := p.adjust(P_Wald, "BH")]
write_result(rbindlist(list(main_dt, strata_dt), fill = TRUE), "table3_cox_models.tsv")
write_result(combine_cox(same_sample), "cox_same_sample_comparators.tsv")
write_result(per_dt, "cox_per_cancer.tsv")
skip_dt <- if (length(skips)) rbindlist(skips, fill = TRUE) else data.table(Model = character(), Reason = character())
write_result(skip_dt, "model_skips.tsv")

write_result(dat[, .(Patients = .N, Deaths = sum(OS), Gained = sum(cn_group == "Gain", na.rm = TRUE),
  Median_OS_months = median(OS_months, na.rm = TRUE)), by = cancer][order(-Patients)], "cohort_summary.tsv")
mat <- dat[, .(VPS72 = vps72_z, MCL1 = mcl1_z, SETDB1 = setdb1_z, Regional_1q21_CN = reg_z)]
cor_mat <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "spearman"))
cor_long <- as.data.table(as.table(cor_mat))
setnames(cor_long, c("Variable_1", "Variable_2", "Spearman_rho"))
write_result(cor_long, "expression_cnv_correlations.tsv")

km_output <- function(data, column, comparison) {
  x <- data[!is.na(get(column)), .(OS_months, OS, Group = as.character(get(column)))]
  if (uniqueN(x$Group) < 2) return(NULL)
  fit <- survfit(Surv(OS_months, OS) ~ Group, data = x)
  s <- summary(fit)
  data.table(Comparison = comparison, Group = sub("^Group=", "", as.character(s$strata)),
    Time_months = s$time, N_risk = s$n.risk, N_event = s$n.event,
    Survival = s$surv, SE = s$std.err)
}
km_parts <- Filter(Negate(is.null), list(
  km_output(dat, "VPS72_group", "VPS72 within-cancer median split"),
  km_output(dat[cn_group %in% c("Neutral", "Gain")], "cn_group", "CN gain vs neutral")))
km <- if (length(km_parts)) rbindlist(km_parts, fill = TRUE) else
  data.table(Comparison = character(), Group = character(), Time_months = numeric(),
    N_risk = integer(), N_event = integer(), Survival = numeric(), SE = numeric())
write_result(km, "kaplan_meier_numerical_curves.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
