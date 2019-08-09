library(magrittr)
library(ggplot2)

contamarkov <- function(sample_table, reports, fdr_threshold=.1,
                        use_nr=FALSE) {
  reports %>%
    dplyr::mutate(n_reads=if (use_nr) NR_r else NT_r) %>%
    {if (use_nr)
       dplyr::mutate(., n_reads=NR_r, rpm=NR_rpm)
      else
        dplyr::mutate(., n_reads=NT_r, rpm=NT_rpm)} %>%
    dplyr::inner_join(
      dplyr::select(
        sample_table, sample_name, compression_ratio)
    ) %>%
    dplyr::mutate(rpm_corrected=rpm * compression_ratio) %>%
    dplyr::select(-compression_ratio) ->
    reports

  reports %>%
    dplyr::inner_join(dplyr::select(sample_table, sample_name, subsampled_fraction,
                                    total_reads, total_ercc_reads)) %>%
    with(all.equal(rpm, 1e6* n_reads / subsampled_fraction / (total_reads-total_ercc_reads), tolerance=1e-5)) %>%
    stopifnot()

  sample_table %>%
    dplyr::filter(sample_name %in% unique(reports$sample_name)) %>%
    dplyr::mutate(
      concentration=ercc_concentration * (total_reads - total_ercc_reads) / total_ercc_reads) ->
    sample_table

  n_water <- sum(sample_table$is_water)
  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(is_water) %>%
    dplyr::group_by(tax_id, name, category_name) %>%
    dplyr::summarize(labome_concentration=sum(
      n_reads * compression_ratio / subsampled_fraction *
        ercc_concentration / total_ercc_reads) / n_water) ->
    labome

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::mutate(tax_concentration=(
      n_reads * compression_ratio / subsampled_fraction *
        ercc_concentration / total_ercc_reads)) %>%
    dplyr::left_join(labome) %>%
    dplyr::mutate(pval=pmin(1, labome_concentration / tax_concentration)) %>%
    dplyr::group_by(category_name) %>%
    dplyr::mutate(padj=p.adjust(pval, method="BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!! c(names(reports), "is_water", "tax_concentration",
                        "labome_concentration", "concentration",
                        "pval", "padj")) ->
    reports

  reports %>%
    dplyr::filter(padj <= .1) %>%
    dplyr::group_by(category_name) %>%
    dplyr::summarize(threshold=max(pval)) ->
    pval_thresholds

  labome %>%
    dplyr::rename(labome=labome_concentration) %>%
    dplyr::select(tax_id, name, category_name, labome) %>%
    dplyr::inner_join(pval_thresholds) %>%
    dplyr::mutate(threshold=labome/threshold) %>%
    dplyr::rename(cutoff=threshold) %>%
    dplyr::rename(labome_expected=labome) %>%
    tidyr::gather(line, concentration, -tax_id, -name, -category_name, -tax_id) %>%
    dplyr::filter(concentration > 0) %>%
    dplyr::mutate(log10_rpm_intercept=log10(concentration*1e6)) %>%
    dplyr::select(tax_id, name, line, log10_rpm_intercept) ->
    log10_rpm_intercept

  list(log10_rpm_intercept=log10_rpm_intercept, reports=reports)
}

contamarkov_ggplot <- function(reports, log10_rpm_intercept, ...) {
  if (is.factor(reports$name)) {
    log10_rpm_intercept %>%
      dplyr::mutate(name=as.character(name)) %>%
      dplyr::filter(name %in% unique(as.character(reports$name))) %>%
      dplyr::mutate(name=factor(name, levels=levels(reports$name))) ->
      log10_rpm_intercept
  }

  # filter log10_rpm_intercept to have same tax_id, name as reports
  reports %>%
    dplyr::select(tax_id, name) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(log10_rpm_intercept) ->
    log10_rpm_intercept

  ggplot(reports, aes(x=log10(concentration), y=log10(rpm_corrected), ...)) +
    geom_abline(data=log10_rpm_intercept, mapping=aes(slope=-1, intercept=log10_rpm_intercept, lty=line))
}
