library(magrittr)
library(ggplot2)

contamarkov <- function(sample_table, reports, fdr_threshold=.1,
                        use_nr=FALSE) {
  reports %>%
    dplyr::mutate(n_reads=if (use_nr) NR_r else NT_r) %>%
    {if (use_nr)
       dplyr::mutate(., n_reads=NR_r, reads_pct2=NR_r_pct)
      else
        dplyr::mutate(., n_reads=NT_r, reads_pct2=NT_r_pct)} %>%
    dplyr::inner_join(
      dplyr::select(sample_table, sample_name,
                    subsampled_fraction, nonhost_reads)
    ) %>%
    dplyr::mutate(reads_pct=NT_r / subsampled_fraction / nonhost_reads * 100) %>%
    dplyr::select(-subsampled_fraction, -nonhost_reads) ->
    reports

  reports %>%
    dplyr::filter(!is.na(reads_pct2)) %>%
    dplyr::arrange(desc(abs(log(reads_pct / reads_pct2)))) %>%
    {all.equal(round(.$reads_pct, 5), .$reads_pct2)} %>%
    stopifnot()

  sample_table %>%
    dplyr::filter(sample_name %in% unique(reports$sample_name)) %>%
    dplyr::mutate(
      host_transcriptome_reads = total_reads - reads_after_star - total_ercc_reads) %>%
    dplyr::mutate(
      host_transcriptome_concentration=(host_transcriptome_reads *
                                          ercc_concentration / total_ercc_reads),
      nonhost_concentration=(nonhost_reads * compression_ratio *
                               ercc_concentration / total_ercc_reads)) ->
    sample_table

  sample_table %>%
    dplyr::filter(!is_water) %>%
    with(sum(host_transcriptome_concentration)) ->
    total_host_transcriptome_mass

  sample_table %>%
    dplyr::filter(is_water) %>%
    with(mean(host_transcriptome_concentration)) / total_host_transcriptome_mass ->
    mean_spill_frac

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(!is_water) %>%
    dplyr::group_by(tax_id, name, category_name) %>%
    dplyr::summarize(spillome_concentration=mean_spill_frac * sum(
      n_reads * compression_ratio / subsampled_fraction *
        ercc_concentration / total_ercc_reads)) %>%
    dplyr::ungroup() ->
    spillome

  n_water <- sum(sample_table$is_water)
  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(is_water) %>%
    dplyr::group_by(tax_id, name, category_name) %>%
    dplyr::summarize(labome_concentration=sum(
      n_reads * compression_ratio / subsampled_fraction *
        ercc_concentration / total_ercc_reads) / n_water) ->
    labome

  dplyr::full_join(spillome, labome) %>%
    tidyr::replace_na(list(spillome_concentration=0, labome_concentration=0)) %>%
    dplyr::mutate(contaminome_concentration=spillome_concentration + labome_concentration) ->
    contaminome

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::mutate(tax_concentration=(
      n_reads * compression_ratio / subsampled_fraction *
        ercc_concentration / total_ercc_reads)) %>%
    dplyr::left_join(contaminome) %>%
    dplyr::mutate(pval=pmin(1, contaminome_concentration / tax_concentration)) %>%
    dplyr::group_by(category_name) %>%
    dplyr::mutate(padj=p.adjust(pval, method="BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!! c(names(reports), "is_water", "tax_concentration",
                        "contaminome_concentration", "labome_concentration", "spillome_concentration", "nonhost_concentration",
                        "pval", "padj")) ->
    reports

  reports %>%
    dplyr::filter(padj <= .1) %>%
    dplyr::group_by(category_name) %>%
    dplyr::summarize(threshold=max(pval)) ->
    pval_thresholds

  contaminome %>%
    dplyr::rename(spillome=spillome_concentration, labome=labome_concentration) %>%
    dplyr::select(tax_id, name, category_name, spillome, labome, contaminome_concentration) %>%
    dplyr::inner_join(pval_thresholds) %>%
    dplyr::mutate(threshold=contaminome_concentration/threshold) %>%
    dplyr::select(-contaminome_concentration) %>%
    dplyr::rename(cutoff=threshold) %>%
    dplyr::rename(spillome_expected=spillome, labome_expected=labome) %>%
    tidyr::gather(line, concentration, -tax_id, -name, -category_name, -tax_id) %>%
    dplyr::filter(concentration > 0) %>%
    dplyr::mutate(log10_pct_intercept=log10(concentration*100)) %>%
    dplyr::select(tax_id, name, line, log10_pct_intercept) ->
    log10_pct_intercept

  list(log10_pct_intercept=log10_pct_intercept, reports=reports)
}

contamarkov_ggplot <- function(reports, log10_pct_intercept, ...) {
  if (is.factor(reports$name)) {
    log10_pct_intercept %>%
      dplyr::mutate(name=as.character(name)) %>%
      dplyr::filter(name %in% unique(as.character(reports$name))) %>%
      dplyr::mutate(name=factor(name, levels=levels(reports$name))) ->
      log10_pct_intercept
  }

  # filter log10_pct_intercept to have same tax_id, name as reports
  reports %>%
    dplyr::select(tax_id, name) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(log10_pct_intercept) ->
    log10_pct_intercept

  ggplot(reports, aes(x=log10(nonhost_concentration), y=log10(reads_pct), ...)) +
    geom_abline(data=log10_pct_intercept, mapping=aes(slope=-1, intercept=log10_pct_intercept, lty=line))
}
