library(magrittr)

contamarkov <- function(sample_table, reports, fdr_threshold=.1) {
  sample_table %>%
    dplyr::mutate(host_reads=total_reads-nonhost_reads-total_ercc_reads) %>%
    dplyr::mutate(host_concentration=host_reads/total_ercc_reads,
                  nonhost_concentration=nonhost_reads/total_ercc_reads) %>%
    dplyr::mutate(total_sample_concentration=(total_reads-total_ercc_reads) / total_ercc_reads *
                    ercc_concentration)->
    sample_table

  sample_table %>%
    dplyr::filter(!is_water) %>%
    with(sum(host_concentration)) ->
    total_host_mass

  sample_table %>%
    dplyr::filter(is_water) %>%
    with(mean(host_concentration)) / total_host_mass ->
    mean_spill_frac

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(!is_water) %>%
    dplyr::group_by(tax_id) %>%
    dplyr::summarize(spillome_concentration=sum(NT_r / total_ercc_reads * ercc_concentration) *
                       mean_spill_frac,
                     category_name=unique(category_name)) ->
    spillome

  n_water <- sum(sample_table$is_water)
  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(is_water) %>%
    dplyr::group_by(tax_id) %>%
    dplyr::summarize(labome_concentration=sum(NT_r / total_ercc_reads * ercc_concentration) / n_water, category_name=unique(category_name)) ->
    labome

  dplyr::full_join(spillome, labome) %>%
    tidyr::replace_na(list(spillome_concentration=0, labome_concentration=0)) %>%
    dplyr::mutate(contaminome_concentration=spillome_concentration + labome_concentration) ->
    contaminome

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::mutate(tax_concentration=NT_r / total_ercc_reads * ercc_concentration) %>%
    dplyr::left_join(contaminome) %>%
    dplyr::mutate(pval=pmin(1, contaminome_concentration / tax_concentration)) %>%
    dplyr::group_by(category_name) %>%
    dplyr::mutate(padj=p.adjust(pval, method="BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!! c(names(reports), "tax_concentration",
                        "total_sample_concentration", "pval", "padj")) ->
    reports

  reports %>%
    dplyr::filter(padj <= .1) %>%
    dplyr::group_by(category_name) %>%
    dplyr::summarize(threshold=max(pval)) ->
    pval_thresholds

  contaminome %>%
    dplyr::rename(spillome=spillome_concentration, labome=labome_concentration) %>%
    dplyr::select(tax_id, category_name, spillome, labome, contaminome_concentration) %>%
    dplyr::inner_join(pval_thresholds) %>%
    dplyr::mutate(threshold=contaminome_concentration/threshold) %>%
    dplyr::select(-contaminome_concentration) %>%
    dplyr::rename(cutoff=threshold) %>%
    dplyr::rename(spillome_expected=spillome, labome_expected=labome) %>%
    tidyr::gather(line, concentration, -tax_id, -category_name, -tax_id) %>%
    dplyr::filter(concentration > 0) %>%
    dplyr::mutate(log10_rpm_intercept=log10(concentration*1e6)) %>%
    dplyr::select(tax_id, line, log10_rpm_intercept) ->
    log10_rpm_intercept

  list(log10_rpm_intercept=log10_rpm_intercept, reports=reports)
}

