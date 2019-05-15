library(magrittr)
library(ggplot2)

# model counts as poisson-gamma (i.e. negative binomial)
empirical_bayes_inverse_rate <- function(counts) {
  # MOM estimator
  q_hat <- mean(counts) / var(counts)
  r_hat <- mean(counts) * q_hat / (1 - q_hat)

  alpha_prior <- r_hat
  beta_prior <- q_hat / (1 - q_hat)

  alpha_posterior <- alpha_prior + counts
  beta_posterior <- beta_prior + 1

  posterior_rate <- alpha_posterior / beta_posterior

  # NOTE requires alpha > 1, else infinity
  stopifnot(alpha_posterior > 1)
  beta_posterior / (alpha_posterior - 1)
}

contamarkov <- function(sample_table, reports, fdr_threshold=.1, use_empirical_bayes=FALSE) {
  if (use_empirical_bayes) {
    sample_table %>%
      dplyr::mutate(inv_ercc = empirical_bayes_inverse_rate(total_ercc_reads)) ->
      sample_table
  } else {
    sample_table %>%
      dplyr::mutate(inv_ercc = 1/total_ercc_reads) ->
      sample_table
  }

  sample_table %>%
    dplyr::mutate(host_reads=total_reads-nonhost_reads-total_ercc_reads) %>%
    dplyr::mutate(host_concentration=host_reads * inv_ercc,
                  nonhost_concentration=nonhost_reads * inv_ercc) %>%
    dplyr::mutate(total_sample_concentration=(total_reads-total_ercc_reads)  * inv_ercc *
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
    dplyr::group_by(tax_id, name, category_name) %>%
    dplyr::summarize(spillome_concentration=sum(NT_r  * inv_ercc * ercc_concentration) *
                       mean_spill_frac) %>%
    dplyr::ungroup() ->
    spillome

  n_water <- sum(sample_table$is_water)
  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(is_water) %>%
    dplyr::group_by(tax_id, name, category_name) %>%
    dplyr::summarize(labome_concentration=sum(NT_r  * inv_ercc * ercc_concentration) /
                       n_water) ->
    labome

  dplyr::full_join(spillome, labome) %>%
    tidyr::replace_na(list(spillome_concentration=0, labome_concentration=0)) %>%
    dplyr::mutate(contaminome_concentration=spillome_concentration + labome_concentration) ->
    contaminome

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::mutate(tax_concentration=NT_r  * inv_ercc * ercc_concentration) %>%
    dplyr::left_join(contaminome) %>%
    dplyr::mutate(pval=pmin(1, contaminome_concentration / tax_concentration)) %>%
    dplyr::group_by(category_name) %>%
    dplyr::mutate(padj=p.adjust(pval, method="BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!! c(names(reports), "is_water", "tax_concentration",
                        "total_sample_concentration", "pval", "padj")) ->
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
    dplyr::mutate(log10_rpm_intercept=log10(concentration*1e6)) %>%
    dplyr::select(tax_id, name, line, log10_rpm_intercept) ->
    log10_rpm_intercept

  list(log10_rpm_intercept=log10_rpm_intercept, reports=reports)
}

contamarkov_ggplot <- function(reports, log10_rpm_intercept, ...) {
  # filter log10_rpm_intercept to have same tax_id, name as reports
  reports %>%
    dplyr::select(tax_id, name) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(log10_rpm_intercept) ->
    log10_rpm_intercept

  ggplot(reports, aes(x=log10(total_sample_concentration), y=log10(NT_rpm), ...)) +
    geom_abline(data=log10_rpm_intercept, mapping=aes(slope=-1, intercept=log10_rpm_intercept, lty=line))
}

plot_contamarkov <- function(contamarkov_list, subset_taxa=NULL, point_aes=aes(color=is_water)) {
  log10_rpm_intercept <- contamarkov_list$log10_rpm_intercept
  reports <- contamarkov_list$reports

  if (!is.null(subset_taxa)) {
    log10_rpm_intercept %>%
      dplyr::filter(tax_id %in% subset_taxa) ->
      log10_rpm_intercept

    reports %>%
      dplyr::filter(tax_id %in% subset_taxa) ->
      reports
  }

  ggplot(reports, mapping=aes(x=log10(total_sample_concentration), y=log10(NT_rpm), group=tax_id)) +
    geom_point(mapping=point_aes) +
    geom_abline(data=log10_rpm_intercept,
                mapping=aes(slope=-1, intercept=log10_rpm_intercept, lty=line)) +
    facet_wrap(~name)
}
