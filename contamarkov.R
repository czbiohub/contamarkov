library(magrittr)

contamarkov <- function(sample_table, reports) {
  sample_table %>%
    dplyr::mutate(host_reads=total_reads-nonhost_reads-total_ercc_reads) %>%
    dplyr::mutate(host_mass=host_reads/total_ercc_reads,
                  nonhost_mass=nonhost_reads/total_ercc_reads) ->
    sample_table

  sample_table %>%
    dplyr::filter(!is_water) %>%
    with(sum(host_mass)) ->
    total_host_mass

  sample_table %>%
    dplyr::filter(is_water) %>%
    with(mean(host_mass)) / total_host_mass ->
    mean_spill_frac

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(!is_water) %>%
    dplyr::group_by(name) %>%
    dplyr::summarize(spillome_concentration=sum(NT_r / total_ercc_reads) *
                       mean_spill_frac) ->
    spillome

  n_water <- sum(sample_table$is_water)
  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::filter(is_water) %>%
    dplyr::group_by(name) %>%
    dplyr::summarize(labome_concentration=sum(NT_r / total_ercc_reads) / n_water) ->
    labome

  dplyr::full_join(spillome, labome) %>%
    tidyr::replace_na(list(spillome_concentration=0, labome_concentration=0)) %>%
    #dplyr::mutate(mass_mean=pmax(spillome_concentration, labome_concentration)) ->
    dplyr::mutate(mass_mean=spillome_concentration + labome_concentration) ->
    contaminome

  reports %>%
    dplyr::inner_join(sample_table) %>%
    dplyr::mutate(mass=NT_r / total_ercc_reads) %>%
    dplyr::left_join(contaminome) %>%
    dplyr::mutate(mass_corrected=mass-mass_mean) %>%
    dplyr::mutate(pval=pmin(1, mass_mean / mass)) %>%
    dplyr::group_by(category_name) %>%
    dplyr::mutate(padj=p.adjust(pval, method="BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample_name, name, NT_r, NT_rpm, category_name, total_ercc_reads,
                  is_water, LRTI_adjudication,
                  mass, mass_mean, mass_corrected,
                  spillome_concentration, labome_concentration,
                  pval, padj) ->
    cleaned_reports

  list(contaminome=contaminome, reports=cleaned_reports)
}
