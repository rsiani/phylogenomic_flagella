### Comparative phylogenomics of 1839 Pseudomonadota flagellar assembly pathway
## Roberto Siani, roberto.siani@helmholtz-munich.de
# 2024

source("scripts/helper.R")

# metadata ----------------------------------------------------------------

# load manually curate metadata along with CDS count from own prediction

metadata <-
  read_tsv("data/meta.tsv") |>
  left_join(read_delim("data/gene_count.fwf",
    col_names = c("n_genes", "Accession"),
    delim = " "
  )) |>
  mutate(
    log_genome_size = log10(Total.Sequence.Length),
    log_num_genes = log10(n_genes)
  )


metadata |>
  summarise(across(where(is.character), ~ n_distinct(.x)))

count(metadata, set)

write_csv(metadata, "figures/fig_tab/SupplementaryData1.csv")

# fig.1a

metadata |>
  count(Class, set) |>
  pivot_wider(names_from = set, values_from = n, values_fill = 0) |>
  gt::gt(row_group_as_column = T) |>
  gt::grand_summary_rows(
    columns = c(free_living, host_associated),
    fns = list(total = ~ sum(.))
  )

# load tree from gtotree

meta_tree <-
  ape::read.tree("data/proto_out/proto_out.tre")

# safety checks

ape::is.binary(meta_tree)
ape::is.rooted(meta_tree)
ape::is.ultrametric(meta_tree)

# tree is made ultrametric, dichotomous and rooted to comply with phylogenetic
# requirements

rooted_tree <-
  meta_tree |>
  phytools::force.ultrametric(method = "extend") |>
  ape::multi2di() |>
  phangorn::midpoint()

ape::is.ultrametric(rooted_tree)
ape::is.binary(rooted_tree)
ape::is.rooted(rooted_tree)

# summary statistics

janitor::tabyl(metadata, set)

metadata |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T)),
    n = n(),
    .by = c(Class, set)
  )

# prevalence ---------------------------------------------------------

# KO definitions

ko_meta <-
  read_tsv("data/proto_out/KO_search_results/target-KOs.tsv") |>
  select(domain_name = knum, definition) |>
  mutate(
    Protein = case_match(
      definition,
      "flagellar hook-associated protein 1" ~ "FlgK",
      "flagellin" ~ "FliC",
      "flagellar protein FliO/FliZ" ~ "FliO/FliZ",
      "flagellar biosynthesis protein FliR/FlhB" ~ "FliR/FlhB",
      "flagellar hook-associated protein 2" ~ "FliD",
      "flagellum-specific ATP synthase [EC:7.4.2.8]" ~ "FliI",
      "L-cystine transport system substrate-binding protein" ~ "FliY",
      "RNA polymerase primary sigma factor" ~ "RpoD",
      "RNA polymerase sigma-54 factor" ~ "RpoN",
      "sigma-54 dependent transcriptional regulator, flagellar regulatory protein" ~
        "FlrA/FleQ/FlaK",
      .default = protein_id(definition)
    ),
    Operon =
      case_when(
        str_detect(Protein, "FlhC|FlhD") ~ "flhDC",
        str_detect(Protein, "FlgA|FlgM|FlgN") ~ "flgANM",
        str_detect(
          Protein,
          "FlgB|FlgC|FlgD|FlgE|FlgF|FlgG|FlgH|FlgI|FlgJ|FlgK|FlgL"
        ) ~
          "flgBCDEFGHIJKL",
        str_detect(Protein, "FlhB|FlhA|FlhE") ~ "flhBAE",
        str_detect(Protein, "FliA|FliZ|FliY") ~ "fliAZY",
        str_detect(Protein, "FliD|FliS|FliT") ~ "fliDST",
        str_detect(Protein, "FliE") ~ "fliE",
        str_detect(Protein, "FliF|FliG|FliH|FliI|FliJ|FliK") ~ "fliFGHIJK",
        str_detect(Protein, "FliL|FliM|FliN|FliO|FliP|FliQ|FliR") ~ "fliLMNOPQR",
        str_detect(Protein, "FlgM|FlgN") ~ "flgMN",
        str_detect(Protein, "FlgK|FlgL") ~ "flgKL",
        str_detect(Protein, "FliC") ~ "fliC",
        str_detect(Protein, "FljB|FljA") ~ "fljBA",
        str_detect(Protein, "MotA|MotB") ~ "motAB",
        str_detect(Protein, "MotC|MotD") ~ "motCD",
        str_detect(Protein, "MotX|MotY") ~ "motXY",
        str_detect(Protein, "FlgO|FlgP|FlgQ|FlgT") ~ "flgOPQT",
        .default = protein_to_gene(Protein)
      ),
    Role_II = case_when(
      str_detect(Protein, "RpoD|RpoN") ~ "Regulator: Sigma factor",
      str_detect(Protein, "FlhC|FlhD") ~ "Regulator: Master",
      str_detect(Protein, "FlrA/FleQ/FlaK|FlrC") ~ "Regulator: Alternative",
      str_detect(Protein, "FliO/FliZ|FliZ|FliY") ~ "Regulator: Accessory",
      str_detect(Protein, "FliA") ~ "Regulator: Flagellar-specific sigma factor",
      str_detect(Protein, "FlgM") ~ "Regulator: Anti-sigma factor",
      str_detect(Protein, "FliK") ~ "Regulator: Hook-lenght control",
      str_detect(Protein, "FlhA|FlhB|FliP|FliQ|FliR") ~ "fT3SS: Export apparatus",
      str_detect(Protein, "FliH|FliI|FliJ") ~ "fT3SS: ATPase and accessory",
      str_detect(Protein, "FlgN") ~ "fT3SS: Hook-assembly chaperon",
      str_detect(Protein, "FlhE") ~ "fT3SS: Export regulator",
      str_detect(Protein, "FliF") ~ "Switch: MS ring",
      str_detect(Protein, "FliG|FliM|FliN") ~ "Switch: C ring",
      str_detect(Protein, "FliL") ~ "Switch: Stabilizer",
      str_detect(Protein, "MotC|MotD|MotX|MotY") ~ "Motor: Motor-associated",
      str_detect(Protein, "MotA|MotB") ~ "Motor: Stator",
      str_detect(Protein, "FliC") ~ "Filament: Flagellin",
      str_detect(Protein, "FliD") ~ "Filament: Cap",
      str_detect(Protein, "FliS|FliT") ~ "Filament: Chaperon",
      str_detect(Protein, "FliE") ~ "Basal-body/Hook: MS-rod junction",
      str_detect(Protein, "FlgT|FlgO|FlgP|FlgQ") ~ "Basal-body/Hook: Accessory",
      str_detect(Protein, "FlgK|FlgL") ~ "Basal-body/Hook: Hook-filament junction",
      str_detect(Protein, "FlgJ") ~ "Basal-body/Hook: Peptidoglycan hydrolysis",
      str_detect(Protein, "FlgH|FlgI") ~ "Basal-body/Hook: L/P ring",
      str_detect(Protein, "FlgE") ~ "Basal-body/Hook: Hook",
      str_detect(Protein, "FlgD") ~ "Basal-body/Hook: Hook-capping",
      str_detect(Protein, "FlgA|FlgB|FlgC|FlgF|FlgG") ~ "Basal-body/Hook: Rod"
    ),
    Tier =
      case_when(
        Operon %in% c("flhDC", "flrA/FleQ/FlaK") ~ "I",
        Operon %in% c("fliE", "flhBAE", "fliFGHIJK", "fliLMNOPQR") ~ "II",
        Operon %in% c("fliAZY", "fliDST", "flgANM", "flgBCDEFGHIJKL", "flrC") ~ "II/III",
        Operon %in% c("rpoN", "rpoD") ~ "I+",
        .default = "III"
      ) |> factor(levels = c("I", "I+", "II", "II/III", "III")),
    Role_I =
      str_remove(Role_II, ":.*") |>
        factor(levels = c(
          "Regulator", "fT3SS", "Switch",
          "Basal-body/Hook", "Motor", "Filament"
        )),
    Gene = protein_to_gene(Protein),
    Gene = case_match(Gene, "flrA/FleQ/FlaK" ~ "flrA/fleQ/flaK",
      "fliO/FliZ" ~ "fliO/fliZ",
      "fliR/FlhB" ~ "fliR/flhB",
      .default = Gene
    )
  ) |>
  arrange(Tier, Role_I, Role_II, Protein) |>
  mutate(Operon = as.factor(Operon) |> fct_inorder())

# prevalence and redundancy data

ko_list <-
  read_tsv("data/proto_out/KO_search_results/KO-hit-counts.tsv") |>
  pivot_longer(
    cols = -c(assembly_id, total_gene_count),
    names_to = "domain_name", values_to = "n_copies"
  ) |>
  mutate(pa = if_else(n_copies > 0, 1, 0)) |>
  left_join(ko_meta) |>
  rename(Accession = assembly_id) |>
  left_join(metadata) |>
  filter(any(pa > 0), .by = Gene)


gene_prevalence <- ko_list |>
  filter(any(pa > 0), .by = Accession) |>
  filter(any(pa > 0), .by = Gene) |>
  mutate(flagellated = if_else(str_detect(Gene, "fliC") & pa == 1, "yes", "no")) |>
  mutate(flagellated = any(flagellated == "yes"), .by = Accession) |>
  filter(flagellated, .by = Accession) |>
  summarise(
    pa = sum(pa) / n(),
    .by = c(Gene, Role_II, Role_I, Operon, Tier, definition)
  ) |>
  arrange(-pa)


gene_prevalence |>
  write_csv("figures/fig_tab/SupplementaryData2.csv")

ggplot(gene_prevalence) +
  geom_histogram(aes(pa)) +
  scale_x_continuous(n.breaks = 20) +
  scale_y_continuous(n.breaks = 20)

rare <- gene_prevalence |>
  filter(pa <= 0.5) |>
  pull(Gene)
accessory <- gene_prevalence |>
  filter(pa > 0.5, pa <= 0.9) |>
  pull(Gene)
core <- gene_prevalence |>
  filter(pa > 0.9) |>
  pull(Gene)


# prevalence -------------------------------------------------------------------


safe_phyloglm <- safely(phylolm::phyloglm)

# iterate over each ortholog and fit model

options(future.globals.maxSize = 1 * 1e10)

res_fit <-
  ko_list |>
  mutate(lGS = log_genome_size - mean(log_genome_size), .by = Gene) |>
  nest(.by = c(Gene, Protein, domain_name, definition, Tier, Role_I, Role_II, Operon)) |>
  mutate(
    data2 = map(data, ~ column_to_rownames(.x, "Accession")),
    fit_logit = map(
      data2,
      ~ safe_phyloglm(
        pa ~ set * lGS,
        data = .x,
        btol = 100,
        log.alpha.bound = 10,
        boot = 999,
        phy = rooted_tree |> ape::keep.tip(phy = _, tip = rownames(.x)),
        method = "logistic_MPLE"
      )
    ),
    fit_pois = map(
      data2,
      ~ safe_phyloglm(
        n_copies ~ set * lGS,
        data = .x,
        boot = 999,
        phy = rooted_tree |> ape::keep.tip(phy = _, tip = rownames(.x)),
        method = "poisson_GEE"
      )
    )
  )

res_fit <- read_rds("bootstrap_phymodel.RDS")

# tidy up results

logit_res <-
  res_fit |>
  mutate(res_logit = map(
    fit_logit,
    ~ pluck(.x, "result") |>
      summary() |>
      pluck("coefficients") |>
      as_tibble(rownames = "term")
  )) |>
  select(-c(data, data2, fit_logit, fit_pois)) |>
  unnest(res_logit) |>
  mutate(
    term = str_remove(term, "set"),
    estimate = Estimate,
    std.error = StdErr,
    statistic = z.value,
    fdr(p.value, 0.05),
    .keep = "unused"
  ) |>
  filter(term != "(Intercept)")

logit_res |> count(term, significant)

filter(logit_res, term %in% "host_associated") |>
  arrange(Tier, Role_I, Role_II) |>
  select(Gene, lowerbootCI, estimate, upperbootCI, std.error, statistic, FDR) |>
  mutate(across(c(estimate, std.error, statistic), ~ round(.x, 2)),
    FDR = format(FDR, scientific = T, digits = 2)
  ) |>
  write_csv("figures/fig_tab/SupplementaryData3.csv")


summarise(ko_list,
  pa = sum(pa) / n(),
  .by = c(Gene, set, Tier, Role_I)
) |>
  ggplot(aes(x = pa, y = Gene)) +
  geom_path(
    aes(group = Gene),
    linewidth = 3,
    color = "#cccccc"
  ) +
  geom_point(
    aes(fill = set),
    shape = 21,
    size = 3,
    stroke = 1,
    color = "#555555"
  ) +
  geom_text(
    data = logit_res |>
      filter(term %in% "host_associated") |>
      mutate(FDR = format(FDR, scientific = T, digits = 2)),
    aes(
      label = if_else(significant, paste0(FDR, "*"), FDR),
      fontface = if_else(significant, "bold", "plain"),
      y = Gene,
      x = 1.1
    ),
    family = "Fira Sans",
    color = "#555555",
    hjust = 0.35
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic"),
    strip.text.y = element_text(
      angle = 0, hjust = 0,
      face = "bold"
    ),
    strip.placement = "outside",
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.major.y = element_line(
      color = "#555555",
      linewidth = 0.1
    ),
    panel.grid.major.x = element_line(
      color = "#555555",
      linewidth = 0.1
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_fill_manual(values = pal_set, labels = list(
    host_associated = "Host Associated",
    free_living = "Free Living"
  )) +
  facet_grid(Role_I ~ .,
    space = "free_y",
    scales = "free_y"
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0.1),
    name = "Prevalence (%)",
    breaks = seq(0, 1, 0.2)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 10)))

my_ggsave("fig_tab/fig1", w = 9, h = 10)

pois_res <-
  res_fit |>
  mutate(res_pois = map(
    fit_pois,
    ~ pluck(.x, "result") |>
      summary() |>
      pluck("coefficients") |>
      as_tibble(rownames = "term")
  )) |>
  select(-c(data, data2, fit_logit, fit_pois)) |>
  unnest(res_pois) |>
  mutate(
    term = str_remove(term, "set"),
    estimate = Estimate,
    std.error = StdErr,
    statistic = z.value,
    fdr(p.value, 0.05),
    .keep = "unused"
  ) |>
  filter(term != "(Intercept)")

pois_res |>
  count(term, significant)

filter(pois_res, term %in% "host_associated") |>
  arrange(Tier, Role_I, Role_II) |>
  select(Gene, estimate, std.error, statistic, FDR) |>
  mutate(across(c(estimate, std.error, statistic), ~ round(.x, 2)),
    FDR = format(FDR, scientific = T, digits = 2)
  ) |>
  write_csv("figures/fig_tab/SupplementaryData4.csv")

summarise(ko_list,
  multiple_copies = sum(n_copies > 1) / n(),
  .by = c(Gene, set)
) |>
  filter(multiple_copies > 0) |>
  ggplot(aes(x = multiple_copies, y = Gene)) +
  geom_col(
    aes(fill = set),
    color = "white",
    position = position_dodge()
  ) +
  geom_text(
    data = pois_res |>
      filter(term %in% "host_associated") |>
      drop_na(significant) |>
      mutate(FDR = format(FDR, scientific = T, digits = 2)),
    aes(
      label = if_else(significant, paste0(FDR, "*"), FDR),
      fontface = if_else(significant, "bold", "plain"),
      y = Gene,
      x = 0.6
    ),
    family = "Fira Sans",
    color = "#555555",
    hjust = 0.35
  ) +
  scale_fill_manual(
    values = pal_set,
    labels = list(
      "free_living" = "Free Living",
      "host_associated" = "Host Associated"
    )
  ) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  scale_x_continuous(
    expand = expansion(mult = 0.2),
    name = "Share of Genomes carrying Multiple Copies",
    breaks = seq(0, 6, 0.2)
  )

my_ggsave("fig_tab/suppFig1", 6, 8)

filter(pois_res, term == "host_associated", significant) |>
  View()

# bitscore ---------------------------------------------------------

hmmsearch_results <-
  read_csv(
    "data/proto_out/KO_search_results/hmmsearch_extended/hmmsearch_results.csv",
    col_names = c(
      "query_name", "domain_name",
      "sequence_evalue", "sequence_score", "sequence_bias",
      "domain_evalue", "domain_score", "domain_bias"
    )
  ) |>
  mutate(Accession = str_remove(query_name, "_[0-9]*$"))

ko_tblout <-
  hmmsearch_results |>
  left_join(metadata) |>
  filter(
    domain_score > 0,
    sequence_score > 0,
    sequence_score > 2 * sequence_bias,
    domain_score > 2 * domain_bias,
    sequence_score < 2 * domain_score,
  ) |>
  mutate(n_hits = n(), .by = c(Accession, domain_name)) |>
  slice_min(sequence_evalue, n = 1, by = c(Accession, domain_name), with_ties = F) |>
  mutate(
    dbs = (sequence_score - mean(sequence_score)) / sd(sequence_score),
    .by = domain_name
  ) |>
  left_join(ko_meta) |>
  mutate(lGS = log_genome_size - mean(log_genome_size), .by = Gene)

ko_tblout |>
  summarise(
    lGS = first(log_genome_size),
    dbs = mean(dbs),
    p_of_gene = n() / 53,
    n_hits = mean(n_hits) |>
      log(), .by = c(Accession, set)
  ) |>
  pivot_longer(c(n_hits, dbs)) |>
  ggplot(aes(lGS, value, color = set, size = p_of_gene)) +
  geom_point(alpha = 2 / 3, ) +
  scale_size(range = c(0.1, 2)) +
  geom_smooth(alpha = 1 / 3, fill = "#d3d3d3", method = "lm") +
  scale_color_manual(values = pal_set) +
  facet_wrap(~name, scales = "free_x")


safe_phylolm <- safely(.f = phylolm::phylolm)

res_fit_linear <-
  ko_tblout |>
  mutate(lGS = log_genome_size - mean(log_genome_size), .by = Gene) |>
  nest(.by = c(Gene, Protein, domain_name, definition, Tier, Role_I, Role_II, Operon)) |>
  mutate(
    data2 = map(data, ~ column_to_rownames(.x, "Accession")),
    fit_lm = map(
      data2,
      ~ safe_phylolm(
        dbs ~ set * lGS,
        data = .x,
        boot = 999,
        model = "BM",
        phy = rooted_tree |> ape::keep.tip(phy = _, tip = rownames(.x))
      )
    )
  )

saveRDS(res_fit_linear, "bootstrap_phymodel2.RDS")

res_fit_linear <- readRDS("bootstrap_phymodel2.RDS")

res_fit_tidy <-
  res_fit_linear |>
  mutate(res_lm = map(
    fit_lm,
    ~ pluck(.x, "result") |>
      summary() |>
      pluck("coefficients") |>
      as_tibble(rownames = "term")
  )) |>
  select(-c(data, data2, fit_lm)) |>
  unnest(res_lm) |>
  mutate(
    term = str_remove(term, "set"),
    estimate = Estimate,
    std.error = StdErr,
    statistic = t.value,
    fdr(p.value, 0.05),
    .keep = "unused"
  ) |>
  filter(term != "(Intercept)")


filter(res_fit_tidy, term %in% "host_associated") |>
  arrange(Tier, Role_I, Role_II) |>
  select(Gene, lowerbootCI, estimate, upperbootCI, std.error, statistic, FDR) |>
  mutate(across(c(estimate, std.error, statistic), ~ round(.x, 2)),
    FDR = format(FDR, scientific = T, digits = 2)
  ) |>
  write_csv("figures/fig_tab/SupplementaryData5.csv")


p3 <-
  res_fit_tidy |>
  filter(term != "lGS") |>
  mutate(text_xpos = max(upperbootCI)) |>
  ggplot() +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
    color = NA, fill = "#555555", alpha = .1
  ) +
  geom_pointrange(
    aes(
      xmin = lowerbootCI,
      x = estimate,
      xmax = upperbootCI,
      y = Gene,
      color = term,
      shape = term
    ),
    position = position_dodge(width = .5),
    size = 0.5,
    stroke = 0.5,
    linewidth = 0.5
  ) +
  geom_text(
    data = ~ filter(.x, term %in% "host_associated") |>
      mutate(FDR = format(FDR, scientific = T, digits = 2)),
    aes(
      label = if_else(significant, paste0(FDR, "*"), FDR),
      y = Gene,
      x = text_xpos,
      fontface = if_else(significant, "bold", "plain")
    ),
    family = "Fira Sans",
    color = "#555555",
    hjust = 0.35
  ) +
  facet_grid(Role_I ~ .,
    scales = "free",
    space = "free_y"
  ) +
  theme(
    axis.title.y = element_blank(),
    strip.text.y = element_text(
      angle = 0, hjust = 0,
      face = "bold"
    ),
    strip.placement = "outside",
    legend.position = "bottom",
    axis.line.y = element_blank(),
    axis.text.y = element_text(face = "italic"),
    axis.ticks.y = element_blank()
  ) +
  scale_color_manual(
    values = c("#CC79A7", "grey50", "grey75"),
    labels = list(
      host_associated = "Host Associated",
      `host_associated:lGS` = "Host Associated & Genome Size",
      lGS = "Genome Size"
    )
  ) +
  scale_shape_manual(
    values = c(19, 15, 17),
    labels = list(
      host_associated = "Host Associated",
      `host_associated:lGS` = "Host Associated & Genome Size",
      lGS = "Genome Size"
    )
  ) +
  guides(alpha = "none", color = guide_legend(override.aes = list(size = 2))) +
  coord_cartesian(clip = "off") +
  labs(x = "Estimated Effect on Sequence Bitscore")

p3

p4 <-
  res_fit_tidy |>
  ggplot() +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
    color = NA, fill = "#555555", alpha = .1,
  ) +
  geom_density(
    aes(x = estimate, fill = term),
    color = "white",
    alpha = .7,
    outline.type = "upper"
  ) +
  scale_fill_manual(values = c("#CC79A7", "grey50", "grey75")) +
  scale_y_continuous(n.breaks = 3) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )

p4

p4 + p3 +
  plot_layout(
    heights = c(1, 9),
    axes = "collect_x"
  )

my_ggsave("fig_tab/fig2", 9, 11)
