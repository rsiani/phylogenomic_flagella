### Comparative phylogenomics of 1839 Pseudomonadota flagellar assembly pathway
## Roberto Siani, roberto.siani@helmholtz-munich.de
# 2024

source("scripts/helper.R")

# metadata ----------------------------------------------------------------

# load manually curate metadata along with CDS count from own prediction

metadata =
  read_tsv("data/meta.tsv") |>
  left_join(read_delim("data/gene_count.fwf",
                       col_names = c("n_genes", "Accession"),
                       delim = " ")) |>
  mutate(log_genome_size = log10(Total.Sequence.Length),
         log_num_genes = log10(n_genes),
         # pathogen = if_else(Organism.Name %in% (read_csv("pathogen_list.csv") |>
         #                                          pluck("species")), "P", "N") |>
         #   factor(levels = c("N", "P")),
         association =
           case_when(
             # pathogen == "P" ~ "pathogen",
             set == "free_living" ~ "free_living",
             str_detect(Organism.Name, "Buchnera|Wigglesworthia|Wolbachia|Carsonella|Portiera|Nardonella|Blochmannia|Baumannia|symbiotic|Fukatsuia|Coxiella|symbiont|Fokinia solitaria|Gromoviella agglomerans|Neorickettsia|Rickettsi|Fonsibacter ubiquis|Ehrlichia|Midichloria|Liberibacter|Anaplasma|Neoehrlichia|Bartonella|Nucleicultrix|Endolissoclinum|Phycorickettsia|Profftella|Kinetoplastibacterium|Vallotia|Azoamicus|endobia|Purcelliella|Schneideria|Riesia|Tachikawaea|Ishikawaella|Arsenophonus|Comchoanobacter|Vesicomyosocius|Ruthturnera|Ruthia|Steffania|Francisella|Annandia|Orientia|Paracaedibacter") &
               Organism.Name != "Photorhabdus asymbiotica" ~ "endosymbionts",
             str_detect(str_c(isolation_source,	env_broad_scale, env_local_scale,	host),
                        "symbiont") ~ "endosymbionts",
             .default = "host_associated") |>
           factor(levels = c("free_living", "host_associated", "endosymbionts")),
         association_degree = as.numeric(association) - 1,
         endosymbiont = if_else(association == "endosymbionts", "E", "N") |>
           factor(levels = c("N", "E")))


metadata |>
  summarise(across(where(is.character), ~ n_distinct(.x)))

count(metadata, set)

write_csv(metadata, "figures/fig_tab/SupplementaryData1.csv")

# fig.1a

metadata |>
  count(Class, association) |>
  pivot_wider(names_from = association, values_from = n, values_fill = 0) |>
  gt::gt(row_group_as_column = T) |>
  gt::grand_summary_rows(
    columns = c(free_living, host_associated, endosymbionts),
    fns = list(total = ~ sum(.)))

# load tree from gtotree

meta_tree =
  ape::read.tree("data/proto_out/proto_out.tre")

# safety checks

ape::is.binary(meta_tree)
ape::is.rooted(meta_tree)
ape::is.ultrametric(meta_tree)

# tree is made ultrametric, dichotomous and rooted to comply with phylogenetic
# requirements

rooted_tree =
  meta_tree |>
  phytools::force.ultrametric(method = "extend") |>
  ape::multi2di() |>
  phangorn::midpoint()

ape::is.ultrametric(rooted_tree)
ape::is.binary(rooted_tree)
ape::is.rooted(rooted_tree)

# summary statistics

janitor::tabyl(metadata, association)

metadata |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T)),
            n = n(),
            .by = c(Class, association))

# fig1a --------------------------------------------------------------------

tree_w_data =
  tidytree::as_tibble(rooted_tree) |>
  left_join(metadata, join_by(label == Accession)) |>
  tidytree::as.treedata()

ggtree::ggtree(tree_w_data,
               layout = "fan",
               branch.length = "none",
               aes(color = set),
               open.angle = 45,
               linewidth = .5,
               show.legend = T) +
  scale_color_manual(values = pal_set) +
  ggtreeExtra::geom_fruit(
    aes(fill = endosymbiont),
    geom = geom_tile,
    width =  3,
    show.legend = T, offset = 0.05) +
  scale_fill_manual(values = c("N" = "#555555",
                               "E" = "#AA4499")) +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom = geom_col,
    aes(x = Total.Sequence.Length,
        fill = Class,
        color = after_scale(fill)),
    pwidth = 0.3, show.legend = F, offset = 0.03,
    axis.params = list(
      axis = "x",
      text.angle = 90,
      text.size = 15/.pt,
      line.size = .5,
      hjust = 1
    )) +
  scale_fill_manual(values = palette_Compartment_Phylum) +
  theme(plot.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.position = "none")

my_ggsave("fig1a", 9, 9)


# fig.1b

metadata |>
  pivot_longer(c(log_genome_size, log_num_genes)) |>
  mutate(name = factor(name, levels = c("log_genome_size",
                                        "log_num_genes"))) |>
  ggplot(aes(x = Class, y = value,
             fill = association)) +
  scale_fill_manual(values = palette_Compartment_Phylum,
                    aesthetics = c("color", "fill")) +
  geom_point(
    shape = 21,
    position = ggbeeswarm::position_quasirandom(dodge.width = 1, width = .25),
    color = "#333333",
    alpha = .8,
    stroke = 1/3,
    size = 1.5) +
  facet_wrap(~ name, scales = "free_y",
             labeller = as_labeller(c(
               `GC.Percent` = "% GC",
               `log_genome_size` = "Genome Size",
               `log_num_genes` = "Number of Genes"))) +
  scale_x_discrete(labels = c(
    Alphaproteobacteria  = "&alpha;",
    Betaproteobacteria = "&beta;",
    Gammaproteobacteria = "&gamma;")) +
  scale_y_continuous(name = "log<sub>10</sub>") +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(angle = 15,
                                               hjust = .5,
                                               vjust = .5)) +
  guides(fill = guide_legend(override.aes = list(size = 5)))


my_ggsave("fig1b", w = 5, h = 3.5)

phylolm::phylolm(log_genome_size ~ association_degree,
                 data = metadata |>
                   column_to_rownames("Accession"),
                 phy = rooted_tree) |>
  summary() |>
  pluck("coefficients") |>
  as_tibble(rownames = "term")


# plot legend

data.frame(set = c("free_living", "host_associated", "endosymbionts",
                   "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria"),
           x = c(1, 1, 1, 2, 2, 2),
           y = c(1, 0.5, 0, 1, 0.5, 0)) |>
  ggplot() +
  geom_point(aes(x, y, color = set), size = 10, show.legend = F) +
  geom_text(aes(x, y, label = set), nudge_x = 0.1, hjust = 0, size = 15/.pt) +
  xlim(0.5, 3) +
  ylim(-1, 2) +
  scale_color_manual(values = palette_Compartment_Phylum) +
  theme_void()

# prevalence ---------------------------------------------------------

# KO definitions

ko_meta =
  read_tsv("data/proto_out/KO_search_results/target-KOs.tsv") |>
  select(domain_name = knum, definition) |>
  mutate(Protein = case_match(
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
    .default = protein_id(definition)),
    Operon =
      case_when(
        str_detect(Protein, "FlhC|FlhD") ~ "flhDC",
        str_detect(Protein, "FlgA|FlgM|FlgN") ~ "flgANM",
        str_detect(Protein,
                    "FlgB|FlgC|FlgD|FlgE|FlgF|FlgG|FlgH|FlgI|FlgJ|FlgK|FlgL") ~
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
        .default = protein_to_gene(Protein)),
    Role_II = case_when(
      str_detect(Protein, "RpoD|RpoN") ~ "Regulator: Sigma factor",
      str_detect(Protein, "FlhC|FlhD") ~ "Regulator: Master",
      str_detect(Protein,"FlrA/FleQ/FlaK|FlrC") ~ "Regulator: Alternative",
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
      str_detect(Protein, "FlgA|FlgB|FlgC|FlgF|FlgG") ~ "Basal-body/Hook: Rod"),
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
      factor(levels = c("Regulator", "fT3SS", "Switch",
                        "Basal-body/Hook", "Motor", "Filament")),
    Gene = protein_to_gene(Protein),
    Gene = case_match(Gene, "flrA/FleQ/FlaK" ~ "flrA/fleQ/flaK",
                 "fliO/FliZ" ~ "fliO/fliZ",
                 "fliR/FlhB" ~ "fliR/flhB",
                 .default = Gene)) |>
  arrange(Tier, Role_I, Role_II, Protein) |>
  mutate(Operon = as.factor(Operon) |> fct_inorder())

# prevalence and redundancy data

ko_list =
  read_tsv("data/proto_out/KO_search_results/KO-hit-counts.tsv") |>
  pivot_longer(cols = -c(assembly_id, total_gene_count),
               names_to = "domain_name", values_to = "n_copies") |>
  mutate(pa = if_else(n_copies > 0, 1, 0)) |>
  left_join(ko_meta) |>
  rename(Accession = assembly_id) |>
  left_join(metadata) |>
  filter(any(pa > 0), .by = Gene)


gene_prevalence = ko_list |>
  filter(any(pa > 0), .by = Accession) |>
  filter(any(pa > 0), .by = Gene) |>
  mutate(flagellated = if_else(str_detect(Gene, "fliC") & pa == 1, "yes", "no")) |>
  mutate(flagellated = any(flagellated == "yes"), .by = Accession) |>
  filter(flagellated, .by = Accession) |>
  summarise(pa = sum(pa)/n(),
            .by = c(Gene, Role_II, Role_I, Operon, Tier, definition)) |>
  arrange(-pa)


gene_prevalence |>
  write_csv("figures/fig_tab/SupplementaryData2.csv")

ggplot(gene_prevalence) + geom_histogram(aes(pa)) +
  scale_x_continuous(n.breaks = 20) +
  scale_y_continuous(n.breaks = 20)

rare = gene_prevalence |> filter(pa <= 0.5) |> pull(Gene)
accessory = gene_prevalence |> filter(pa > 0.5, pa <= 0.9) |> pull(Gene)
core = gene_prevalence |>  filter(pa > 0.9) |> pull(Gene)

simple_tree =
  ko_list |>
  mutate(Order_set = str_c(Order, n_distinct(Accession), association, sep = " - "),
         .by = c(Order, association)) |>
  select(Class, Order, Order_set) |>
  distinct() |>
  drop_na() |>
  mutate(across(everything(), as.factor)) |>
  ape::as.phylo.formula(~ Class/Order/Order_set, data = _)

tree_dat =
  ko_list |>
  mutate(
    Order_set = str_c(replace_na(Order, "Undefined"), n_distinct(Accession), association, sep = " - "),
    Accession = list(unique(Accession)),
    .by = c(Order, association)) |>
  summarise(pa = mean(pa), .by = c(Order_set, Accession, Operon, Tier))

collapsed_tree = rooted_tree


for (order in unique(tree_dat$Order_set)) {
  tips = tree_dat |>
    select(Order_set, Accession) |>
    distinct() |>
    filter(Order_set %in% order) |>
    pluck("Accession", 1)

  collapsed_tree = ape::drop.tip(collapsed_tree, tips[-1])
  collapsed_tree$tip.label[collapsed_tree$tip.label == tips[1]] = order

  }

plot(collapsed_tree)

nodes_class =
  metadata |>
  mutate(
    Order_set = str_c(replace_na(Order, "Undefined"), n_distinct(Accession), association, sep = " - "),
    .by = c(Order, association)) |>
  select(Class, Order_set) |>
  distinct() |>
  filter(str_detect(Order_set, "Undefined", negate = T)) |>
  group_by(Class) |>
  summarise(node = ggtree::MRCA(collapsed_tree, Order_set))


ggtree::ggtree(collapsed_tree, branch.length = T) +
  scale_color_manual(name = "set", values = palette_Compartment_Phylum,
                     aesthetics = c("fill", "color")) +
  ggtree::geom_highlight(data = nodes_class,
                         aes(node = node, fill = Class),
                         alpha = .5,
                         align = "right",
                         extend = 0.1,
                         type = "roundrect",
                         show.legend = F) +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom = geom_tile,
    data = tree_dat,
    aes(fill = pa,
        x = Operon,
        y = Order_set),
    offset = 0,
    pwidth = 1,
    linewidth = 0,
    axis.params = list(
      axis = "x",
      text.angle = 90,
      text.size = 10/.pt,
      line.size = 0,
      hjust = 1
    )) +
  colorspace::scale_fill_continuous_sequential(palette = "Purple-Blue", name = "Prevalence") +
  theme(plot.margin = margin(0, 0, 70, 0),
        legend.position = "bottom",
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(clip = "off") +
  ggtree::geom_tiplab(aes(label = str_remove(label, "- host_associated| - free_living| - endosymbionts| - pathogen"),
                          color = case_when(str_detect(label, "host") ~ "host_associated",
                                            str_detect(label, "free") ~ "free_living",
                                            str_detect(label, "endo") ~ "endosymbionts",
                                            str_detect(label, "path") ~ "pathogen")), align = T,
                      offset = 4.7,
                      linesize = 0,
                      show.legend = F) +
  xlim(NA, 15)

my_ggsave("fig2", 10, 10)

# mk  ---------------------------------------------------------------------

safe_fitDiscrete = safely(geiger::fitDiscrete)

# estimate transition rates

pa_fitDiscrete =
  ko_list |>
  mutate(present = if_else(pa == 0, "absent", "present") |>
           as.factor()) |>
  nest(.by = c(domain_name, definition, Protein, Operon, Role_II,
               Tier, Role_I, Gene)) |>
  mutate(
    data2 = map(data, ~ .x[match(rooted_tree$tip.label, .x$Accession),]),
    fit = map(data2,
              .progress = T,
              ~ safe_fitDiscrete(
                phy = ape::keep.tip(phy = rooted_tree,
                                    tip = pluck(.x, "Accession")),
                dat = pluck(.x, "present") |> set_names(pluck(.x, "Accession")),
                transform = "none",
                model = "ARD",
                ncores = 4)))

pa_fitDiscrete = readRDS("data/pa_fitDiscrete.RDS")

res_mk =
  pull(pa_fitDiscrete, "fit") |>
  set_names(pull(pa_fitDiscrete, "Gene")) |>
  discard_at(c("fliR/flhB", "flgQ"))

transition_matrix =
  res_mk |>
  map(~ pluck(.x, "result") |>
       geiger::as.Qmatrix.gfit() |>
       as.table() |> as.data.frame()) |>
  list_rbind(names_to = "Gene") |>
  rename(from = Var1, to = Var2, rate = Freq) |>
  mutate(transition_type = case_when(
    from == "present" ~ "loss",
    from == "absent" ~ "acquisition")) |>
  filter(from != to)

write_csv(transition_matrix, "figures/fig_tab/SupplementaryData3.csv")

transition_matrix |>
  summarise(median = median(rate),
            mad =  mad(rate),
            .by = c(transition_type)) |>
  mutate(across(where(is.numeric), ~ round(.x, 4))) |>
  arrange(median)

transition_matrix |>
  ggplot(aes(x = rate, fill = transition_type)) +
  geom_histogram(alpha = .5, bins = 55) +
  scale_fill_manual(values = c("loss" = "tan", "acquisition" =  "slategrey")) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size = 7)))

my_ggsave("suppFig.1", 3.5, h = 3.5)


# prevalence -------------------------------------------------------------------

safe_phyloglm = safely(phylolm::phyloglm)

# iterate over each ortholog and fit model

res_fit =
  ko_list |>
  mutate(lGS = log_genome_size - mean(log_genome_size), .by = Gene) |>
  nest(.by = c(Gene, Protein, domain_name, definition, Tier, Role_I, Role_II, Operon)) |>
  mutate(data2 = map(data, ~ column_to_rownames(.x, "Accession")),
         fit_logit = map(data2,
                         ~ safe_phyloglm(
                           pa ~  set * lGS,
                           data = .x,
                           btol = 100,
                           log.alpha.bound = 10,
                           phy = rooted_tree |> ape::keep.tip(phy = _, tip = rownames(.x)),
                           method = "logistic_MPLE")),
         fit_pois = map(data2,
                        ~ safe_phyloglm(
                          n_copies ~ set * lGS,
                          data = .x,
                          phy = rooted_tree |> ape::keep.tip(phy = _, tip = rownames(.x)),
                          method = "poisson_GEE")))

res_fit = read_rds("data/logistic_poisson.RDS")

# tidy up results

logit_res =
  res_fit |>
  mutate(res_logit = map(fit_logit,
                         ~ pluck(.x, "result") |>
                           summary() |>
                           pluck("coefficients") |>
                           as_tibble(rownames = "term"))) |>
  select(-c(data, data2, fit_logit, fit_pois)) |>
  unnest(res_logit) |>
  mutate(
    term = str_remove(term, "set"),
    estimate = Estimate,
    std.error = StdErr,
    statistic = z.value,
    fdr(p.value, 0.05),
    .keep = "unused") |>
  filter(term != "(Intercept)")

logit_res |> count(term, significant)

filter(logit_res, term %in% "host_associated") |>
  arrange(Tier, Role_I, Role_II) |>
  select(Gene, estimate, std.error, statistic, FDR) |>
  mutate(across(c(estimate, std.error, statistic), ~ round(.x, 2)),
         FDR = format(FDR, scientific = T, digits = 2)) |>
  write_csv("figures/fig_tab/SupplementaryData4.csv")


summarise(ko_list,
          pa = sum(pa)/n(),
          .by = c(Gene, set, Tier, Role_I)) |>
  ggplot(aes(x = pa, y = Gene)) +
  geom_path(
    aes(group = Gene),
    linewidth = 3,
    color = "#cccccc") +
  geom_point(
    aes(fill = set),
    shape = 21,
    size = 3,
    stroke = 1,
    color = "#555555") +
  geom_text(
    data = logit_res |>
      filter(term %in% "host_associated"),
    aes(label = format(FDR, scientific = T, digits = 2),
        fontface = if_else(significant, "bold", "plain"),
        y = Gene,
        x = 1.1),
    color = "#555555",
    hjust = 0.35) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0, hjust = 0,
                                    face = "bold"),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.y = element_line(color = "#555555",
                                          linewidth = 0.1),
        panel.grid.major.x = element_line(color = "#555555",
                                          linewidth = 0.1),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_manual(values = pal_set, labels = list(
    host_associated = "Host Associated",
    free_living = "Free Living")) +
  facet_grid(Tier + Role_I ~ ., space = "free_y",
             scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = 0.1),
                     name = "Prevalence (%)",
                     breaks = seq(0, 1, 0.2)) +
  guides(fill = guide_legend(override.aes = list(size = 10)))

my_ggsave("fig3", w = 9, h = 10)

pois_res =
  res_fit |>
  mutate(res_pois = map(fit_pois,
                        ~ pluck(.x, "result") |>
                          summary() |>
                          pluck("coefficients") |>
                          as_tibble(rownames = "term"))) |>
  select(-c(data, data2, fit_logit, fit_pois)) |>
  unnest(res_pois) |>
  mutate(
    term = str_remove(term, "set"),
    estimate = Estimate,
    std.error = StdErr,
    statistic = z.value,
    fdr(p.value, 0.05),
    .keep = "unused") |>
  filter(term != "(Intercept)")

pois_res |>
  count(term, significant)

filter(pois_res, term %in% "host_associated") |>
  arrange(Tier, Role_I, Role_II) |>
  select(Gene, estimate, std.error, statistic, FDR) |>
  mutate(across(c(estimate, std.error, statistic), ~ round(.x, 2)),
         FDR = format(FDR, scientific = T, digits = 2)) |>
  write_csv("figures/fig_tab/SupplementaryData5.csv")

summarise(ko_list,
          multiple_copies = sum(n_copies > 1)/n(),
          .by = c(Gene, set)) |>
  left_join(pois_res) |>
  filter(multiple_copies > 0, term %in% "host_associated") |>
  ggplot(aes(x = multiple_copies, y = Gene, fill = set)) +
  geom_col(
    color = "white",
    position = position_dodge()) +
  geom_text(
    data = ~ filter(.x, set %in% "host_associated"),
    aes(x = .6, y = Gene, label = format(FDR, scientific = T, digits = 2))) +
  scale_fill_manual(values = pal_set,
                    labels = list("free_living" = "Free Living",
                                  "host_associated" = "Host Associated")) +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  scale_x_continuous(expand = expansion(mult = 0.2),
                     name = "Share of Genomes carrying Multiple Copies",
                     breaks = seq(0, 6, 0.2))

my_ggsave("suppFig2", 6, 8)

filter(pois_res, term == "host_associated", significant) |>
  View()

# bitscore ---------------------------------------------------------

hmmsearch_results =
  read_csv(
    "data/proto_out/KO_search_results/hmmsearch_extended/hmmsearch_results.csv",
    col_names = c("query_name", "domain_name",
                  "sequence_evalue", "sequence_score", "sequence_bias",
                  "domain_evalue", "domain_score", "domain_bias")) |>
  mutate(Accession = str_remove(query_name, "_[0-9]*$"))

ko_tblout =
  hmmsearch_results |>
  left_join(metadata) |>
  filter(domain_score > 0,
         sequence_score > 0,
         sequence_score > 2 * sequence_bias,
         domain_score > 2 * domain_bias,
         sequence_score < 2 * domain_score) |>
  mutate(n_hits = n(), .by = c(Accession, domain_name)) |>
  slice_min(sequence_evalue, n = 1, by = c(Accession, domain_name), with_ties = F) |>
  mutate(
    dbs =  (sequence_score - mean(sequence_score))/sd(sequence_score),
    .by = domain_name) |>
  left_join(ko_meta) |>
  mutate(lGS = log_genome_size - mean(log_genome_size), .by = Gene)

ko_tblout |>
  summarise(lGS = first(log_genome_size),
            dbs = mean(dbs),
            p_of_gene = n()/53,
            n_hits = mean(n_hits) |>
              log(), .by = c(Accession, association)) |>
  pivot_longer(c(n_hits, dbs)) |>
  ggplot(aes(lGS, value, color = association, size = p_of_gene)) +
  geom_point(alpha = 2/3, ) +
  scale_size(range = c(0.1, 2)) +
  # geom_smooth(alpha = 1/3, fill = "#d3d3d3", method = "lm") +
  scale_color_manual(values = pal_set) +
  facet_wrap(~name, scales = "free_x")

# # bitscore model ----------------------------------------------------------
#

# n. of hits is correlated with the genome size, so we can exclude it

safe_gls = safely(.f = nlme::gls)

brownian_gls =
  ko_tblout |>
  nest(.by = c(Gene, Protein, domain_name, definition, Tier, Role_I, Role_II, Operon)) |>
  mutate(
    fit = map(data,
              ~ safe_gls(
                dbs ~  set * lGS,
                data = .x,
                correlation = ape::corBrownian(
                  value = 1,
                  phy = ape::keep.tip(phy = rooted_tree,
                                      tip = .x$Accession),
                  form = ~ Accession),
                verbose = T,
                control =
                  nlme::glsControl(maxIter = 1e3,
                                   msMaxIter = 1e3,
                                   tolerance = 1e-5,
                                   singular.ok = T))))

# results -----------------------------------------------------------------

mod_tidy =
  # choose_model |>
  # filter(model == "brownian") |>
  brownian_gls |>
  hoist(fit, "result") |>
  mutate(tidied = map(result,
                      tidy)) |>
  unnest(tidied) |>
  filter(term != "(Intercept)") |>
  mutate(fdr(p.value, 0.05),
         term = str_remove(term, "set"))

filter(mod_tidy, term %in% "host_associated") |>
  arrange(Tier, Role_I, Role_II) |>
  select(Gene, estimate, std.error, statistic, FDR) |>
  mutate(across(c(estimate, std.error, statistic), ~ round(.x, 2)),
         FDR = format(FDR, scientific = T, digits = 2)) |>
  write_csv("figures/fig_tab/SupplementaryData6.csv")

  # mutate(term =
  #          case_match(
  #            term,
  #            "sethost_associated" ~ "host_associated",
  #            "setfree_living:lGS" ~ "Genome Size (free_living)",
  #            "setfree_living:log(n_hits)" ~ "N. of hits (free_living)",
  #            "sethost_associated:lGS" ~ "Genome Size (host_associated)",
  #            "sethost_associated:log(n_hits)" ~ "N. of hits (host_associated)") |>
  #          factor(levels = c("host_associated",
  #                            "Genome Size (host_associated)",
  #                            "N. of hits (host_associated)",
  #                            "Genome Size (free_living)",
  #                            "N. of hits (free_living)")),
  #        fdr(p.value, 0.05))

# compare results from logistic and continuous models

mod_tidy |>
  count(term, significant)

mod_tidy |>
  summarise(mean(estimate), sd(estimate), .by = term)

left_join(mod_tidy |>
            filter(term %in% "host_associated"),
          logit_res |>
            filter(term %in% "host_associated"),
          by = names(ko_meta)) |>
  mutate(sign = case_when(
    significant.x & significant.y ~ "both",
    significant.x ~ "dbs",
    significant.y ~ "prev",
    .default = "none"
  )) |>
  ggplot() +
  geom_text(aes(estimate.x, estimate.y, label = Gene, color = sign)) +
  theme(legend.position = "bottom")

p3 =
  mod_tidy |>
  mutate(text_xpos = max(estimate + std.error)) |>
  ggplot() +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           color = NA, fill = "#555555", alpha = .1) +
  geom_pointrange(
    aes(xmin = estimate - std.error,
        x = estimate,
        xmax = estimate + std.error,
        y = Gene,
        color = term,
        shape = term),
    position = position_dodge(width = .5),
    size = 0.5,
    stroke = 0.5,
    linewidth = 0.5) +
  geom_text(
    data = ~ filter(.x, term %in% "host_associated"),
    aes(label = format(FDR, scientific = T, digits = 2),
        y = Gene,
        x = text_xpos),
    color = "#555555",
    hjust = 0.35) +
  facet_grid(Tier + Role_I ~ .,
             scales = "free",
             space = "free_y")  +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0,
                                    face = "bold"),
        strip.placement = "outside",
        legend.position = "bottom",
        axis.line.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.ticks.y = element_blank()) +
  scale_color_manual(values = c("#CC79A7", "grey50", "grey75"),
                     labels = list(
                       host_associated = "Host Associated",
                       `host_associated:lGS` = "Host Associated & Genome Size",
                       lGS = "Genome Size")) +
  scale_shape_manual(values = c(19, 15, 17),
                     labels = list(
                       host_associated = "Host Associated",
                       `host_associated:lGS` = "Host Associated & Genome Size",
                       lGS = "Genome Size")) +
  guides(alpha = "none", color = guide_legend(override.aes = list(size = 2))) +
  coord_cartesian(clip = "off") +
  labs(x = "Estimated Effect on Sequence Bitscore")

p3

p4 =
  mod_tidy |>
  ggplot() +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           color = NA, fill = "#555555", alpha = .1, ) +
  geom_density(
    aes(x = estimate, fill = term),
    color = "white",
    alpha = .7,
    outline.type = "upper") +
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
        legend.position = "none")

p4

p4 + p3 +
  plot_layout(heights = c(1, 9),
              axes = "collect_x")

my_ggsave("fig_tab/fig4", 9, 11)


list_rbind(
  list("prevalence" = logit_res,
       "conservation" = mod_tidy),
  names_to = "effect") |>
  filter(term %in% "host_associated", significant) |>
  count(Gene, effect) |>
  pivot_wider(values_from = n, names_from = effect) |> View()

