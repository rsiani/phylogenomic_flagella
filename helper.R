### Helper functions and other stuff
## Roberto Siani
# 2024

# load libraries

pacman::p_load(pacman, tidyverse, patchwork,
              broom, broom.mixed, nlme, ggraph, tidygraph)

# ggplot theme for visualization

theme_set(
  theme_minimal() +
    theme(
      plot.margin = margin(3, 3, 3, 3),
      text = element_text(size = 15,
                          family = "Fira Sans",
                          color = "#555555"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "#555555", linewidth = .5),
      axis.ticks = element_line(color = "#555555", linewidth = .5),
      axis.title = element_text(hjust = 1),
      legend.position = "none",
      legend.title = element_blank(),
      panel.spacing.x = unit(.5, "lines"),
      panel.spacing.y = unit(.5, "lines"),
    )
)

# read HMMER tblout, modified from https://github.com/arendsee/rhmmer

read_tblout <- function(file){

  col_types <-
    readr::cols(
      query_name         = readr::col_character(),
      query_accession    = readr::col_character(),
      domain_name          = readr::col_character(),
      domain_accession     = readr::col_character(),
      sequence_evalue     = readr::col_double(),
      sequence_score      = readr::col_double(),
      sequence_bias       = readr::col_double(),
      best_domain_evalue  = readr::col_double(),
      best_domain_score   = readr::col_double(),
      best_domain_bis     = readr::col_double(),
      domain_number_exp   = readr::col_double(),
      domain_number_reg   = readr::col_integer(),
      domain_number_clu   = readr::col_integer(),
      domain_number_ov    = readr::col_integer(),
      domain_number_env   = readr::col_integer(),
      domain_number_dom   = readr::col_integer(),
      domain_number_rep   = readr::col_integer(),
      domain_number_inc   = readr::col_character()
    )

  N <- length(col_types$cols)

  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy=FALSE, progress=FALSE)

  table <- sub(
    pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
    replacement = '\\1',
    x=lines,
    perl = TRUE
  ) %>%
    gsub(pattern="  *", replacement="\t") %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(
      col_names=names(col_types$cols),
      comment='#',
      na='-',
      col_types = col_types,
      lazy=FALSE,
      progress=FALSE
    )

  descriptions <- lines[!grepl("^#", lines, perl=TRUE)] %>%
    sub(
      pattern = sprintf("%s *(.*)", paste0(rep('\\S+', N), collapse=" +")),
      replacement = '\\1',
      perl = TRUE
    )

  table$description <- descriptions[!grepl(" *#", descriptions, perl=TRUE)]

  table
}

# extract gene ID

protein_id = \(.x) str_extract(.x,
                                 "\\p{alpha}\\p{lower}{2}\\d?\\p{upper}\\b")

protein_to_gene =
  function(.x) {
    str_sub(.x, 1, 1) <- str_to_lower(str_sub(.x, 1, 1)); .x
    }

# palette for phylogenetic comparison

pal_set = c("Other" = "#DDDDDD", "free_living" = "#44AA99", "host_associated" = "#CC6677",
            "endosymbionts" = "#882255")

palette_Compartment_Phylum =
  c("#DDDDDD", "#44AA99", "#CC6677", "#AA4499","#DDCC77",
    "grey30", "grey90", "grey60") |>
  set_names(c("Other", "free_living", "host_associated",
              "endosymbionts",
              "pathogen",
              "Alphaproteobacteria",
              "Betaproteobacteria",
              "Gammaproteobacteria"))

# preconfigured figure saving

my_ggsave =
  function(.x, w, h) {
    ggsave(filename = str_c("figures/", .x, ".svg", sep = ""),
           width = w,
           height = h,
           device = svglite::svglite,
           units = "in",
    )
  }

# simple fdr adjustment

fdr = \(.x, fdr_level) tibble(
  p.value = .x,
  FDR = p.adjust(.x, "fdr"),
  significant = FDR <= fdr_level,
  s.value = -log2(.x))
#
# fisher.Z =
#   function(.x, .y, n.x, n.y) {
#     Z.x = 0.5 * log((1 + .x)/(1 - .x))
#     Z.y = 0.5 * log((1 + .y)/(1 - .y))
#     res = data.frame(statistic = (Z.x - Z.y) / sqrt((1 / (n.x - 3)) + (1 / (n.y - 3)))) |>
#       mutate(p.value = 2 * (1 - pnorm(abs(statistic))))
#     return(res)
#     }

okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",  "#F0E442")

lifestyle_change =
  function(.state1, .state2) {
    .state1 = str_remove(.state1, ".[0-1]")
    .state2 = str_remove(.state2, ".[0-1]")
    .transition = if_else(.state1 == .state2, .state1, paste0(.state1, " -> ", .state2))
    return(.transition)
    }
