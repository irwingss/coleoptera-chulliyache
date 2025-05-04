# ──────────────────────────────
# 1. Load required libraries
# ──────────────────────────────
library(tidyverse)
library(vegan)
library(factoextra)
library(ape)
library(rstatix)
library(ggstatsplot)
library(ggthemes)

# ──────────────────────────────
# 2. Load and split the dataset
# ──────────────────────────────
# The file contains ecological sampling data on coleopteran species found on vertebrate carcasses
data <- openxlsx::read.xlsx("data/Datos_Chu.xlsx")

# Metadata includes spatial info (Zone), vertebrate species, decomposition stage, family, and coordinates
metadata <- data[, 1:6]

# Species abundance data (coleoptera counts) starts from column 7 onward
species_data <- data[, 7:15]

# Standardize column names for consistency
colnames(species_data) <- c("Dermestes ater", "Dermestes frischii", "Dermestes maculatus", 
                            "Necrobia rufipes", "Psammetichus sp.", "Aleochara sp.",
                            "Paraneda pallidula guticollis", "Phaleria psammatea", "Hypocaccus gaudens")

# Check total abundance per species
colSums(species_data)

# ──────────────────────────────
# 3. Principal Component Analysis (PCA)
# ──────────────────────────────

# Apply Hellinger transformation to abundance data (suitable for community composition analysis)
species_hellinger <- decostand(species_data, method = "hellinger")

# PCA on standardized data
pca_model <- prcomp(species_hellinger, center = TRUE, scale. = TRUE)

# Summary of explained variance
summary(pca_model)

# Optional: Rename rows for species
rownames(pca_model$rotation) <- colnames(species_data)

# ──────────────────────────────
# 4. Custom PCA plots using user-defined function
# ──────────────────────────────

# PCA by Zone
PCA_Zone <- Plot_pca(pca_model, vector_grupos = metadata$Zona,
                     scale_factor = 8, use_ellipse = FALSE,
                     ellipse_alpha = 0.2, ellipse_lwd = 0,
                     name_groups = "Zone", show_group_shape = TRUE, point_alpha = 1) +
  scale_fill_tableau() + scale_color_tableau()

# PCA by Decomposition Stage
PCA_Decomposition <- Plot_pca(pca_model, vector_grupos = metadata$Condición,
                              scale_factor = 8, use_ellipse = FALSE,
                              ellipse_alpha = 0.2, ellipse_lwd = 0,
                              name_groups = "Decomposition Stage", show_group_shape = TRUE, point_alpha = 1) +
  scale_fill_tableau() + scale_color_tableau()

# PCA by Vertebrate Species
PCA_Species <- Plot_pca(pca_model, vector_grupos = metadata$Especie.de.vertebrado,
                        scale_factor = 8, use_ellipse = FALSE,
                        ellipse_alpha = 0.2, ellipse_lwd = 0,
                        name_groups = "Vertebrate Species", show_group_shape = FALSE, point_alpha = 1) +
  scale_fill_tableau() + scale_color_tableau() +
  theme(legend.text = element_text(face = 3))

# Save PCA plots
ggsave("figures/PCA_Zone.png", plot = PCA_Zone, width = 20, height = 20, units = "cm", dpi = 600, bg = "white")
ggsave("figures/PCA_Decomposition.png", plot = PCA_Decomposition, width = 20, height = 20, units = "cm", dpi = 600, bg = "white")
ggsave("figures/PCA_Species.png", plot = PCA_Species, width = 20, height = 20, units = "cm", dpi = 600, bg = "white")

# ──────────────────────────────
# 5. Diversity Indices: Shannon and Simpson
# ──────────────────────────────

diversity(data %>% select(., 7:15), index = "shannon")

diversity_data <- data %>%
  mutate(
    N = rowSums(select(., 7:15)),
    Shannon = diversity(select(., 7:15), index = "shannon"),
    Simpson = diversity(select(., 7:15), index = "simpson"),
    Zona = as.factor(Zona),
    Condición = as.factor(Condición),
    Especie.de.vertebrado = as.factor(Especie.de.vertebrado)
  )

# Ranges of diversity (mean and sd)
range(diversity_data$Shannon)
mean(diversity_data$Shannon)
sd(diversity_data$Shannon)

range(diversity_data$Simpson)
mean(diversity_data$Simpson)
sd(diversity_data$Simpson)

# ──────────────────────────────
# 6. ANOVA by Zone
# ──────────────────────────────

diversity_data_sha %>% group_by(Zona) %>% 
  summarise(Shanon_med = median(Shannon),
            Shannon_q1 = 
              quantile(Shannon, probs = c(0.25)),
            Shannon_q2 = 
              quantile(Shannon, probs = c(0.75)),
            Simpson_med = median(Simpson),
            Simpson_q1 = 
              quantile(Simpson, probs = c(0.25)),
            Simpson_q2 = 
              quantile(Simpson, probs = c(0.75))
  )

# Shannon 
diversity_data_sha <- diversity_data %>% filter(Shannon != 0)

levene_test(Shannon ~ Zona, data = diversity_data_sha) #Homokedasticity
diversity_data_sha %>% group_by(Zona) %>%  shapiro_test(Shannon) #Non-normality
kruskal_test(Shannon ~ Zona, data = diversity_data_sha)
dunn_test(Shannon ~ Zona, data = diversity_data_sha, p.adjust.method = "holm")

gg_box_zone_shannon <- ggstatsplot::ggbetweenstats(
  data = diversity_data_sha, x = Zona, y = Shannon,
  type = "np", plot.type = "box", p.adjust.method = "holm",
  bf.message = FALSE, pairwise.comparisons = FALSE,
  centrality.point.args = list(color="black", size=4, alpha=0.5),
  boxplot.args = list(aes(color=Zona), fill="transparent"),
  violin.args = list(color="transparent", fill="transparent"),
  centrality.plotting = TRUE, results.subtitle = TRUE
) +
  ggthemes::scale_color_tableau(direction = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_blank(), 
    legend.position = "none",
    text = element_text(size = 13)) +
  labs(x = "Zone", 
       y = bquote("Shannon diversity index (" ~ H~ ")"))
gg_box_zone_shannon

ggsave("figures/NEW_Boxplot_shannon_Zone.png", gg_box_zone_shannon, width = 20, height = 16, units = "cm", dpi = 600)


# Simpson
diversity_data_simp <- diversity_data %>% filter(Simpson != 0)

levene_test(Simpson ~ Zona, data = diversity_data_simp) #Homokedasticity
diversity_data_simp  %>% group_by(Zona) %>%  shapiro_test(Simpson) #Non-normality
kruskal_test(Simpson ~ Zona, data = diversity_data_simp )
dunn_test(Simpson ~ Zona, data = diversity_data_simp )

gg_box_zone_simpson <- ggstatsplot::ggbetweenstats(
  data = diversity_data_simp , x = Zona, y = Simpson,
  type = "np", plot.type = "box", p.adjust.method = "holm",
  bf.message = FALSE, pairwise.comparisons = FALSE,
  centrality.point.args = list(color="black", size=4, alpha=0.5),
  boxplot.args = list(aes(color=Zona), fill="transparent"),
  violin.args = list(color="transparent", fill="transparent"),
  centrality.plotting = TRUE, results.subtitle = TRUE
) +
  ggthemes::scale_color_tableau(direction = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_blank(), 
    legend.position = "none",
    text = element_text(size = 13)) +
  labs(x = "Zone", 
       y = bquote("Simpson diversity index" * ~ (1 - D)))
gg_box_zone_simpson

ggsave("figures/NEW_Boxplot_Simpson_Zone.png", gg_box_zone_simpson, width = 20, height = 16, units = "cm", dpi = 600)

# Unified plot

library(ggpubr)
ZONE <- ggarrange(gg_box_zone_shannon, gg_box_zone_simpson, labels = c("A","B"))
ggsave("figures/NEW_Boxplot_ZONE.png", ZONE, 
       width = 40, height = 16, units = "cm", dpi = 600)


# ──────────────────────────────
# 7. ANOVA by Decomposition Stage
# ──────────────────────────────

diversity_data_sha %>% group_by(Condición) %>% 
  summarise(Shanon_med = median(Shannon),
            Shannon_q1 = 
              quantile(Shannon, probs = c(0.25)),
            Shannon_q2 = 
              quantile(Shannon, probs = c(0.75)),
            Simpson_med = median(Simpson),
            Simpson_q1 = 
              quantile(Simpson, probs = c(0.25)),
            Simpson_q2 = 
              quantile(Simpson, probs = c(0.75))
  )


# Shannon 
diversity_data_sha <- diversity_data %>% filter(Shannon != 0)

levene_test(Shannon ~ Condición, data = diversity_data_sha) #Homokedasticity
diversity_data_sha %>% group_by(Condición) %>%  shapiro_test(Shannon) #Non-normality
kruskal_test(Shannon ~ Condición, data = diversity_data_sha)
dunn_test(Shannon ~ Condición, data = diversity_data_sha)

gg_box_compo_shannon <- ggstatsplot::ggbetweenstats(
  data = diversity_data_sha, x = Condición, y = Shannon,
  type = "np", plot.type = "box", p.adjust.method = "holm",
  bf.message = FALSE, pairwise.comparisons = FALSE,
  centrality.point.args = list(color="black", size=4, alpha=0.5),
  boxplot.args = list(aes(color=Condición), fill="transparent"),
  violin.args = list(color="transparent", fill="transparent"),
  centrality.plotting = TRUE, results.subtitle = TRUE
) +
  ggthemes::scale_color_tableau(direction = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_blank(), 
    legend.position = "none",
    text = element_text(size = 13)) +
  labs(x = "Decay Stage", 
       y = bquote("Shannon diversity index (" ~ H~ ")"))
gg_box_compo_shannon

ggsave("figures/NEW_Boxplot_shannon_Compo.png", gg_box_zone_shannon, width = 20, height = 16, units = "cm", dpi = 600)


# Simpson
diversity_data_simp <- diversity_data %>% filter(Simpson != 0)

levene_test(Simpson ~ Condición, data = diversity_data_simp) #Homokedasticity
diversity_data_simp  %>% group_by(Condición) %>%  shapiro_test(Simpson) #Non-normality
kruskal_test(Simpson ~ Condición, data = diversity_data_simp )
dunn_test(Simpson ~ Condición, data = diversity_data_simp )

gg_box_compo_simpson <- ggstatsplot::ggbetweenstats(
  data = diversity_data_simp , x = Condición, y = Simpson,
  type = "np", plot.type = "box", p.adjust.method = "holm",
  bf.message = FALSE, pairwise.comparisons = FALSE,
  centrality.point.args = list(color="black", size=4, alpha=0.5),
  boxplot.args = list(aes(color=Condición), fill="transparent"),
  violin.args = list(color="transparent", fill="transparent"),
  centrality.plotting = TRUE, results.subtitle = TRUE
) +
  ggthemes::scale_color_tableau(direction = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_blank(), 
    legend.position = "none",
    text = element_text(size = 13)) +
  labs(x = "Decay Stage", 
       y = bquote("Simpson diversity index" * ~ (1 - D)))
gg_box_compo_simpson

ggsave("figures/NEW_Boxplot_Simpson_Compo.png", gg_box_compo_simpson, width = 20, height = 16, units = "cm", dpi = 600)

# Unified plot

library(ggpubr)
COMPO <- ggarrange(gg_box_compo_shannon, gg_box_compo_simpson, labels = c("A","B"))
ggsave("figures/NEW_Boxplot_COMPO.png", COMPO, 
       width = 40, height = 16, units = "cm", dpi = 600)

# ──────────────────────────────
# 8. ANOVA by Vertebrate Species
# ──────────────────────────────

diversity_data_sha %>% group_by(Especie.de.vertebrado) %>% 
  summarise(Shanon_med = median(Shannon),
            Shannon_q1 = 
              quantile(Shannon, probs = c(0.25)),
            Shannon_q2 = 
              quantile(Shannon, probs = c(0.75)),
            Simpson_med = median(Simpson),
            Simpson_q1 = 
              quantile(Simpson, probs = c(0.25)),
            Simpson_q2 = 
              quantile(Simpson, probs = c(0.75))
  ) %>% as.data.frame()

# Removing species with just one record

freq_spe <- diversity_data %>% count(Especie.de.vertebrado) %>% arrange(n)
del_spe <- freq_spe[1:4,"Especie.de.vertebrado"]

filtered_data <- diversity_data %>%
  filter(!Especie.de.vertebrado %in% del_spe)

# Shannon 
diversity_data_sha <- filtered_data %>% filter(Shannon != 0)

levene_test(Shannon ~ Especie.de.vertebrado, data = diversity_data_sha) #Homokedasticity
diversity_data_sha %>% group_by(Especie.de.vertebrado) %>%  shapiro_test(Shannon) #Non-normality
kruskal_test(Shannon ~ Especie.de.vertebrado, data = diversity_data_sha)
dunn_test(Shannon ~ Especie.de.vertebrado, data = diversity_data_sha)

gg_box_EspVert_shannon <- ggstatsplot::ggbetweenstats(
  data = diversity_data_sha, x = Especie.de.vertebrado, y = Shannon,
  type = "np", plot.type = "box", p.adjust.method = "holm",
  bf.message = FALSE, pairwise.comparisons = FALSE,
  centrality.point.args = list(color="black", size=4, alpha=0.5),
  boxplot.args = list(aes(color=Especie.de.vertebrado), fill="transparent"),
  violin.args = list(color="transparent", fill="transparent"),
  centrality.plotting = TRUE, results.subtitle = TRUE
) +
  ggthemes::scale_color_tableau(direction = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_blank(), 
        legend.position = "none",
        text = element_text(size = 13)) +
  labs(x = "Decay Stage", 
       y = bquote("Shannon diversity index (" ~ H~ ")"))
gg_box_EspVert_shannon

ggsave("figures/NEW_Boxplot_shannon_EspVert.png", gg_box_EspVert_shannon, width = 20, height = 16, units = "cm", dpi = 600)


# Simpson
diversity_data_simp <- filtered_data %>% filter(Simpson != 0)

levene_test(Simpson ~ Especie.de.vertebrado, data = diversity_data_simp) #Homokedasticity
diversity_data_simp  %>% group_by(Especie.de.vertebrado) %>%  shapiro_test(Simpson) #Non-normality
kruskal_test(Simpson ~ Especie.de.vertebrado, data = diversity_data_simp )
dunn_test(Simpson ~ Especie.de.vertebrado, data = diversity_data_simp )

gg_box_EspVert_simpson <- ggstatsplot::ggbetweenstats(
  data = diversity_data_simp , x = Especie.de.vertebrado, y = Simpson,
  type = "np", plot.type = "box", p.adjust.method = "holm",
  bf.message = FALSE, pairwise.comparisons = FALSE,
  centrality.point.args = list(color="black", size=4, alpha=0.5),
  boxplot.args = list(aes(color=Especie.de.vertebrado), fill="transparent"),
  violin.args = list(color="transparent", fill="transparent"),
  centrality.plotting = TRUE, results.subtitle = TRUE
) +
  ggthemes::scale_color_tableau(direction = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_blank(), 
        legend.position = "none",
        text = element_text(size = 13)) +
  labs(x = "Decay Stage", 
       y = bquote("Simpson diversity index" * ~ (1 - D)))
gg_box_EspVert_simpson

ggsave("figures/NEW_Boxplot_Simpson_EspVert.png", gg_box_EspVert_simpson, width = 20, height = 16, units = "cm", dpi = 600)

# Unified plot

library(ggpubr)
ESPVERT <- ggarrange(gg_box_EspVert_shannon, gg_box_EspVert_simpson, labels = c("A","B"))
ggsave("figures/NEW_Boxplot_ESPVERT.png", ESPVERT, 
       width = 40, height = 16, units = "cm", dpi = 600)
