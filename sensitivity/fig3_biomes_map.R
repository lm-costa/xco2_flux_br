suppressMessages({library(dplyr); library(ggplot2); library(sf)}); sf::sf_use_s2(FALSE)
setwd("D:/R/xco2_flux_br")
b <- geobr::read_biomes(year = 2019) |> filter(name_biome != "Sistema Costeiro")
sa <- sf::read_sf(list.files("data-raw/South_America/", pattern = "shp", full.names = TRUE)[1])
lv <- c("Amazônia" = "Amazon (AMZ)", "Caatinga" = "Caatinga (CAAT)", "Cerrado" = "Cerrado (CER)", "Mata Atlântica" = "Atlantic Forest (AF)", "Pampa" = "Pampa (PMP)", "Pantanal" = "Pantanal (PNT)")
b$bio <- factor(lv[as.character(b$name_biome)], levels = unname(lv))
cols <- c("Amazon (AMZ)" = "#5fa35b", "Caatinga (CAAT)" = "#f4a65c", "Cerrado (CER)" = "#f6e36b", "Atlantic Forest (AF)" = "#6f9ccb", "Pampa (PMP)" = "#e48ad2", "Pantanal (PNT)" = "#8fdcf0")
has_sp <- requireNamespace("ggspatial", quietly = TRUE); cat("ggspatial:", has_sp, "\n")
p <- ggplot() + geom_sf(data = sa, fill = "grey93", colour = "grey70", linewidth = 0.25) +
  geom_sf(data = b, aes(fill = bio), colour = "grey25", linewidth = 0.35, alpha = 0.95) +
  scale_fill_manual(values = cols, name = "Biomes") +
  coord_sf(xlim = c(-75, -33), ylim = c(-35.5, 6), expand = FALSE) +
  scale_x_continuous(breaks = seq(-70, -35, 5)) + scale_y_continuous(breaks = seq(-35, 5, 5)) +
  theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = "#a9d3e6", colour = NA), panel.grid = element_line(colour = "white", linewidth = 0.2), panel.ontop = FALSE,
        legend.position = c(0.2, 0.2), legend.background = element_rect(fill = "white", colour = "grey60"), legend.title = element_text(face = "bold"), axis.title = element_blank())
if (has_sp) p <- p + ggspatial::annotation_scale(location = "bl", width_hint = 0.25, pad_x = unit(0.5, "cm"), pad_y = unit(0.4, "cm")) + ggspatial::annotation_north_arrow(location = "tr", which_north = "true", style = ggspatial::north_arrow_minimal(), height = unit(1.2, "cm"), width = unit(1.2, "cm"))
ggsave("figs_manuscript/fig3_biomes.png", p, width = 7.3, height = 7.4, dpi = 200, bg = "white"); cat("ok\n")
