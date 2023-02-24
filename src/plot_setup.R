# Aesthetic options used in plots ----
# standardize visuals across scripts

## basic elements ----
theme_set(theme_bw())
alpha_plt <- 0.4
alpha_plt_sensitivity <- 0.6

## legend sizing ----
### main analyses ----
theme_single_panel <- theme(axis.title = element_text(size = 8),
                            axis.text = element_text(size = 7),
                            legend.title = element_text(size = 7),
                            legend.text = element_text(size = 7),
                            legend.key.size = unit(.5, "cm"),
                            legend.position = "right")

# put shape before colour and make legend key bigger
guides_man_fig1.2 <- guides(
  shape = guide_legend(order = 1, override.aes = list(size = 4)),
  colour = guide_legend(order = 2, override.aes = list(size = 1.5)),
  fill = guide_legend(order = 2)
)

theme_man_fig1.2 <- theme(
  # make axis size legible
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 24, hjust = 1),
  
  # put legend in empty space, make legend key bigger and text legible
  legend.position = c(.92, .4),
  # legend.box.just = "left",
  # legend.key.height = unit(1, "lines"),
  legend.key.width = unit(2, "lines"),
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 18),
  
  # remove strips and add space between panel rows and columns
  strip.text = element_blank(),
  panel.spacing.x = unit(1, "lines"),
  panel.spacing.y = unit(2, "lines")
)

# settings for plots of rate distributions (supplementary)
theme_pointcloud <- theme(
  # resize title tag
  plot.title = element_text(hjust = -0.1, size = 8),
  
  # make axis and facet text legible
  axis.title = element_text(size = 8),
  axis.text.y = element_text(size = 6.5),
  axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1),
  strip.text = element_text(size = 6.5, margin = margin(0.08, 0, 0.08, 0, "cm")),
  
  # make legend key bigger and text legible
  legend.key.width = unit(0.7, "line"),
  legend.text = element_text(size = 6.5),
  legend.position = "none",
  
  plot.background = element_blank()
)

# put shape before colour and make legend key bigger
guides_pointcloud <- guides(
  shape = guide_legend(order = 1, override.aes = list(size = 2)),
  colour = guide_legend(order = 2, override.aes = list(size = 1))
)

### colour schemes ----
# counterfactual/observed colour scheme
# met.brewer("Hiroshige", 10, "discrete")
# met.brewer("Hiroshige", 10, "discrete")[1:10]
met_hiroshige <- met.brewer("Hiroshige", 10, "discrete")
colour_model <- scale_colour_manual(values = c(met_hiroshige[4], met_hiroshige[7]))
fill_model <- scale_fill_manual(values = c(met_hiroshige[3], met_hiroshige[8]))
