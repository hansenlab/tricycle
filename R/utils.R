

.gg_theme <- theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 0.25), axis.ticks = element_line(colour = "black", 
    linetype = 1, size = 0.25), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), 
    plot.title = element_text(face = "plain", size = 8, hjust = 0.5), plot.margin = unit(c(5, 5, 6, 6), "pt"), legend.background = element_blank(), legend.position = "right", 
    legend.justification = c(0.5, 0.5), legend.key.size = unit(7, "pt"), legend.key = element_blank(), legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"), legend.title = element_text(size = 5), 
    legend.title.align = 0.5, legend.text.align = 0.5, legend.text = element_text(size = 5), strip.text = element_text(color = "white"), strip.background = element_rect(fill = "black", 
        color = "black"))
