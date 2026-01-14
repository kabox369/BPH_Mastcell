#Figure S5A
p <- ggplot(proportion_plot_primary, aes(x = Name, y = Value, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) + 
  facet_grid(
    . ~ Group, 
    scales = "free_x", 
    space = "free_x", 
    labeller = labeller(Group = label_value)
  ) +
  scale_fill_manual(
    values = my_cols_primary,
    breaks = legend_primary_order,  
    drop = FALSE           
  ) +
  labs(
    title = "Primary Cell Type Composition by Samplename",
    x = "Samplename",
    y = "Percentage (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),  
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14),  
    panel.spacing = unit(1, "lines")  
  )
ggsave(
  filename = "cell_type_composition_primary_type_Stacked_Bar_Chart.pdf",
  plot = p,
  width = 15,
  height = 12,
  dpi = 300
)

#Figure S5B
p <- ggplot(proportion_plot_secondary, aes(x = Name, y = Value, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) + 
  facet_grid(
    . ~ Group, 
    scales = "free_x", 
    space = "free_x", 
    labeller = labeller(Group = label_value)
  ) +
  scale_fill_manual(
    values = my_cols_secondary,
    breaks = legend_secondary_order,  
    drop = FALSE           
  ) +
  labs(
    title = "Cell Type Composition by Samplename",
    x = "Samplename",
    y = "Percentage (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),  
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14),  
    panel.spacing = unit(1, "lines")  
  )
ggsave(
  filename = "cell_type_composition_secondary_type_Stacked_Bar_Chart.pdf",
  plot = p,
  width = 15,
  height = 12,
  dpi = 300
)

#Figure S5C
p <- ggplot(group_proportion_secondary, 
            aes(x = Group, y = Percent, fill = Celltype)) +
  geom_col(position = "stack", width = 0.7, 
           color = "white", linewidth = 0.2, key_glyph = "rect") +
  coord_flip() +
  scale_fill_manual(values = my_cols_secondary, 
                    breaks = legend_secondary_order,
                    guide = guide_legend(
                      reverse = TRUE,
                      keywidth = 0.5,  
                      keyheight = 0.8   
                    )) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = scales::label_number(suffix = "%")) +
  labs(
    title = "Secondary Cell Type Composition Across Groups",
    x = NULL,
    y = NULL,
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(color = "black", size = 10, family = "Arial"),
    axis.text.y = element_text(color = "black", size = 12, face = "bold", margin = margin(r = 10)),  
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),  
    legend.spacing.y = unit(0.05, "cm"),  
    legend.title = element_text(face = "bold", size = 9),  
    legend.text = element_text(size = 8),  
    legend.margin = margin(0, 0, 0, -5),  
    plot.margin = margin(1, 1.5, 1, 1, "cm"),  
    axis.line.x = element_line(color = "black", linewidth = 0.5)
  )
ggsave(
  "Secondary Cell Type Composition Across Groups.pdf", 
  plot = p,
  width = 18,  
  height = 6,   
  device = cairo_pdf,
  dpi = 300
)