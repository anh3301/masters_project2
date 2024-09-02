# Author:       Anh Nguyen Phuong
# Description:  ggplot theme for thesis
#=========================================================================

library (ggplot2)
library (ggpubr)

theme_manuscript <- function(base_size=21, base_family="arial") {
  library(grid)
  (ggthemes::theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",hjust = 0.5,size = base_size),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title.y = element_text(angle=90,vjust =2,size = base_size),
            axis.title.x = element_text(vjust = -0.2,size = base_size),
            axis.text = element_text(size = base_size), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(size=base_size)
    ))
}
