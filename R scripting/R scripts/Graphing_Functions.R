#########################
## 1. Load in packages ##
#########################
library( tidyverse )

mytheme <- theme(
  panel.background = element_rect( fill = "white" ),
  text = element_text( color = "black", face = "bold", family = "sans" ),
  axis.text = element_text( color = "black" ),
  axis.ticks = element_line( color = "black" ),
  plot.margin = unit( c( 0.25, 0.25, 0.25, 0.25 ), "cm" ),
  plot.title = element_text( vjust = 2 ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.y = element_text( vjust = 0 ),
  axis.title.x = element_text( vjust = 0 ),
  panel.border = element_blank(),
  axis.line = element_line()
)