library(ggplot2)

# ggsavae alias to automatically save the plot into png, pdf, svg, and tiff
ggsave2 <- function(filename, plot, width = 7, height = 7, dpi = 300) {
  ggsave(
    paste0(filename, ".png"),
    plot, width = width, height = height, dpi = dpi
  )
  ggsave(
    paste0(filename, ".svg"),
    plot, width = width, height = height, dpi = dpi
  )
}