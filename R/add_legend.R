add_legend <- function(...) {
  opar <- graphics::par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(graphics::par(opar))
  graphics::plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  graphics::legend(...)
}