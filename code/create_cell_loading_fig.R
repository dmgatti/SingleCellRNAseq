library(tidyverse)
library(gridExtra)

x = data.frame(mult = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6),
               load = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
               recov = c(500, 1:10 * 1000))


p1 = x %>%
      select(load, mult) %>%
      ggplot(aes(load, mult)) +
        geom_point() +
        geom_line() +
        labs('Number Cells Loaded vs. Multiplet Rate',
             x = 'Number of Cells Loaded',
             y = '% Multiplets') +
        lims(x = c(500, 16000))

p2 = x %>%
      select(load, recov) %>%
      ggplot(aes(load, recov)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
        geom_line() +
        labs('Number Cells Loaded vs. Cells Recovery',
             x = 'Number of Cells Loaded',
             y = 'Number of Cells Recovered') +
        lims(x = c(500, 16000), y = c(500, 16000))



png('~/classes/JAX/scRNAseq/SingleCellRNAseq/fig/cell_loading.png',
    width = 600, height = 1200, res = 192)
print(grid.arrange(p1, p2, nrow = 2))
dev.off()

