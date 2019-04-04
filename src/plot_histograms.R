library(data.table)
library(ggplot2)

hist_files <- list.files("output/040_stats",
                         pattern = "_readlength.txt",
                         full.names = TRUE)
names(hist_files) <- sub("_readlength.txt", "", basename(hist_files))

hist_data_list <- lapply(hist_files,  fread, skip = 9)
hist_data <- rbindlist(hist_data_list, idcol = "barcode")

hist_data[`#Length` == 400, col := 'a']
hist_data[`#Length` != 400, col := 'b']

ggplot(hist_data, aes(x = `#Length`, y = reads)) +
    scale_x_log10() +
    xlim(c(0, 1400)) +
    facet_wrap(~ barcode, scales = "free_y") +
    geom_vline(xintercept = 434) +
    geom_point(colour = alpha("black", 0.5)) +
    geom_smooth(n = 1000, se = FALSE, span = 0.1)

fs <- RColorBrewer::brewer.pal(3, "Set1")[1:2]

ggplot(hist_data, aes(x = `#Length`, y = reads, fill = col)) +
    scale_fill_manual(values = fs, guide = FALSE) +
    xlim(c(0, 1400)) +
    facet_wrap(~ barcode, scales = "free_y", ncol = 12) +
    geom_col(position = "identity")
