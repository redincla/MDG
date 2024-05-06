#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
library(plotly)

args <- commandArgs(trailingOnly = TRUE)
data <- read.table(args, sep = "\t", fill=TRUE, header = TRUE)

# Select even-numbered columns, where the ratio is stored
ratios <- data[, seq(2, ncol(data), by = 2)]
SMN_copies<-2*(as.data.frame(t(ratios)))
SMN_copies<-SMN_copies[,1:8]

# perform k-means on ratios, using 7 clusters
predefined_centers<- matrix(c(
     1, 1, 1, 1, 2, 2, 2, 2,
     2, 2, 2, 2, 0, 0, 0, 0,
     2, 2, 2, 2, 1, 1, 1, 1,
     2, 2, 2, 2, 2, 2, 2, 2,
     2, 2, 2, 2, 3, 3, 3, 3,
     3, 3, 3, 3, 1, 1, 1, 1,
     3, 3, 3, 3, 2, 2, 2, 2), ncol = 8, byrow = TRUE)


#kmeans_result <- kmeans(SMN_copies, centers = predefined_centers, nstart=25)
kmeans_result <- kmeans(SMN_copies, 6, nstart=25)
clustered_data <- data.frame(PC1 = SMN_copies[,2], PC2 = SMN_copies[,6], Cluster = as.factor(kmeans_result$cluster))


# plot clustering displaying sample IDs
p <- ggplot(clustered_data, aes(x = PC1, y = PC2, color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN PCA plot with K-means Clustering", x = "SMN1 copies", y = "SMN2 copies", color = "Cluster") +
theme_minimal() + 
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2")

# save plot to PNG file
png("SMN_PCA_plot.png", width = 800, height = 600)
print(p)
dev.off()

# Save interactive plot to HTML file
p_plotly <- ggplotly(p, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_PCA_plot.html", selfcontained = TRUE)

# save tables to file
write.table(ratios, file = "ratios_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kmeans_result$centers, file = "kmeans_result_table.txt", sep = "\t", quote = FALSE, row.names = TRUE)