#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
library(plotly)

args <- commandArgs(trailingOnly = TRUE)
data <- as.matrix(read.table(args, sep = "\t", header = TRUE))
GBA_all_copies<-2*(as.data.frame(data))

GBA_table<-t(GBA_all_copies)

# perform k-means on ratios, using 2 clusters
kmeans_result <- kmeans(GBA_table, 2, nstart=25)
clustered_data <- data.frame(PC1 = GBA_table[,1], PC2 = GBA_table[,36], Cluster = as.factor(kmeans_result$cluster))

# plot clustering displaying sample IDs
p1 <- ggplot(clustered_data, aes(x = GBA_table[,1], y = GBA_table[,36], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "GBA loci CN plot with K-means Clustering", x = "GBAP1 copies", y = "GBA copies", color = "Cluster") +
	theme_minimal() +
	theme(
		text = element_text(family = "sans"),
		plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
		axis.title = element_text(face = "bold", size = 14),
		legend.title = element_text(face = "bold", size = 12),
		legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2")
  

# plot GBA loci CN copies 
colors <- c("blue", "red", "green", "purple", "orange", "cyan", "magenta", "black","darkgreen")
line_types <- c(1, 2, 3, 4)
p<-matplot(GBA_all_copies, type="l", xlab="", ylab="copy number", main="GBA/GBAP1 copies", col=colors, lty=line_types, lwd=2, ylim=c(0, 4),xaxt="n")
legend(legend=colnames(GBA_all_copies),col=colors, lty=line_types, lwd=2,x="topright",cex=0.8, ncol=3)


# save plot to PNG file
png("SMN_plot_i7.png", width = 800, height = 600)
print(p1)
png("SMN_plot_E8.png", width = 800, height = 600)
print(p2)
dev.off()

# Save interactive plot to HTML file
p_plotly <- ggplotly(p1, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_i7.html", selfcontained = TRUE)
p_plotly <- ggplotly(p2, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_E8.html", selfcontained = TRUE)

# save tables to file
write.table(ratios, file = "ratios_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kmeans_result$centers, file = "kmeans_result_table.txt", sep = "\t", quote = FALSE, row.names = TRUE)
