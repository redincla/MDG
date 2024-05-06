#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
library(plotly)

args <- commandArgs(trailingOnly = TRUE)
ratios <- read.table(args, sep = "\t", fill=TRUE, header = TRUE)

SMN_copies<-2*(as.data.frame(t(ratios)))
#SMN_copies<-SMN_copies[,4:11]

# perform k-means on ratios, using 7 clusters
predefined_centers<- matrix(c(
     2, 2, 2, 2, 1, 1, 1, 1, 
     0, 0, 0, 0, 2, 2, 2, 2, 
     1, 1, 1, 1, 2, 2, 2, 2, 
     2, 2, 2, 2, 2, 2, 2, 2, 
     3, 3, 3, 3, 2, 2, 2, 2, 
     1, 1, 1, 1, 3, 3, 3, 3, 
     2, 2, 2, 2, 3, 3, 3, 3), ncol = 8, byrow = TRUE)


#kmeans_result <- kmeans(SMN_copies, centers = predefined_centers, nstart=25)
nclus<- round(sqrt(nrow(SMN_copies)))-1
kmeans_result <- kmeans(SMN_copies, nclus, nstart=25)
clustered_data <- data.frame(PC1 = SMN_copies[,6], PC2 = SMN_copies[,2], Cluster = as.factor(kmeans_result$cluster))


# plot clustering displaying sample IDs
p1 <- ggplot(clustered_data, aes(x = SMN_copies[,5], y = SMN_copies[,1], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-i7 CN plot with K-means Clustering", x = "SMN1-i7 copies", y = "SMN2 copies-i7", color = "Cluster") +
theme_minimal() + 
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2")
  
p2 <- ggplot(clustered_data, aes(x = SMN_copies[,6], y = SMN_copies[,2], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-E8 CN plot with K-means Clustering", x = "SMN1-E8 copies", y = "SMN2 copies-E8", color = "Cluster") +
theme_minimal() + 
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2")
  
p3 <- ggplot(clustered_data, aes(x = SMN_copies[,7], y = SMN_copies[,3], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-i8a CN plot with K-means Clustering", x = "SMN1-i8a copies", y = "SMN2 copies-i8a", color = "Cluster") +
theme_minimal() + 
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2")
  
p4 <- ggplot(clustered_data, aes(x = SMN_copies[,8], y = SMN_copies[,4], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-i8b CN plot with K-means Clustering", x = "SMN1-i8b copies", y = "SMN2 copies-i8b", color = "Cluster") +
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
png("SMN_plot_i7.png", width = 800, height = 600)
print(p1)
png("SMN_plot_E8.png", width = 800, height = 600)
print(p2)
png("SMN_plot_i8a.png", width = 800, height = 600)
print(p3)
png("SMN_plot_i8b.png", width = 800, height = 600)
print(p4)
dev.off()

# Save interactive plot to HTML file
p_plotly <- ggplotly(p1, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_i7.html", selfcontained = TRUE)
p_plotly <- ggplotly(p2, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_E8.html", selfcontained = TRUE)
p_plotly <- ggplotly(p3, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_i8a.html", selfcontained = TRUE)
p_plotly <- ggplotly(p4, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_i8b.html", selfcontained = TRUE)

# save tables to file
write.table(ratios, file = "ratios_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kmeans_result$centers, file = "kmeans_result_table.txt", sep = "\t", quote = FALSE, row.names = TRUE)
