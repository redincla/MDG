#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
library(plotly)

args <- commandArgs(trailingOnly = TRUE)
ratios <- read.table(args, sep = "\t", fill=TRUE, header = TRUE)

SMN_copies<-2*(as.data.frame(t(ratios)))
#SMN_copies<-SMN_copies[,2:5]


#kmeans_result <- kmeans(SMN_copies, centers = predefined_centers, nstart=25)
nclus<- round(sqrt(nrow(SMN_copies)))-1
kmeans_result <- kmeans(SMN_copies, nclus, nstart=25)
clustered_data <- data.frame(PC1 = SMN_copies[,3], PC2 = SMN_copies[,1], Cluster = as.factor(kmeans_result$cluster))


# plot clustering displaying sample IDs
p1 <- ggplot(clustered_data, aes(x = SMN_copies[,3], y = SMN_copies[,1], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-i7 CN plot with K-means Clustering", x = "SMN1-i7 copies", y = "SMN2 copies-i7", color = "Cluster") +
theme_minimal() + 
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 3.5) + ylim(0, 3.5)
  
p2 <- ggplot(clustered_data, aes(x = SMN_copies[,4], y = SMN_copies[,2], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-E8 CN plot with K-means Clustering", x = "SMN1-E8 copies", y = "SMN2 copies-E8", color = "Cluster") +
theme_minimal() + 
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 3.5) + ylim(0, 3.5)
  
#p3 <- ggplot(clustered_data, aes(x = SMN_copies[,6], y = SMN_copies[,3], color = Cluster)) + geom_point(size = 3, alpha = 0.8) + geom_text(aes(label=rownames(clustered_data)), size = 3, nudge_y = 0.1) + labs(title = "SMN-i8 CN plot with K-means Clustering", x = "SMN1-i8 copies", y = "SMN2 copies-i8", color = "Cluster") +
#theme_minimal() + 
#  theme(
#    text = element_text(family = "sans"),
#    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
#    axis.title = element_text(face = "bold", size = 14),
#    legend.title = element_text(face = "bold", size = 12),
#    legend.text = element_text(size = 10)
#  ) + 
#  scale_color_brewer(palette = "Dark2") +
#  xlim(0, 3.5) + ylim(0, 3.5)
  

# save plot to PNG file
png("SMN_plot_i7.png", width = 800, height = 600)
print(p1)
png("SMN_plot_E8.png", width = 800, height = 600)
print(p2)
#png("SMN_plot_i8.png", width = 800, height = 600)
#print(p3)
dev.off()

# Save interactive plot to HTML file
p_plotly <- ggplotly(p1, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_i7.html", selfcontained = TRUE)
p_plotly <- ggplotly(p2, tooltip = "text")
htmlwidgets::saveWidget(p_plotly, "SMN_plot_E8.html", selfcontained = TRUE)
#p_plotly <- ggplotly(p3, tooltip = "text")
#htmlwidgets::saveWidget(p_plotly, "SMN_plot_i8.html", selfcontained = TRUE)

# save tables to file
write.table(ratios, file = "ratios_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kmeans_result$centers, file = "kmeans_result_table.txt", sep = "\t", quote = FALSE, row.names = TRUE)
