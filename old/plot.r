library(ggplot2)
plot_loadings <- function(csv_path, plot_name, trait_levs, cat_levs){
    df  <- read.csv(csv_path, header=TRUE)

    df$trait <- factor(df$trait, levels=trait_levs)
    df$cat <- factor(df$cat, levels=cat_levs)
    df$L <- sapply(df$L, as.numeric)
    ggplot(df, aes(x=trait, y=row, fill=L)) +
            facet_grid(~ cat, scales="free_x", space="free_x") +
            geom_tile() +
            scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
            scale_x_discrete(position = "top") +
            scale_y_discrete(limits=fact) +
      labs(y="Loadings Row", x="Trait") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust=0),
            strip.text.x = element_blank()
            )
            # axis.title.y = element_text())

    ggsave(plot_name, width=15, height=10, units="in")
    gc()
}

plot_factors <- function(csv_path, plot_name){
    df  <- read.csv(csv_path, header=TRUE)

    ggplot(df, aes(x=f1, y=f2, color=class)) +
            geom_point() +
            labs(y="f2", x="f1") +
      theme_minimal()

    ggsave(plot_name, width=6, height=4, units="in")

gc()
}

csv_path <- "old_test.csv"
data_path <- "../data/yeast_continuous.csv"
plot_name <- "old_test.pdf"
labels_df <- read.csv("../data/yeast_labels.csv")
fact <- 1:2

plot_loadings(csv_path, plot_name, labels_df$pretty, unique(labels_df$cat))
