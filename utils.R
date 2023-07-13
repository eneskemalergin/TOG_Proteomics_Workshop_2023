impute_with_downshifted_normal <- function(
    data,
    intensity_log2,
    prctl = 0.05,
    downshift_mag = 1.8,
    downshift_min = 0.1
) {

    # Find complete distribution
    complete_dist <- data %>%
        dplyr::filter(!is.na({{ intensity_log2 }})) %>%
        dplyr::pull({{ intensity_log2 }}) %>%
        as.vector()

    # Find the 5th percentile of the complete distribution
    downshift_threshold <- quantile(
        x = complete_dist,
        probs = c(prctl)
    )
    # Downshift the lower part of the whole distribution
    dist2choose <- (
        complete_dist[complete_dist <= downshift_threshold]
    ) / (2**downshift_mag)
    # Get mean from dist2choose
    mu <- mean(dist2choose)
    # Get the sd from the whole distribution
    sigma <- sd(complete_dist)
    # Impute missing values in the data with the downshifted normal
    data <- data %>%
        dplyr::mutate(
            imputed_intensity_log2 = {{ intensity_log2 }}
        ) %>%
        dplyr::mutate(
            imputed_intensity_log2 = ifelse(
                is.na(imputed_intensity_log2),
                rnorm(
                    n = sum(is.na(imputed_intensity_log2)),
                    mean = mu,
                    sd = sigma
                ),
                imputed_intensity_log2
            )
        ) %>%
        dplyr::mutate(
            imputed_intensity_log2 = ifelse(
                imputed_intensity_log2 < downshift_min,
                downshift_min,
                imputed_intensity_log2
            )
        )

    return(data)
}

plot_volcano <- function(
    data,
    pval_thr = 0.05,
    log2_fc_thr = 1,
    title = "Volcano Plot",
    point_size = 2
) {
    # NOTE: This function expects the following columns
    # in the data frame as they appear.
    # - log2_fc
    # - adj_pvalues
    # - significance


    p <- ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = log2_fc,
            y = -log10(adj_pvalues),
            color = significance,
            alpha = significance
        )
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_vline(
        xintercept = log2_fc_thr,
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::geom_vline(
        xintercept = -log2_fc_thr,
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::geom_hline(
        yintercept = -log10(pval_thr),
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "Up regulated" = "#e63946",
            "Down regulated" = "#1d3557",
            "no significance" = "#b1a7a6"
        )
    ) +
    ggplot2::scale_alpha_manual(
        values = c(
            "Up regulated" = 1.0,
            "Down regulated" = 1.0,
            "no significance" = 0.2
        )
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::labs(
        x = "log2(Fold-Change)",
        y = "-log10(adjusted p-value)"
    )
    return(p)
}