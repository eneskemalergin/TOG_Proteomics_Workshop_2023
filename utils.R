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
                is.na( imputed_intensity_log2 ),
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