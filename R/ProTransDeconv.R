#' Process and Report CV with Optional Deconvolution
#'
#' Computes CV across multiple normalization methods, generates a report,
#' and performs EDec and Rodeo deconvolution if cell type proportions are provided.
#'
#' @param data Expression matrix (genes × samples)
#' @param type One of "intensity", "ratio", or "spectra count"
#' @param marker_genes Optional vector of marker gene names
#' @param use_markers_only Logical, filter CV computation to marker_genes
#' @param cv_threshold Threshold for flagging high CV
#' @param output_html Path to save the generated HTML report
#' @param cell_proportion Optional cell type proportion matrix for deconvolution
#'
#' @return A list containing transformed data, CV summary, ridge plot data,
#'         deconvolution results, and marker gene tables
#' @export
process_and_report_cv <- function(data,
                                  type = c("intensity", "ratio", "spectra count"),
                                  marker_genes = NULL,
                                  use_markers_only = FALSE,
                                  cv_threshold = 0.25,
                                  output_html = NULL,
                                  cell_proportion = NULL) {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(ggridges)
        library(dplyr)
        library(knitr)
        library(rmarkdown)
    })
    
    type <- match.arg(type)
    
    if (is.null(output_html)) {
        output_html <- file.path(getwd(), "cv_summary_report.html")
    }
    
    dir_path <- dirname(output_html)
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    
    output_html <- tryCatch(
        normalizePath(output_html, mustWork = FALSE),
        error = function(e) output_html
    )
    
    compute_cv <- function(mat) {
        apply(mat, 1, function(x) {
            m <- mean(x, na.rm = TRUE)
            if (m <= 0) return(NA)
            sd(x, na.rm = TRUE) / m
        })
    }
    
    min_score_normalize <- function(data) {
        data <- as.matrix(data)
        t(apply(data, 1, function(x) {
            min_val <- min(x, na.rm = TRUE)
            sd_val <- sd(x, na.rm = TRUE)
            normed <- (x - min_val) / sd_val
            normed[is.na(x)] <- NA
            return(normed)
        }))
    }
    
    min_max_normalize <- function(data, a = 0, b = 1) {
        data <- as.matrix(data)
        apply(data, 2, function(x) {
            rng <- range(x, na.rm = TRUE)
            scaled <- (x - rng[1]) / (rng[2] - rng[1]) * (b - a) + a
            scaled[is.na(x)] <- NA
            return(scaled)
        })
    }
    
    logistic <- function(x) {
        result <- rep(NA, length(x))
        if (all(is.na(x))) return(result)
        a <- 1 / max(x, na.rm = TRUE)
        result[!is.na(x)] <- 1 - exp(-a * x[!is.na(x)])
        return(result)
    }
    
    tanh_scaled <- function(x) {
        result <- rep(NA, length(x))
        result[!is.na(x)] <- (tanh(x[!is.na(x)]) + 1) / 2
        return(result)
    }
    
    quantile_normalize <- function(mat) {
        if (!requireNamespace("preprocessCore", quietly = TRUE)) {
            stop("❌ Package 'preprocessCore' is required for quantile normalization.")
        }
        mat <- as.matrix(mat)
        normed <- preprocessCore::normalize.quantiles(mat)
        dimnames(normed) <- dimnames(mat)
        return(normed)
    }
    
    ratio_shift <- function(data) {
        data <- as.matrix(data)
        t(apply(data, 1, function(x) {
            min_val <- min(x, na.rm = TRUE)
            shifted <- x - min_val
            shifted[is.na(x)] <- NA
            return(shifted)
        }))
    }
    
    transform_methods <- list(
        min_score = min_score_normalize,
        min_max = min_max_normalize
    )
    
    if (type == "intensity") {
        transform_methods$logistic <- function(data) {
            dat <- as.matrix(data)
            original_colnames <- colnames(dat)
            original_rownames <- rownames(dat)
            
            dat <- t(apply(t(dat), 2, logistic))
            valid_rows <- rowSums(dat, na.rm = TRUE) > 0
            dat <- dat[valid_rows, , drop = FALSE]
            dat[!is.na(dat)] <- dat[!is.na(dat)] / max(dat, na.rm = TRUE)
            
            colnames(dat) <- original_colnames
            rownames(dat) <- original_rownames[valid_rows]
            return(dat)
        }
        transform_methods$quantile <- quantile_normalize
        transform_methods$inverse <- function(data) 2 ^ data
        transform_methods$original <- as.matrix
    }
    
    if (type == "ratio") {
        transform_methods$ratio_shift <- ratio_shift
        transform_methods$inverse <- function(data) 2 ^ data
        transform_methods$tanh <- function(data) {
            dat <- as.matrix(data)
            dat <- t(apply(dat, 1, tanh_scaled))
            return(dat)
        }
    }
    
    if (type == "spectra count") {
        transform_methods$logistic <- function(data) {
            dat <- as.matrix(data)
            original_colnames <- colnames(dat)
            original_rownames <- rownames(dat)
            
            dat <- t(apply(t(dat), 2, logistic))
            dat <- dat[rowSums(dat, na.rm = TRUE) > 0, , drop = FALSE]
            dat[!is.na(dat)] <- dat[!is.na(dat)] / max(dat, na.rm = TRUE)
            
            colnames(dat) <- original_colnames
            rownames(dat) <- original_rownames[rowSums(dat, na.rm = TRUE) > 0]
            return(dat)
        }
        transform_methods$quantile <- quantile_normalize
        transform_methods$original <- as.matrix
    }
    
    transformed_list <- lapply(transform_methods, function(f) f(data))
    
    compute_cv_summary <- function(mat) {
        if (use_markers_only && !is.null(marker_genes)) {
            mat <- mat[rownames(mat) %in% marker_genes, , drop = FALSE]
        }
        vals <- compute_cv(mat)
        vals[is.finite(vals)]
    }
    
    prop_col_name <- paste0("Prop_GT_", cv_threshold)
    
    cv_summary <- lapply(names(transformed_list), function(method) {
        vals <- compute_cv_summary(transformed_list[[method]])
        qtiles <- quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
        row <- data.frame(
            Method = method,
            Group = type,
            Mean_CV = round(mean(vals), 3),
            Median_CV = round(median(vals), 3),
            Min = round(min(vals), 3),
            Max = round(max(vals), 3),
            Q1 = round(qtiles[1], 3),
            Q3 = round(qtiles[2], 3),
            stringsAsFactors = FALSE
        )
        row[[prop_col_name]] <- round(mean(vals > cv_threshold), 3)
        row
    }) %>% bind_rows()
    
    rownames(cv_summary) <- seq_len(nrow(cv_summary))
    cv_df <- lapply(names(transformed_list), function(method) {
        vals <- compute_cv_summary(transformed_list[[method]])
        data.frame(Method = method, CV = vals, Group = type)
    }) %>% bind_rows() %>% na.omit()
    
    method_order <- cv_summary %>% arrange(desc(Median_CV)) %>% pull(Method)
    cv_df$Method <- factor(cv_df$Method, levels = method_order)
    cv_summary$Method <- factor(cv_summary$Method, levels = method_order)
    
    tmp_rmd <- tempfile(fileext = ".Rmd")
    writeLines(c(
        "---",
        "title: \"CV Summary Report\"",
        "output: html_document",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "library(ggplot2)",
        "library(ggridges)",
        "library(dplyr)",
        "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
        "```",
        "",
        "## 📊 CV Summary Table",
        "```{r}",
        "knitr::kable(cv_summary, digits = 3, row.names = FALSE)",
        "```",
        "",
        "## 📈 CV Ridge Plot",
        "```{r}",
        "ggplot(cv_df, aes(x = CV, y = Method, fill = Method)) +",
        "  geom_density_ridges(scale = 1.5, alpha = 0.65, rel_min_height = 0.01,",
        "                      color = 'gray20', size = 0.4) +",
        "  facet_grid(. ~ Group, scales = 'fixed') +",
        "  coord_cartesian(xlim = c(0, 4)) +",
        "  scale_fill_brewer(palette = 'Set2') +",
        "  theme_ridges(font_size = 18, grid = TRUE) +",
        "  labs(title = 'CV Distribution by Method',",
        "       x = 'Coefficient of Variation (CV)', y = NULL) +",
        "  theme(legend.position = 'none')",
        "```"
    ), tmp_rmd)
    
    rmarkdown::render(tmp_rmd, output_file = basename(output_html), output_dir = dirname(output_html), quiet = TRUE)
    message("✅ HTML Report Finished: ", output_html)
    
    edec_results <- NULL
    edec_marker_table <- NULL
    rodeo_results <- list()
    rodeo_marker_table <- list()
    
    compute_specificity <- function(mat) {
        mat <- as.matrix(mat)
        row_sums <- rowSums(mat, na.rm = TRUE)
        spec_mat <- sweep(mat, 1, row_sums, "/")
        spec_mat[!is.finite(spec_mat)] <- 0
        return(spec_mat)
    }
    
    if (!is.null(cell_proportion)) {
        edec_results <- deconvolve_transformed_list(transformed_list, cell_proportion)
        edec_marker_table <- lapply(names(edec_results), function(method) {
            res <- edec_results[[method]]
            if (is.null(res) || is.null(res$means)) return(NULL)
            spec <- compute_specificity(res$means)
            genes <- rownames(spec)
            cells <- colnames(spec)
            bin_mat <- matrix("no", nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
            for (cell in cells) bin_mat[spec[, cell] > 0.5, cell] <- "yes"
            colnames(spec) <- paste0(colnames(spec), "_specificity")
            colnames(bin_mat) <- paste0(colnames(bin_mat), "_marker")
            combined <- cbind(spec[rownames(bin_mat), , drop = FALSE], bin_mat)
            data.frame(Gene = rownames(combined), combined, row.names = NULL)
        })
        names(edec_marker_table) <- names(edec_results)
    }
    
    if (requireNamespace("Rodeo", quietly = TRUE)) {
        for (method in names(transformed_list)) {
            mat <- transformed_list[[method]]
            if (is.null(mat) || is.null(cell_proportion)) next
            try({
                result <- Rodeo::Rodeo(mat, t(cell_proportion))
                if (!is.null(result)) {
                    spec <- compute_specificity(result)
                    genes <- rownames(spec)
                    cells <- colnames(spec)
                    bin_mat <- matrix("no", nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
                    for (cell in cells) bin_mat[spec[, cell] > 0.5, cell] <- "yes"
                    colnames(spec) <- paste0(colnames(spec), "_specificity")
                    colnames(bin_mat) <- paste0(colnames(bin_mat), "_marker")
                    combined <- cbind(spec[rownames(bin_mat), , drop = FALSE], bin_mat)
                    marker_table <- data.frame(Gene = rownames(combined), combined, row.names = NULL)
                    rodeo_results[[method]] <- result
                    rodeo_marker_table[[method]] <- marker_table
                }
            }, silent = TRUE)
        }
    } else {
        warning("⚠️ Package 'Rodeo' is not installed. Skipping Rodeo deconvolution.")
    }
    
    return(invisible(list(
        transformed_list = transformed_list,
        cv_summary = cv_summary,
        cv_plot_data = cv_df,
        edec = list(edec_results = edec_results, edec_marker_table = edec_marker_table),
        rodeo = list(rodeo_results = rodeo_results, rodeo_marker_table = rodeo_marker_table)
    )))
}
