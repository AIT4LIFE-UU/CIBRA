#' edgeR_test: differential expression analysis of RNASeq data using edgeR
#'
#' @param data data matrix of counts data
#' @param gene_names vector of gene symbols to read to the result dataframe
#' @param group column name of the group of interest
#' @param reference reference group to use in the DE analysis
#'
#' @return EdgeR result object with the pvalues, adjusted pvalues and
#'   foldchanges of the DE analysis
#'
#'
#' @examples
#' # Internal function
#'
#' @import edgeR
#' @importFrom stats model.matrix
edgeR_test <- function(data, gene_names, group, reference){

  # select groups of interest from the design matrix
  group = factor(group[colnames(data),])
  group = relevel(group, ref = reference)

  y <- DGEList(data, group = group)
  # remember the gene name
  rownames(y$counts) <- gene_names

  # filter out low expressed genes
  keep <- filterByExpr(y, min.prop = 0.3)#, group=y$samples$group)
  y <- y[keep, , keep.lib.sizes=FALSE]
  # compute normalization
  y <- calcNormFactors(y)
  design <- model.matrix(~group)

  # fit model quasi-likelihood test
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit,coef=2)
  res <- topTags(qlf, n = Inf)
  res <- res$table

  colnames(res) <- c("log2FoldChange", "logCPM", "F",
                           "pvalue", "padj")

  return(res)
}

#' Limma_v_test: differential expression analysis of RNASeq data using Limma + Voom
#'
#' @param data data matrix of counts data
#' @param column column name of the groups to be used for the definitions
#' @param groups design matrix containing the definitions to be compared.
#' @param reference reference group to use in the DE analysis
#'
#' @return Limma-Voom result object with the pvalues, adjusted pvalues and
#'   foldchanges of the DE analysis
#'
#'
#' @examples
#' # Internal function
#'
#' @import limma
#' @importFrom stats relevel model.matrix
Limma_v_test <- function(data, column, groups, reference = "WT") {

  # select groups of interest from the design matrix
  group = factor(groups[colnames(data),])
  group = relevel(group, ref = reference)

  # create DGE object for analysis
  y <- DGEList(data, group = group)
  # calculate the normalization factors
  y <- calcNormFactors(y)

  # filter the genes by expression
  keep.exprs <- filterByExpr(y)
  y.filt <- y[keep.exprs,]

  # setup the model matrix
  des <- model.matrix(~ group, data=y.filt$samples)
  # fit the model
  y <- voom(y.filt, design = des)
  fit <- eBayes(lmFit(y, des))

  # extract the results
  res_table = topTable(fit, n=Inf, sort = "none")
  colnames(res_table) <- c("log2FoldChange", "AveExpr", "t",
                           "pvalue", "padj", "B")
  return(res_table)
}

#' DESeq2_test: differential expression analysis of RNASeq data using DESeq2
#'
#' @param data data matrix of counts data
#' @param column column name of the groups to be used for the definitions
#' @param groups design matrix containing the definitions to be compared.
#' @param reference Reference group to use in the DE analysis
#' @param covariates List of column names in the design matrix to be used as
#'   covariates
#' @param speedup Boolean parameter to indicate if the faster version of DESeq2
#'   should be used with fitType = "glmGamPoi" (note: it may not work depending
#'   on the system)
#' @param parallel Boolean parameter to indicate if DESeq2 should be run in
#'   parallel mode
#'
#' @return DESeq2 result object with the pvalues, adjusted pvalues and
#'   foldchanges of the DE analysis
#'
#'
#' @examples
#' # Internal function
#'
#' @import DESeq2
#' @import dplyr
#' @importFrom stats relevel formula
DESeq2_test <-
  function(data, column,
           groups,
           reference = "WT",
           covariates = c(),
           speedup = FALSE,
           parallel = FALSE) {
    # convert data to matrix
    data = as.matrix(data)
    # set groups as factor
    groups = dplyr::mutate_if(groups, is.character, as.factor)

    # create DESeq2 data matrix with the data and groups data
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                          colData = groups,
                                          design = formula(paste("~", column)))
    # filter out low count genes to speed up the process
    keep <- rowSums(DESeq2::counts(dds)) >= 10
    dds <- dds[keep, ]

    # set the reference as reference in the factor
    dds$condition <- factor(dds[[column]])
    dds$condition <- relevel(dds$condition, ref = reference)

    dds$condition <- droplevels(dds$condition)

    if (length(covariates) == 0) {
    DESeq2::design(dds) <- ~ condition
    } else {
      DESeq2::design(dds) <-
        formula(paste(c("~ condition", covariates), collapse = "+"))
    }

    # run DESeq2 in either normal or speedup mode
    if (speedup) {
      dds <- DESeq2::DESeq(dds, parallel = parallel, fitType = "glmGamPoi")
    } else {
      dds <- DESeq2::DESeq(dds, parallel = parallel)
    }
    # extract results
    res <- DESeq2::results(dds, contrast = c("condition",
                                             rev(levels(dds$condition))))
    return(res)
  }

#' Randomized CIBRA signal measure calculation
#'
#' @param data RNA count dataframe with genes as rows and samples as columns
#'   (dataframe)
#' @param n_cases number of cases in the data (num)
#' @param n_control number of controls in the data (num)
#' @param iterations number of iterations to run the permutation (num)
#' @param confidence confidence (Tau) for the proportion calculation (num)
#' @param case case definition (str)
#' @param control control definition (str)
#' @param covariates list of column names from the definition matrix to use as covariates (supported only with DESeq2)
#' @param covariate_matrix design dataframe of the covariates,
#'   columns to take along as covariate values and samples as rownames.
#' @param parallel boolean value indicating if the method should be run in
#'   parallel (boolean)
#' @param speedup boolean value if the DESeq2 sould be run in speedup mode
#'   (boolean)
#' @param column column name to give to the permutated sample column (string)
#' @param method DE analysis method to use (options: DESeq2, edgeR and
#'   limma-voom)
#' @param permutation permutatin appraoch to use, either sample or full (string)
#'
#' @return list containing 6 variables: the mean random proportion (float), the
#'   standard deviation of the calculated proportions (float), the mean random
#'   significant area (float), the standard deviation of the calculated
#'   significant area (float), signal_data: dataframe of the results containing
#'   the proportion, and significant area for each iteration, pvalue matrix of
#'   the differential expression analysis for each iteration as a column and
#'   genes as rows, adjusted pvalue matrix of the differential expression
#'   analysis for each iteration as a column and genes as rows, foldchange
#'   matrix of the differential expression analysis for each iteration as a
#'   column and genes as rows
#'
#' @examples
#' # Internal function
#'
#' @import dplyr
#' @importFrom stats sd
randomization <-
  function(data,
           n_cases,
           n_control,
           iterations,
           confidence,
           case,
           control,
           covariates = c(),
           covariate_matrix = NULL,
           parallel = FALSE,
           speedup = FALSE,
           column = "rand", method = "DESeq2",
           permutation = "sample") {
    # set state for first iteration to create the result dataframe
    state = TRUE
    # store the random signal measures
    proportion_rnd <- c()
    sign_area_rnd <- c()

    # store the genenames
    tryCatch({
      data_nogenename = dplyr::select(data, !(GeneName))
      genename = data$GeneName
    }, error = function(e){
      data_nogenename = data
      genename = rownames(data)
    }, finally = {
      data_nogenename = data
      genename = rownames(data)
    })

    # generate permutation signal measures for the number of iterations
    for (i in 1:iterations) {
      tryCatch({

        # select random cases and controls
        print(paste0("cases: ", n_cases, ", controls: ",
                     n_control, ", iteration: ", i))
        random_samples <- sample(colnames(data_nogenename))
        group <-
          as.data.frame(cbind(sample = random_samples[1:n_cases], group = case))
        group_tmp <-
          as.data.frame(cbind(sample = random_samples[(n_cases + 1):(n_control +
                                                                       n_cases)],
                              group = control))
        group <- dplyr::bind_rows(group, group_tmp)

        #ToDO: add the covariate_matrix to the group dataframe if it is not NULL
        # match covariate matrix with random group matrix
        if (!is.null(covariate_matrix)) {
          # limit the covariate matrix to only the samples in the group matrix
          covariate_matrix$sample = make.names(rownames(covariate_matrix))
          covariate_tmp = covariate_matrix[which(covariate_matrix$sample %in% group$sample),]
          covariate_tmp = covariate_tmp[match(covariate_tmp$sample, group$sample),]

          group = cbind(group, covariate_tmp)
        }

        rownames(group) <- group$sample

        data_filt <- dplyr::select(data, as.character(group$sample))

        # incase of full permutation permute also the genes
        if (permutation == "full") {
          rows <- sample(dim(data_filt)[1])
          shuffled_data <- data_filt[rows,]
          rownames(shuffled_data) <- genename
        } else if (permutation == "sample") {
          shuffled_data <- data_filt
        }

        group$sample <- NULL

        # perform the DE analysis with the chosen method
        if (tolower(method) == "deseq2") {
          # add option for covariates
          res <-
            DESeq2_test(shuffled_data,
                        "group",
                        group,
                        reference = control,
                        parallel = parallel,
                        covariates = covariates)
        } else if (tolower(method) == "limma") {
          res <- Limma_v_test(shuffled_data,
                              "group",
                              group,
                              reference = control)
        } else if (tolower(method) == "edgeR") {
          res <- edgeR_test(shuffled_data,
                            gene_names = genename,
                            group = group,
                            reference = control)
        }

        # store the results of the DE analysis
        pval_res <- data.frame(res$pvalue, row.names = rownames(res))
        colnames(pval_res) <- paste0(column, "_", i)
        padj_res <- data.frame(res$padj, row.names = rownames(res))
        colnames(padj_res) <- paste0(column, "_", i)
        fc_res <-
          data.frame(res$log2FoldChange, row.names = rownames(res))
        colnames(fc_res) <- paste0(column, "_", i)

        res_pval <- res$pvalue

        # filter non finite pvalues
        finite_values = is.finite(res_pval)
        res_pval = res_pval[finite_values]
        names(res_pval) = rownames(res[finite_values,])

        # calculate the signal measures
        sign_measures_rand <- signal_measures(res_pval, confidence)

        # store the results of the analysis
        if (state) {
          signal_data = data.frame(
            n_case = n_cases,
            n_control = n_control,
            iteration = i,
            proportion = sign_measures_rand[[1]],
            significant_area = sign_measures_rand[[2]]
          )
          pval_mat <- pval_res
          fc_mat <- fc_res
          padj_mat <- padj_res
          state = FALSE
        } else {
          signal_rand = data.frame(
            n_case = n_cases,
            n_control = n_control,
            iteration = i,
            proportion = sign_measures_rand[[1]],
            significant_area = sign_measures_rand[[2]]
          )
          signal_data = rbind(signal_data,  signal_rand)
          pval_mat = merge(pval_mat, pval_res, by = 0, all = TRUE)
          rownames(pval_mat) = pval_mat$Row.names
          pval_mat$Row.names = NULL
          fc_mat <- merge(fc_mat, fc_res, by = 0, all = TRUE)
          rownames(fc_mat) = fc_mat$Row.names
          fc_mat$Row.names = NULL
          padj_mat <- merge(padj_mat, padj_res, by = 0, all = TRUE)
          rownames(padj_mat) = padj_mat$Row.names
          padj_mat$Row.names = NULL
        }
      }, error = function(e) {
        print(e)
      })
    }

    # calculate the average measures
    proportion_rnd_avg <- mean(signal_data$proportion)
    proportion_rnd_sd <- sd(signal_data$proportion)
    sign_area_rnd_avg <- mean(signal_data$significant_area)
    sign_area_rnd_sd <- sd(signal_data$significant_area)

    return(
      list(
        proportion_rnd_avg,
        proportion_rnd_sd,
        sign_area_rnd_avg,
        sign_area_rnd_sd,
        signal_data,
        pval_mat,
        padj_mat,
        fc_mat
      )
    )

  }

#' Beta Uniform Mixture density function
#'
#' @param x x values from the fitted BUM model
#' @param lambda lambda from the fitted BUM model
#' @param a alpha from the fitted BUM model
#'
#' @return density value of the Beta Uniform Mixture density function
#'
#'
#' @examples
#' # Internal function
BUM_dens_sign <- function(x, lambda=lambda, a=a){
  # integral formula between the Beta and Uniform component of the Beta-Uniform
  # mixture model
  return((1-lambda)*a*(x^(a-1)-1))
}

#' calculate the CIBRA signal measures
#' Calculate the significant area and proportion from a p-value distribution
#'
#' @param pvalues vector of pvalues generated from DE analysis
#' @param confidence confidence value to use to calculate the proportion
#' @param gene vector of gene name
#' @param cond vector of condition name
#' @param outputDir path for directory of result visualizations
#'
#' @return list containing the poportion and significant area
#'
#'
#' @examples
#' # Internal function
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline
#' @importFrom stats integrate
#' @import BioNet
signal_measures <-
  function(pvalues,
           confidence,
           gene = NULL,
           cond = NULL,
           outputDir = FALSE) {
    tab_prop <-
      as.numeric(prop.table(table(factor(
        pvalues < confidence, levels = c(F, T)
      ))))
    pval <- as.matrix(pvalues)
    #rownames(pval) <- res$genename
    pval <- BioNet::aggrPvals(pval, plot = FALSE)
    # fit the beta uniform mixture model and save the parameters
    fb <- BioNet::fitBumModel(pval, plot = FALSE)
    a <- fb$a
    lambda <- fb$lambda
    # compute the significant area
    sign_area <-
      integrate(BUM_dens_sign, lambda, a, lower = 0, upper = 1)$value

    if (outputDir != FALSE) {
      dir.create(paste0(outputDir, gene), recursive = T)
      pdf(paste0(outputDir, gene, "/", gene,"_", cond, "_signal_plot.pdf"))
      BioNet::hist.bum(fb,
                       main = paste(
                         "Gene:",
                         gene,
                         "\nProportion of <=",
                         confidence,
                         ":",
                         (1 - round(tab_prop[1], 4)),
                         "\nSignificant area:",
                         (round(sign_area, 4))
                       ))
      abline(v = confidence, col = "red", lty = "dashed")
      dev.off()
    }

    # return some significance measures:
    # the proportion of pvalues > than a given confidence
    # the significant component of the BUM model
    return(list(1 - round(tab_prop[1], 4), round(sign_area, 4)))

  }

#' `scale_distance` calculates the similarity scores
#' Uses the diffential expression states to calculate the similarity and anti-similarity scores.
#'
#' @param data1 dataframe of the foldchange and adjusted pvalue of condition 1
#'   with the columns: fc (foldchange), pval (adjusted pvalue), genes (gene
#'   identifiers)
#' @param data2 data1 for condition 2
#' @param genes vector of the intersect of genes between data1 and data2
#'
#' @return list containing the similarity scores dplus and dmin and the corresponding genes associated with the dplus and dmin scores, dplusgenes and dmingenes.
#'
#'
#' @examples
#' # Internal function
scale_distance <- function(data1, data2, genes){
  # calculate similarity and anti-similarity distance
  Dplus = 0 # similarity distance
  Dplus_genes = c()
  Dmin = 0 # anti-similarity distance
  Dmin_genes = c()
  gene_n = length(data1)

  # iterate through the genes
  for (i in 1:gene_n) {
    # if the signs of the data scores are equal it is a similarity distance
    if (data1[i] != 0 & data2[i] != 0) {
      if (sign(data1[i]) == sign(data2[i])) {
        # if the scores are equal add 1 (diagonal)
        if (data1[i] == data2[i]) {
          Dplus = Dplus + 1
          Dplus_genes = c(Dplus_genes, genes[i])
          # if the scores are deviating by one neighboring cell add 0.5
          # (next to the diagonal)
        } else if (abs(data1[i] - data2[i]) == 0.5) {
          Dplus = Dplus + 0.5
          Dplus_genes = c(Dplus_genes, genes[i])
        }
        # repeat with the same logic for the anti diagonal
      } else {
        if (abs(data1[i]) == abs(data2[i])) {
          Dmin = Dmin - 1
          Dmin_genes = c(Dmin_genes, genes[i])
        } else if (abs(abs(data1[i]) - abs(data2[i])) == 0.5) {
          Dmin = Dmin - 0.5
          Dmin_genes = c(Dmin_genes, genes[i])
        }
      }
    }
  }

  return(list(dplus = Dplus, dmin = Dmin, dplusgenes = Dplus_genes,
              dmingenes = Dmin_genes))
}

#' Assign score to the genes depending on in which box they lay in the volcano plot
#'
#' @param pval vector of adjusted pvalues from DE analysis
#' @param fc vector of foldchanges from DE analysis
#'
#' @return vector of score assignments
#'
#'
#' @examples
#' # Internal function
score_assignment <- function(pval, fc) {
  if (pval > -log10(0.05)) {
    # check if box 1 or box 2
    if (fc < -1) {
      # box 1
      return(-1)
    } else if (fc > 1) {
      # box 2
      return(1)
    } else if (fc > 0) {
      # indeterminate_plus_box
      return(0.5)
    } else {
      # indeterminate_negative_box
      return(-0.5)
    }
  } else if (pval > -log(0.1)) {
    if (fc > 0) {
      # indeterminate_plus_box
      return(0.5)
    } else {
      # indeterminate_negative_box
      return(-0.5)
    }
  } else {
    # no box
    return(0)
  }
}

#' `dissimilarity_correlation`
#'
#' @param data1 dataframe of the foldchange and adjusted pvalue of condition 1
#'   with the columns: fc (foldchange), pval (adjusted pvalue), genes (gene
#'   identifiers)
#' @param data2 the same as data1 for condition 2
#' @param outputName filename for the output visualization
#' @param cond1 name of condition 1
#' @param cond2 name of condition 2
#' @param perm_dist vector of permutation of the dissimilarity scores
#' @param perm_dplus vector of permutation of the dplus values
#' @param perm_dmin vector of permutation of the dmin values
#' @param perm_mode if TRUE, permutation mode is used to calculate the
#'   permutation distributions. No need to supply the permutation distributions.
#'   No visualizations will be made
#'
#' @return list of the correlation measure (rho), the similarity score (dplus)
#'   and the anti-similarity score (dmin).
#'
#'
#' @examples
#' # Internal function
#'
#' @importFrom stats cor
#' @importFrom grDevices dev.off pdf
#' @import tidyr
#' @import tibble
#' @import grid
dissimilarity_correlation <- function(data1, data2, outputName = "",
                                      cond1 = "data1",
                                      cond2 = "data2",
                                      perm_dist = NULL, perm_dplus = NULL,
                                      perm_dmin = NULL, perm_mode = FALSE) {

  # remove na values and calculate the score assignment
  data1 = data1 %>% filter(!(is.na(fc) | is.na(pval))) %>% rowwise() %>%
    mutate(score_assignment = score_assignment(pval, fc))
  data2 = data2 %>% filter(!(is.na(fc) | is.na(pval))) %>% rowwise() %>%
    mutate(score_assignment = score_assignment(pval, fc))

  # select overlapping genes
  rownames(data1) <- data1$genes
  rownames(data2) <- data2$genes
  genes <- intersect(rownames(data1), rownames(data2))

  # exclude no box genes in both conditions (not significant in both conditions)
  no_box <- data1[genes,]$score_assignment == 0 &
    data2[genes,]$score_assignment == 0
  genes <- genes[!no_box]

  # calculate the spearman correlation measure between the score assignments between the two conditions
  spearman_rho <- cor(data1[intersect(data1$genes,
                                      data2$genes),]$score_assignment,
                      data2[intersect(data1$genes,
                                      data2$genes),]$score_assignment,
                      method = "spearman")

  # calculate similarity count measure for both positive and negative similarity
  scale_d = scale_distance(data1[intersect(data1$genes,
                                           data2$genes),]$score_assignment,
                           data2[intersect(data1$genes,
                                           data2$genes),]$score_assignment,
                           intersect(data1$genes, data2$genes))

  if (!perm_mode) {

    # calculate pvalues
    perm_pval = sum(abs(spearman_rho) < perm_dist) / length(perm_dist)
    perm_dvalplus = sum(scale_d$dplus < perm_dplus) / length(perm_dplus)
    perm_dvalmin = sum(scale_d$dmin > perm_dmin) / length(perm_dmin)

    # volcano plot of one of the conditions colored by the Dplus genes and Dmin genes

    # colors
    blue = "#998ec3"
    red = "#f1a340"

    # set parameters for the plot
    plot_data = data1[intersect(data1$genes, data2$genes),]
    plot_data$color = "-"
    plot_data$alpha = 0.5
    plot_data[plot_data$genes %in% scale_d$dplusgenes, "color"] = "D+"
    plot_data[plot_data$genes %in% scale_d$dmingenes, "color"] = "D-"
    plot_data[plot_data$genes %in% scale_d$dplusgenes, "alpha"] = 1
    plot_data[plot_data$genes %in% scale_d$dmingenes, "alpha"] = 1

    fc_lim = abs(max(plot_data$fc))

    # create dataframe for plotting
    plot_data <- data.frame(data1 = data1[intersect(data1$genes,
                                                    data2$genes),]$score_assignment,
                              data2 = data2[intersect(data1$genes,
                                                    data2$genes),]$score_assignment)
    plot_data <- plot_data %>% group_by(data1, data2) %>% summarise(events = n())

    # set parameters for the plot
    arrow_size = 0.1*1
    arrow_space = 0.3*arrow_size
    max_data1 <- max_data2 <- 0.9
    min_data1 <- min_data2 <- -0.9

    # reformat data for the heatmap
    mat_data <- plot_data %>% mutate(events = log(events)) %>%
    tidyr::pivot_wider(names_from = data2, values_from = events) %>%
    tibble::column_to_rownames("data1")

    # color range
    f1 <- colorRamp2(seq(min(mat_data, na.rm = T), max(mat_data, na.rm = T),
                         length = 2),
                      c("#f7fcfd", "#377eb8"))

    # rounding summary measures of the similarity between the two conditions
    plot_spearman = round(spearman_rho, 3)
    plot_pval = round(perm_pval, 4)
    plot_pdplus = round(perm_dvalplus, 4)
    plot_pdmin = round(perm_dvalmin, 4)

    # heatmap of the correlation pattern between the two conditions
    mat_data = mat_data[c("1", "0.5", "0", "-0.5", "-1"),
                      c("-1", "-0.5", "0", "0.5", "1")]
    rownames(mat_data) <- c("HU", "MU", "NS", "MD", "HD")
    colnames(mat_data) <- c("HD", "MD", "NS", "MU", "HU")
    try({
      ht_list = Heatmap(mat_data,
                        name = "mat", na_col = "white", cluster_rows = F, col = f1,
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if (i == 3 | j == 3) {
                            if (!is.na(mat_data[i,j])) {
                              grid::grid.points(x, y, pch = 21, size = unit(12, "mm"),
                                          gp = grid::gpar(col = "black", fill = "darkgrey", lwd = 3,
                                                    alpha = .7))
                            }
                          } else if(i == j) {
                            if (!is.na(mat_data[i,j])) {
                              grid::grid.points(x, y, pch = 21, size = unit(12, "mm"),
                                          gp = grid::gpar(col = "black", fill = "#984ea3", lwd = 3,
                                                    alpha = .7))
                            }
                          } else if ((i == j+1 | j == i+1)) {
                            if (!is.na(mat_data[i,j])) {
                              grid::grid.points(x, y, pch = 21, size = unit(12, "mm"),
                                          gp = grid::gpar(col = "black", fill = "#984ea3", lwd = 3,
                                                    alpha = .4))
                            }
                          } else if (i == (6-j)) {
                            if (!is.na(mat_data[i,j])) {
                              grid::grid.points(x, y, pch = 21, size = unit(12, "mm"),
                                          gp = grid::gpar(col = "black", fill = "#ff7f00", lwd = 3,
                                                    alpha = .7))
                            }
                          } else if (i == (6-j+1) | (6-j) == (i+1)) {
                            if (!is.na(mat_data[i,j])) {
                              grid::grid.points(x, y, pch = 21, size = unit(12, "mm"),
                                          gp = grid::gpar(col = "black", fill = "#ff7f00", lwd = 3,
                                                    alpha = .4))
                            }
                          }


                        },
                        cluster_columns = F,
                        row_title = cond1, column_title = cond2, row_names_side = "left",
                        column_title_side = "bottom", column_names_rot = 0,
                        rect_gp = grid::gpar(col = "black", lwd = 3),
                        heatmap_legend_param = list(title_position = "leftcenter-rot",
                                                    title = "log(nr. events)",
                                                    title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                                                    labels_gp = grid::gpar(fontsize = 12, fontface = "bold")),
                        row_title_gp = grid::gpar(fontsize = 22, fontface =
                                              "bold"),
                        column_title_gp = grid::gpar(fontsize = 22, fontface = "bold"),
                        row_names_gp = grid::gpar(fontsize = 16, fontface = "bold"),
                        column_names_gp = grid::gpar(fontsize = 16, fontface = "bold"))

        pdf(paste0(outputName, "heatmap.pdf"))
        draw(ht_list, padding = unit(c(2, 2, 10, 2), "mm"),
            column_title = expr(paste(bolditalic(rho), ": ", !!plot_spearman,
                                     " D+: ", !!scale_d$dplus,
                                     " p-value: ", !!plot_pdplus,
                                     " D-: ", !!scale_d$dmin,
                                     " p-value: ", !!plot_pdmin, sep="")))
        dev.off()
      })
    }

    return(list(rho = spearman_rho, dplus = scale_d$dplus,
                dmin = scale_d$dmin, dplus_genes = scale_d$dplusgenes,
                dmin_genes = scale_d$dmingenes))
}

#' Volcano differential expression state plot
#'
#' @param fc vector of foldchanges from condition 1 and 2
#' @param pval vector of adjusted p-values from condition 1 and 2
#' @param link_column gene names
#' @param condition_column vector of condition status of the rows
#' @param state_column Differential expression states from the
#'   dissimilarity_correlation function
#' @param score score to include as annotation for the plot
#' @param pvalue pvalue to include as annotation for the plot
#'
#' @return generates a figure
#'
#' @examples
#' # Internal function
volcano_de_state_plot <- function(fc, pval, link_column, condition_column, state_column, score = NULL, pvalue = NULL) {

  # create dataframe for plotting
  plot_data = data.frame(fc = fc, pval = pval, link_column = link_column,
                         condition_column = condition_column, state_column = state_column)
  plot_data = plot_data[sample(1:nrow(plot_data)), ]

  # remove na values
  plot_data = plot_data[!is.na(plot_data$fc) & !is.na(plot_data$pval),]


  # create a volcano plot with links between the two conditions
  fc_lim = abs(max(plot_data$fc))

  p <- ggplot(plot_data, aes(x = fc, y = pval, fill = condition_column)) +
    geom_line(data = plot_data[plot_data$state_column == "DEdiss",],
              aes(group = link_column, color = state_column, alpha = state_column), size = 1) +
    geom_point(aes(size = pval), alpha = .6, shape = 21, color = "black") +
    geom_line(data = plot_data[plot_data$state_column %in% c("DEplus"),],
              aes(group = link_column, color = state_column, alpha = state_column), size = 1) +
    geom_line(data = plot_data[plot_data$state_column %in% c("DEmin"),],
              aes(group = link_column, color = state_column, alpha = state_column), size = 1) +
    xlab(bquote(~bold(log[2] ~ '(fold change)'))) +
    ylab(bquote(~bold(-log[10] ~"(adjusted p-value)"))) + #ggprism::theme_prism() +
    xlim(c(-1*fc_lim, fc_lim)) +
    scale_size_continuous(range = c(1, 2.5)) +
    scale_color_manual(name = "similarity state", values = c("DEdiss" = "black", "DEplus" = "#f1a340",
                                                             "DEmin" = "#998ec3"),
                       labels = c("DEdiss" = "Dissimilarity", "DEplus" = "Similarity",
                                  "DEmin" = "Anti-similarity")) +
    scale_alpha_manual(name = "similarity state", values = c("DEdiss" = 0.3, "DEplus" = .5, "DEmin" = .5),
                       labels = c("DEdiss" = "Dissimilarity", "DEplus" = "Similarity",
                                  "DEmin" = "Anti-similarity")) +
    scale_fill_discrete(na.translate = F) +
    labs(fill = "condition")  +
    theme(axis.text = element_text(face="bold", size = 16),
          axis.title = element_text(face="bold", size = 22),
          legend.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 14),
          panel.grid.major = element_line(colour = "lightgrey"),
          legend.position = "top",
          legend.box = "vertical",
          legend.spacing.y = unit(-0.5, "cm")) +
    guides(size = "none",
           fill = guide_legend(override.aes = list(size=3)))
  if (!is.null(score)) {
    if (!is.null(pvalue)) {
      p = p + geom_text(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0("score: ", score, ", pval: ", pvalue),
                        color = "black")
    } else {
      p = p + geom_text(x=Inf,y=Inf,hjust=1,vjust=1,label=paste0("score: ", score), color = "black")
    }
  }
  print(p)
  return(p)
}
