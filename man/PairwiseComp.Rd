% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_PairwiseComp.r
\name{PairwiseComp}
\alias{PairwiseComp}
\title{Generate Crossed Lists from Nested Results}
\description{
A class representing a nested representation of pairwise
comparisons given a pairwise design.
}
\details{
Pairwise Comparison Class

This class serves as a superclass for all methods of pairwise comparisons
considered in this package. It defines the shared public and private methods
common to all subclass.

Pairwise comparisons are used in experimental designs with test vs control
samples, comparing their expression tags to determine deregulation between
conditions.

This function filters the results table based on specified groups and
batches. Note that the term "filter" in this context does not refer to a
statistical filter.
}
\section{Attributes}{

\itemize{
\item \code{expr_data}: An ExprData object.
\item \code{results}: A nested data.frame with columns "batch", "group", and a nested
"data" column containing:
\itemize{
\item baseMean: Numeric vector with the base mean of expression value
\item log2FoldChange: Numeric vector with the log2FC
\item lfcSE: Numeric vector with the log2FC standard error
\item pvalue: Numeric vector with the p-values
\item padj: Numeric vector with the adjusted p-values
\item status: Factor vector with four levels: undetected, filtered, outlier,
analyzed
\item tag_id: Character vector with the tag id
}
\item \code{opts}: Important options defined during initialization.
}
}

\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-PairwiseComp-get_results}{\code{PairwiseComp$get_results()}}
\item \href{#method-PairwiseComp-filter_and_get_results}{\code{PairwiseComp$filter_and_get_results()}}
\item \href{#method-PairwiseComp-generate_a_list}{\code{PairwiseComp$generate_a_list()}}
\item \href{#method-PairwiseComp-cross_args_and_generate_lists}{\code{PairwiseComp$cross_args_and_generate_lists()}}
\item \href{#method-PairwiseComp-write_to_xlsx}{\code{PairwiseComp$write_to_xlsx()}}
\item \href{#method-PairwiseComp-generate_summary}{\code{PairwiseComp$generate_summary()}}
\item \href{#method-PairwiseComp-plot_summary}{\code{PairwiseComp$plot_summary()}}
\item \href{#method-PairwiseComp-extract_data_for_plot}{\code{PairwiseComp$extract_data_for_plot()}}
\item \href{#method-PairwiseComp-plot_ma}{\code{PairwiseComp$plot_ma()}}
\item \href{#method-PairwiseComp-plot_vulcano}{\code{PairwiseComp$plot_vulcano()}}
\item \href{#method-PairwiseComp-plot_lfc_per_group}{\code{PairwiseComp$plot_lfc_per_group()}}
\item \href{#method-PairwiseComp-plot_lfc_per_group_facet_tags}{\code{PairwiseComp$plot_lfc_per_group_facet_tags()}}
\item \href{#method-PairwiseComp-plot_de}{\code{PairwiseComp$plot_de()}}
\item \href{#method-PairwiseComp-build_tag_hclust}{\code{PairwiseComp$build_tag_hclust()}}
\item \href{#method-PairwiseComp-plot_heatmap}{\code{PairwiseComp$plot_heatmap()}}
\item \href{#method-PairwiseComp-clone}{\code{PairwiseComp$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-get_results"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-get_results}{}}}
\subsection{Method \code{get_results()}}{
Get results
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$get_results()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A nested dataframe.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-filter_and_get_results"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-filter_and_get_results}{}}}
\subsection{Method \code{filter_and_get_results()}}{
Retrieve the results table after filtering for specific
groups/batches and with the option to add new columns to extend the tag
ID.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$filter_and_get_results(
  in_batch,
  in_group,
  verbose = FALSE,
  add_ids = c("gid", "symbol", "uniprot", "protein_names", "type", "tax_name"),
  safe_translate = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{in_batch}}{Vector of batch codes; select or report only samples in these
batches}

\item{\code{in_group}}{Vector of group codes; select or report only samples in these
groups}

\item{\code{verbose}}{If TRUE, the column names are replaced with more
meaningful names}

\item{\code{add_ids}}{Character vector; see \code{to} argument of method
\code{generate_translate_dict} from class \link{Annotation}
Available options are: "gid", "symbol", "uniprot", "protein_names",
"type", and "tax_name"}

\item{\code{safe_translate}}{If TRUE, the function will keep original ids
if the translation fails}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.frame containing filtered results with optional additional
ID columns
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-generate_a_list"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-generate_a_list}{}}}
\subsection{Method \code{generate_a_list()}}{
Generate a list of deregulated genes for a specific comparison based
on given criteria.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$generate_a_list(
  in_batch,
  in_group,
  id = "symbol",
  use_padj = TRUE,
  type = "deregulated",
  lfc_abs_lim = 1,
  min_signif = 0.05
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{in_batch}}{Vector of batch codes; select or report only samples in these
batches}

\item{\code{in_group}}{Vector of group codes; select or report only samples in these
groups}

\item{\code{id}}{Character vector; see \code{to} argument of method
\code{generate_translate_dict} from class \link{Annotation}}

\item{\code{use_padj}}{Whether to use the adjusted p-value or not}

\item{\code{type}}{Type of deregulation: "deregulated", "upregulated", or
"downregulated"}

\item{\code{lfc_abs_lim}}{Threshold for log fold change}

\item{\code{min_signif}}{Threshold for significance}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character vector
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-cross_args_and_generate_lists"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-cross_args_and_generate_lists}{}}}
\subsection{Method \code{cross_args_and_generate_lists()}}{
Generate crossed lists of genes from all comparisons in
nested results using multiple filtering criteria.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$cross_args_and_generate_lists(
  cross_id = c("symbol", "uniprot"),
  use_padj = TRUE,
  cross_type = c("deregulated", "upregulated", "downregulated"),
  cross_lfc_abs_lim = c(log2(1.5), 1),
  cross_min_signif = c(0.01, 0.05)
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{cross_id}}{Character vector; see \code{id} argument in
\code{\link{generate_a_list}} method}

\item{\code{use_padj}}{Whether to use the adjusted p-value or not}

\item{\code{cross_type}}{Vector of deregulation types; see \code{type} argument in
\code{\link{generate_a_list}} method
Available options are: "deregulated", "upregulated", "downregulated"
(default: c("deregulated", "upregulated", "downregulated"))}

\item{\code{cross_lfc_abs_lim}}{Vector of log fold change thresholds; see
\code{lfc_abs_lim} argument in \code{\link{generate_a_list}} method}

\item{\code{cross_min_signif}}{Vector of significance thresholds; see
\code{min_signif} argument in \code{\link{generate_a_list}} method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A nested data.frame containing crossed lists of genes for each
combination of filtering criteria
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-write_to_xlsx"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-write_to_xlsx}{}}}
\subsection{Method \code{write_to_xlsx()}}{
Write comparison results to xlsx
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$write_to_xlsx(
  file_suffix = "_pairwiseComp.xlsx",
  output_folder = ".",
  in_batch = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file_suffix}}{suffix for output file}

\item{\code{output_folder}}{output folder}

\item{\code{in_batch}}{Vector of batch codes; select or report only samples in these
batches}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
NULL
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-generate_summary"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-generate_summary}{}}}
\subsection{Method \code{generate_summary()}}{
Get the summary in number of deregulated genes from the comparisons
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$generate_summary(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{passed to \code{cross_args_and_generate_lists} method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
tibble with summary stats
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_summary"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_summary}{}}}
\subsection{Method \code{plot_summary()}}{
Plot the summary in number of deregulated genes from the comparisons
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_summary(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{passed to \code{generate_summary} method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
ggplot2 object with summary plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-extract_data_for_plot"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-extract_data_for_plot}{}}}
\subsection{Method \code{extract_data_for_plot()}}{
Return data for plots
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$extract_data_for_plot(
  use_padj = TRUE,
  lfc_abs_lim = 1,
  min_signif = 0.05,
  tag_id_select = "symbol",
  tag_id_show = "symbol",
  select_ids = NULL,
  select_batches = NULL,
  max_tags = 15,
  in_batches = NULL,
  hard_select = FALSE,
  safe_translate = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{use_padj}}{Whether to use the adjusted p-value or not}

\item{\code{lfc_abs_lim}}{Threshold for log fold change}

\item{\code{min_signif}}{Threshold for significance}

\item{\code{tag_id_select}}{Tag ID used for selection. Tag id can be "gid", "tgid",
"txid", "symbol" or "uniprot"}

\item{\code{tag_id_show}}{Tag ID used for plots. Tag id can be "gid", "tgid",
"txid", "symbol" or "uniprot"}

\item{\code{select_ids}}{Vector of IDs (\code{tag_id_select}) to select the genes}

\item{\code{select_batches}}{vector of batch to mark for selection}

\item{\code{max_tags}}{maximum number of tags to show}

\item{\code{in_batches}}{vector of batches code to keep (default all)}

\item{\code{hard_select}}{whether to remove unselected data}

\item{\code{safe_translate}}{If TRUE, the function will keep original ids
if the translation fails}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
long data frame for plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_ma"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_ma}{}}}
\subsection{Method \code{plot_ma()}}{
Return MA plot(s)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_ma(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{options passed to plot_de method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_vulcano"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_vulcano}{}}}
\subsection{Method \code{plot_vulcano()}}{
Return VULCANO plot(s)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_vulcano(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{options passed to plot_de method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_lfc_per_group"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_lfc_per_group}{}}}
\subsection{Method \code{plot_lfc_per_group()}}{
Return LFC per group plot(s)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_lfc_per_group(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{options passed to plot_de method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_lfc_per_group_facet_tags"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_lfc_per_group_facet_tags}{}}}
\subsection{Method \code{plot_lfc_per_group_facet_tags()}}{
Return LFC per group with tags as facet plot(s)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_lfc_per_group_facet_tags(geoms = c("bar", "errorbar"), ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{geoms}}{vector with geoms see method \code{plot_de}}

\item{\code{...}}{options passed to plot_de method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_de"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_de}{}}}
\subsection{Method \code{plot_de()}}{
Return plot(s)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_de(
  plot_type,
  lfc_limits = NULL,
  geoms = "point",
  show_selected_ids = TRUE,
  text_in_box = FALSE,
  max_nrow = 5,
  tag_ids_size = 4,
  tag_ids_alpha = 0.5,
  log2_expr = TRUE,
  batch_layout = "facet_grid_x",
  facet_scales = "fixed",
  facet_space = "fixed",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{plot_type}}{plot to draw either:
\itemize{
\item "ma" : MA-plot
\item "vulcano" : vulcano-plot
\item "lfc_per_group" : lfc per group
\item "lfc_per_group_facet_tags" : lfc per group with tag facet
}}

\item{\code{lfc_limits}}{values to draw hline and or vline for lfc}

\item{\code{geoms}}{vector with geoms to add:
\itemize{
\item point : compatible with all
\item bar : only for "lfc_per_group_facet_tags"
\item line : only for "lfc_per_group*"
\item errorbar : only for "lfc_per_group_facet_tags"
}}

\item{\code{show_selected_ids}}{whether or not to show the selected ids
(not compatible with \verb{lfc_per_group*})}

\item{\code{text_in_box}}{whether or not to use box to plot text
(show_selected_ids)}

\item{\code{max_nrow}}{maximum of facet row for the nested wrap facetting}

\item{\code{tag_ids_size}}{size when tag ids are plotted in plot area}

\item{\code{tag_ids_alpha}}{alpha when tag ids are plotted in plot area}

\item{\code{log2_expr}}{Whether to use log2(x+2) transform for expression data}

\item{\code{batch_layout}}{Define the type of layout for batch
\itemize{
\item "facet_nested_wrap"
\item "facet_grid_x"
\item "facet_grid_y"
\item "color_selected" represent all batch on the same plot with a focus on
a selected batch;only with plot_type "lfc_per_group_facet_tags"
with \code{color_selected} we must use \code{select_batches} from
\code{extract_data_for_plot}
}}

\item{\code{facet_scales}}{scales = "fixed" by default}

\item{\code{facet_space}}{space = "fixed" by default}

\item{\code{...}}{options passed to extract_data_for_plot method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-build_tag_hclust"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-build_tag_hclust}{}}}
\subsection{Method \code{build_tag_hclust()}}{
Return tag hclust
"NA" values are replaced with 0
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$build_tag_hclust(
  data = NULL,
  meth_dist = "minkowski",
  meth_value = "z",
  meth_clust = "centroid",
  max_abs_value = 100,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{Optional data similar to extract_data_for_plot output}

\item{\code{meth_dist}}{method to compute distance, see \code{dist()} function}

\item{\code{meth_value}}{value to compute distance
\itemize{
\item "z"
\item "log2FoldChange"
}}

\item{\code{meth_clust}}{method to compute clustering, see \code{hclust()} function}

\item{\code{max_abs_value}}{maximum absolute value to consider}

\item{\code{...}}{options passed to extract_data_for_plot method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-plot_heatmap"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-plot_heatmap}{}}}
\subsection{Method \code{plot_heatmap()}}{
Return HEATmap for a batch
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$plot_heatmap(
  meth_tag_dist = "minkowski",
  meth_tag_clust = "centroid",
  meth_tag_value = "z",
  max_abs_value_hclust = 100,
  plot_value = "z",
  max_abs_value_plot = 5,
  cut_tree_k = NULL,
  show_selected_ids = FALSE,
  meth_tag_dist2 = "minkowski",
  meth_tag_value2 = "z",
  meth_tag_clust2 = "ward.D2",
  max_abs_value_hclust2 = 100,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{meth_tag_dist}}{\code{meth_dist} for \code{build_tag_hclust()}}

\item{\code{meth_tag_clust}}{\code{meth_clust} for \code{build_tag_hclust()}}

\item{\code{meth_tag_value}}{\code{meth_value} for \code{build_tag_hclust()}}

\item{\code{max_abs_value_hclust}}{\code{max_abs_value} for \code{build_tag_hclust()}}

\item{\code{plot_value}}{value to plot}

\item{\code{max_abs_value_plot}}{maximum value for values to plot}

\item{\code{cut_tree_k}}{if provided with cut the heatmap into pieces based
on the hclust
\itemize{
\item "log2FoldChange"
\item "z" : z-score
}}

\item{\code{show_selected_ids}}{whether to show or not the selected ids.}

\item{\code{meth_tag_dist2}}{\code{meth_dist} for \code{build_tag_hclust()} when finding
gene cluster if cut_tree_k is not null}

\item{\code{meth_tag_value2}}{\code{meth_value} for \code{build_tag_hclust()} when finding
gene cluster if cut_tree_k is not null}

\item{\code{meth_tag_clust2}}{\code{meth_clust} for \code{build_tag_hclust()} when finding
gene cluster if cut_tree_k is not null}

\item{\code{max_abs_value_hclust2}}{\code{max_abs_value} for \code{build_tag_hclust()}
when finding gene cluster if cut_tree_k is not null}

\item{\code{...}}{options passed to extract_data_for_plot method}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The plot
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseComp-clone"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseComp-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseComp$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
