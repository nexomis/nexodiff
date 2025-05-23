% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_PairwiseDearseq.r
\name{PairwiseDearseq}
\alias{PairwiseDearseq}
\title{Dearseq Results from Pairwise Design}
\description{
A class representing DESeq2 results from a pairwise design.
}
\details{
This class generates results from expression data associated with a pairwise
design using the dearseq package.

This class inherits from the "PairwiseComp" class.
}
\section{Super class}{
\code{\link[nexodiff:PairwiseComp]{nexodiff::PairwiseComp}} -> \code{PairwiseDearseq}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-PairwiseDearseq-new}{\code{PairwiseDearseq$new()}}
\item \href{#method-PairwiseDearseq-clone}{\code{PairwiseDearseq$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="build_tag_hclust"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-build_tag_hclust'><code>nexodiff::PairwiseComp$build_tag_hclust()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="cross_args_and_generate_lists"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-cross_args_and_generate_lists'><code>nexodiff::PairwiseComp$cross_args_and_generate_lists()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="extract_data_for_plot"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-extract_data_for_plot'><code>nexodiff::PairwiseComp$extract_data_for_plot()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="filter_and_get_results"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-filter_and_get_results'><code>nexodiff::PairwiseComp$filter_and_get_results()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="generate_a_list"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-generate_a_list'><code>nexodiff::PairwiseComp$generate_a_list()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="generate_summary"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-generate_summary'><code>nexodiff::PairwiseComp$generate_summary()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="get_results"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-get_results'><code>nexodiff::PairwiseComp$get_results()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_de"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_de'><code>nexodiff::PairwiseComp$plot_de()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_heatmap"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_heatmap'><code>nexodiff::PairwiseComp$plot_heatmap()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_lfc_per_group"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_lfc_per_group'><code>nexodiff::PairwiseComp$plot_lfc_per_group()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_lfc_per_group_facet_tags"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_lfc_per_group_facet_tags'><code>nexodiff::PairwiseComp$plot_lfc_per_group_facet_tags()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_ma"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_ma'><code>nexodiff::PairwiseComp$plot_ma()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_summary"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_summary'><code>nexodiff::PairwiseComp$plot_summary()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="plot_vulcano"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-plot_vulcano'><code>nexodiff::PairwiseComp$plot_vulcano()</code></a></span></li>
<li><span class="pkg-link" data-pkg="nexodiff" data-topic="PairwiseComp" data-id="write_to_xlsx"><a href='../../nexodiff/html/PairwiseComp.html#method-PairwiseComp-write_to_xlsx'><code>nexodiff::PairwiseComp$write_to_xlsx()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseDearseq-new"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseDearseq-new}{}}}
\subsection{Method \code{new()}}{
Initialize \code{PairwiseDearseq} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseDearseq$new(
  expr_data,
  add_vars = c(),
  auto_filtering = TRUE,
  ncpus = 1,
  seed = 1234567,
  force_deseq2_norm = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{expr_data}}{ExprData : expression data}

\item{\code{add_vars}}{Character vactor with additional variable name in deisgn.
Those given variable will be add to the analysis formula. However the
fold-change in results will be the coefficient from the variable \code{group}}

\item{\code{auto_filtering}}{whether or not to perform independent filtering.}

\item{\code{ncpus}}{number of cpus to use for computation (default = 1)}

\item{\code{seed}}{seed for random computation in Dearseq}

\item{\code{force_deseq2_norm}}{renorm with deseq2 in place of internorm}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PairWiseCompDESeq2} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PairwiseDearseq-clone"></a>}}
\if{latex}{\out{\hypertarget{method-PairwiseDearseq-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PairwiseDearseq$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
