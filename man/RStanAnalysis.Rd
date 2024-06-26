% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_RstanAnalysis.r
\name{RStanAnalysis}
\alias{RStanAnalysis}
\title{Class runing the bayesian model for the distribution parameters inference}
\description{
This class generates results using rStan in order to inference diferent
distribution parameters from a PairwiseLFC object
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-RStanAnalysis-new}{\code{RStanAnalysis$new()}}
\item \href{#method-RStanAnalysis-get_results}{\code{RStanAnalysis$get_results()}}
\item \href{#method-RStanAnalysis-get_options}{\code{RStanAnalysis$get_options()}}
\item \href{#method-RStanAnalysis-get_formatted_results}{\code{RStanAnalysis$get_formatted_results()}}
\item \href{#method-RStanAnalysis-clone}{\code{RStanAnalysis$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-RStanAnalysis-new"></a>}}
\if{latex}{\out{\hypertarget{method-RStanAnalysis-new}{}}}
\subsection{Method \code{new()}}{
Initialize \code{RStanAnalysis} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RStanAnalysis$new(
  lfc_data,
  option_data = data.frame(distribution = rep(1, 2), base_alpha_equals_0 = c(1, 1),
    alpha_shift_equals_0 = c(1, 1), base_mu_equals_0 = c(1, 0), mu_shift_equals_0 = c(0,
    0), sigma_ratio_equals_1 = c(0, 0)),
  ncpus = 1,
  seed = 1234567,
  iter = 3000,
  refresh = 0,
  chains = 2,
  estimation_approach = c("sampling", "optimizing"),
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lfc_data}}{PairwiseLFC : log fold change daa}

\item{\code{option_data}}{is a dataframe with all the parametrization
of the models that must be run. Each row correpond to a model.
\itemize{
\item "distribution": 0 (gaussian), 1 (cauchy) or 2 (logistic);
default (1,1)
\item "base_alpha_equals_0": 0 or 1; default (1,1)
\item "alpha_shift_equals_0" 0 or 1; default (1,1)
\item "base_mu_equals_0": 0 or 1; default (1,0)
\item "mu_shift_equals_0": 0 or 1; default (0,0)
\item "sigma_ratio_equals_1": 0 or 1; default (0,0)
}}

\item{\code{ncpus}}{number of cpus to use for computation (default = 1)}

\item{\code{seed}}{seed for random computation in rStan (default=1234567)}

\item{\code{iter}}{number of iteration in the rStan sampling (default=3000)}

\item{\code{refresh}}{can be used to control how often the progress of the
sampling is reported (i.e. show the progress every refresh
iterations). By default, refresh = max(iter/10, 1). The progress
indicator is turned off if refresh <= 0 (default=0).}

\item{\code{chains}}{number of chains to run the rStan sampling (default=2)}

\item{\code{estimation_approach}}{approach to find one estimates.}

\item{\code{...}}{other parameters pass to rStan sampling (see the
documentation)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RStanAnalysis} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-RStanAnalysis-get_results"></a>}}
\if{latex}{\out{\hypertarget{method-RStanAnalysis-get_results}{}}}
\subsection{Method \code{get_results()}}{
Get the raw results (full fit object) from the bayesian
inference using rStan.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RStanAnalysis$get_results()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A nested list of lists: one list per batch that contains
one list per model. Each one of the model lists corresponds to the
fit object of rStan.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-RStanAnalysis-get_options"></a>}}
\if{latex}{\out{\hypertarget{method-RStanAnalysis-get_options}{}}}
\subsection{Method \code{get_options()}}{
Get the options used to parametrize the models.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RStanAnalysis$get_options()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A dataframe with all the parametrization
of the models used to crete the object. Each row correpond to a model
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-RStanAnalysis-get_formatted_results"></a>}}
\if{latex}{\out{\hypertarget{method-RStanAnalysis-get_formatted_results}{}}}
\subsection{Method \code{get_formatted_results()}}{
Get the formated data for the prediction
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RStanAnalysis$get_formatted_results(batchs = NULL, probs = c(0.25, 0.75))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{batchs}}{vector with the batchs names that should
be formatted. If NULL all the batchs are included in the
final output (default=NULL).}

\item{\code{probs}}{vector of the quantiles vector of
probabilities for which the quantiles should be return.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
list of dataframes (one per model in the options
data dataframe) with all the selected batchs as rows and
the the predictors (descriptors of the hyperparameters
of the models) as columns. The default descriptors are:
\itemize{
\item mean
\item sd
\item quantiles associeted whit the selected probabilities (\code{probs}).
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-RStanAnalysis-clone"></a>}}
\if{latex}{\out{\hypertarget{method-RStanAnalysis-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RStanAnalysis$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
