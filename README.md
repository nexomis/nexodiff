[[_TOC_]]

# 1. Installation

First clone the repository to `PATH` and then install on R:
```plantuml
@startuml
[*] --> PairwiseDesignWithAnnotation
PairwiseDesignWithAnnotation: * Design of analysis 
PairwiseDesignWithAnnotation: * Paths to expression files 
PairwiseDesignWithAnnotation: * Transcriptome annotation
note left
 **""design_filename""**
 **""annotation_filename""**
 **""idmapping_filename""**
end note
PairwiseDesignWithAnnotation -right-> ExprDataTranscript : **Import data at transcript level**
ExprDataTranscript : Expression data at transcript level
ExprDataTranscript --> ExprDataTranscript : **Filtering**
ExprDataTranscript --> ExprDataTranscript : **Normalization**
ExprDataTranscript --> ExprDataGene : **Summarization**
ExprDataGene : Expression data at gene level
ExprDataGene --> ExprDataGene : **Filtering**
ExprDataGene --> ExprDataGene : **Normalization**
ExprDataGene --> PairwiseLFC
PairwiseLFC: Pairwise Log2FC values
PairwiseLFC: * group VS control
PairwiseLFC: * control VS control
PairwiseLFC --> PairwiseDistBayMod
note right
 **""statistical_models""**
end note
PairwiseDistBayMod: Bayesian inference of parameters
PairwiseDistBayMod: based on statistical models
PairwiseDistBayMod: modelling LFC distributions
ExprDataGene -left-> PairwiseDESeq2
PairwiseDESeq2: Pairwise comparisons
PairwiseDESeq2: using DESeq2
PairwiseDESeq2 -left[dashed]-> PairwiseNestedGeneLists : **""$get_crossed_list()""**
PairwiseDesignWithAnnotation -down-> PairwiseNestedGeneLists
PairwiseNestedGeneLists: Nested gene list with design
PairwiseNestedGeneLists --> PairwiseNestedEnrichments
note left
 **""annotation_files""**
end note
PairwiseNestedEnrichments : Per list per functional database database
PairwiseNestedEnrichments : Gene set enrichments
PairwiseNestedEnrichments --> PairwiseNestedNetwork
PairwiseNestedNetwork: Per list per functional database database
PairwiseNestedNetwork: Network representation of enrichments
@enduml
```

## 1.1 Build documentation

```
devtools::build_vignettes(pkg = PATH)
roxygen2::roxygenize(package.dir = PATH)
```

## 1.2 Build and install the library

```
devtools::install_local(path=PATH,force=T,upgrade = 'never')
```

# Getting started 

See the vignette `Examples`


# check dependencies 

```
cat R/* | grep "::" | sed -re "s/.*[\(\[ ]([a-zA-Z0-9_]+)\:\:.*/\1/g" | sort -u | xargs | sed 's/ /, /g'
```