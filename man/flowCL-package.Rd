\name{flowCL-package}
\alias{flowCL-package}
\docType{package}
\title{
    Semantic labelling of flow cytometric cell populations.
}

\description{
    flowCL uses SPARQL to match a phenotype to a cell type from the cell ontology. If the match is not unique, then the best alternative is returned.
}
\details{
\tabular{ll}{
Package: \tab flowCL\cr
Type: \tab Package\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
Depends: \tab Rgraphviz, SPARQL\cr
}

flowCL executes queries against the Cell Ontology (CL), available at http://cellontology.org. The CL file is hosted on a triplestore, i.e., a database for storage and retrieval of Resource Description Framework (RDF) triples. The SPARQL endpoint at http://cell.inference.me:7200/repositories/CL (old endpoint: http://cell.ctde.net:8080/openrdf-sesame/repositories/CL) is used to execute the SPARQL queries retrieving the correct matches from the CL. While other SPARQL endpoints can be used, users should be aware that in our case the CL file has been reasoned upon, and resulting extra inferred axioms have been added to the triplestore, providing a more complete result set.

}

\author{
    Maintainer: Justin Meskas <justinmeskas@gmail.com>

    Authors: Justin Meskas, Radina Droumeva
}

\keyword{ FlowCytData }

\examples{

# # Load a pre-loaded archive. Skipping this chuck will cause flowCL to
# ## slowly build a new one.
# flowCL("archive")
# 
# # Simple two marker example
# Res <- flowCL("CCR7+CD45RA+")
# tmp <- Res$'CCR7+CD45RA+'
# plot(tmp[[1]], nodeAttrs=tmp[[2]], edgeAttrs=tmp[[3]], attrs=tmp[[4]])
# Res$Table
# 
# # Exact match example
# Res <- flowCL("CD3-CD19-CD20-CD14+", Verbose = TRUE, OntolNamesTD = TRUE)
# tmp <- Res$'CD3-CD19-CD20-CD14+'
# plot(tmp[[1]], nodeAttrs=tmp[[2]], edgeAttrs=tmp[[3]], attrs=tmp[[4]])
# Res$Table
# 
# # Cell Label Example
# x <-"CD3-CD19-CD20-CD14+"
# Res <- flowCL(x)
# Res$Cell_Label[[x]][[1]]
# 
# # As a secondary way to view the results,
# ## see "[current directory]/flowCL_results/".
# # Figures created called tree_(phenotype).pdf give the cell hierarchy
# ## dependent on the markers in the phenotype.
# # A list of results from Res$Table are stored in listPhenotypes.csv.
}
