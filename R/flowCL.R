#' @importFrom magrittr %>% %<>%
#' @importFrom grDevices dev.off pdf
#' @importFrom methods new
#' @importFrom utils data read.csv write.csv write.table

#' @export
flowCL <- function (
    MarkerList = "HIPC",
    ExpMrkrLst = NULL,
    Indices = NULL,
    Verbose = FALSE,
    KeepArch = TRUE,
    MaxHitsPht = 5,
    OntolNamesTD = FALSE,
    ResetArch = FALSE,
    VisualSkip = FALSE
) {

# message <- "The endpoint no longer exists. \nPlease email Justin at justinmeskas@gmail.com \nto let him know that you rely on flowCL \nand that you would like it working again. \nThank you.\n"
# cat(message)
# return()

# flowCL (Semantic labelling of flow cytometric cell populations)
# Authors: Justin Meskas (justinmeskas@gmail.com), Radina Droumeva (radina.droumeva@gmail.com )

initialTime <- Sys.time ( )

VisualizationSkip <- VisualSkip       # Skips the making of the visualization (tree diagram)
OntolNamesTD      <- OntolNamesTD     # Shows either the common marker name or the ontology marker name in the tree diagram

finalResults <- list ( )

# Load some required data from supportingFunctions.
listPhenotypes_flowCL <- listPhenotypes_flowCL()
listColours_flowCL    <- listColours_flowCL()

# Load querying data from supportingFunctions.
prefix.info                           <- flowCL_query_data_prefix.info()
que.getParentClasses                  <- flowCL_query_data_getParentClasses()
que.hasPlasmaMembranePart             <- flowCL_query_data_hasPlasmaMembranePart()
que.hasProperLabel                    <- flowCL_query_data_hasProperLabel()
# que.hasPMPsingle                      <- flowCL_query_data_hasPMPsingle()
que.hasProperSynonym                  <- flowCL_query_data_hasProperSynonym()
# que.lacksPMPsingle                    <- flowCL_query_data_lacksPMPsingle()
que.lacksPlasmaMembranePart           <- flowCL_query_data_lacksPlasmaMembranePart()
que.hasLowPlasmaMembraneAmount        <- flowCL_query_data_hasLowPlasmaMembraneAmount()
que.hasHighPlasmaMembraneAmount       <- flowCL_query_data_hasHighPlasmaMembraneAmount()
que.celllabel_lacksPMP                <- flowCL_query_data_celllabel_lacksPMP()
que.celllabel_hasPMP                  <- flowCL_query_data_celllabel_hasPMP()
que.celllabel_lowPMA                  <- flowCL_query_data_celllabel_lowPMA()
que.celllabel_highPMA                 <- flowCL_query_data_celllabel_highPMA()

# Function to test if the user has input cd in non-capitals. If so, the code will terminate.
if ( cdTest ( MarkerList ) == TRUE)
    return()

# Define the cell.ctde.net SPARQL endpoint
# endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL" # original
# endpoint <- "http://75.127.15.173:8080/openrdf-sesame/repositories/CL" # Jonathan April 18 2017 (different way to get to Alan's)
# endpoint <- "http://cell.inference.me:7200/repositories/CL" # Jonathan Nov 17 2017
endpoint <- "http://Tchernobog:7200/repositories/flowCL" # JJJ test 8 Sept 2021
print(endpoint)


# Check date of Ontology update
if ( length(MarkerList) == 1 ) {
    if ( tolower(MarkerList) == "date") {
        que.getDate <- flowCL_query_date()
        res <- SPARQL ( url = endpoint, paste ( que.getDate, collapse = "\n" ) )$results
        return ( c ( paste ( strsplit(res[1,2], split="/" )[[1]][7], ": year-month-day" ) ) )
    }
}

# Check to create archive or not
if ( length(MarkerList) == 1 ) {
    if ( tolower(MarkerList) == "archive") {

        flowCL_archive <- NULL; remove(flowCL_archive) # load and remove so 'check' doesnt give any notes/warnings
        data(flowCL_archive) # load data
        # seperate data into different variables for understanding purposes
        Parents_query_archive <- flowCL_archive[[1]]; hasPMP_archive  <- flowCL_archive[[3]];  lacksPMP_archive <- flowCL_archive[[5]];  highPMA_archive <- flowCL_archive[[7]];
        lowPMA_archive        <- flowCL_archive[[9]]; results_archive <- flowCL_archive[[11]];  Parents_Names    <- flowCL_archive[[2]];  hasPMP_Names    <- flowCL_archive[[4]];
        lacksPMP_Names        <- flowCL_archive[[6]]; highPMA_Names   <- flowCL_archive[[8]]; lowPMA_Names     <- flowCL_archive[[10]]; results_Names   <- flowCL_archive[[12]];
        FolderNames <- c("/parents_query/", "/reverse_query/hasPMP/", "/reverse_query/lacksPMP/", "/reverse_query/highPMA/", "/reverse_query/lowPMA/", "/results/")
        tempArchive <- list(Parents_query_archive, hasPMP_archive, lacksPMP_archive, highPMA_archive, lowPMA_archive, results_archive)
        tempNames   <- list(Parents_Names,         hasPMP_Names,   lacksPMP_Names,   highPMA_Names,   lowPMA_Names,   results_Names)
        for ( k1 in 1 : 6 ) {
            tempPath <- paste(getwd(), "/flowCL_results", FolderNames[k1], sep="")
            dir.create ( tempPath, showWarnings=FALSE, recursive=TRUE )
            for ( k2 in 1 : length(tempNames[[k1]] ) ) {
                # create each of the .csv files
                write.table(tempArchive[[k1]][[k2]],paste(getwd(),"/flowCL_results", FolderNames[k1], tempNames[[k1]][[k2]], sep=""), sep=",", row.names = FALSE)
            }
        }
        return ( cat ( "Archive created in \"", paste ( getwd ( ) , "/flowCL_results", sep=""), "\"\n", sep="" ) )
    }
}

# Load phenotypes.
suppressWarnings ( if ( MarkerList == "HIPC" ) {
    listPhenotypes <- listPhenotypes_flowCL
} else {
    listPhenotypes <- MarkerList
} )
listPhenotypes <- as.matrix ( listPhenotypes )

# Calculate the start and end point of the iterations.
if ( is.null ( Indices ) ) {
    IterStart = 1
    IterEnd   = length ( listPhenotypes )
    markersToQuery <- IterStart:IterEnd
} else {
    listPhenotypes <- listPhenotypes[Indices]
    markersToQuery <- 1:length(Indices)
}

# Send a warning that ExpMrkrLst was not defined.
if ( is.null ( ExpMrkrLst ) ) {
    cat("ExpMrkrLst was not defined. Defining ExpMrkrLst as MarkerList.\n")
    listExpPhenotypes <- listPhenotypes # updated MarkerList
} else {
    # If the user used a different number of phenotypes as experimental phenotypes then a warning is shown.
    if ( length(ExpMrkrLst) != length(listPhenotypes) ) {
        if ( length(ExpMrkrLst) == 1 ) {
            message("ExpMrkrLst is of length 1. ExpMrkrLst is assumed to be the experimental marker list for all phenotypes.")
            listExpPhenotypes <- rep(ExpMrkrLst,length(MarkerList))
        } else {
            message("ExpMrkrLst is not of the same length as MarkerList nor of length 1.")
        }
    } else { # if there are as many phenotypes as experimental phenotypes then do nothing.
        listExpPhenotypes <- ExpMrkrLst
    }
}

# Preallocate lists.
listPhenotypeUpdate  <- listPhenotypeOriginal <- listPhenotypeID <- listCellID   <- listMarkers <- listRanking <- listPhenotypes
listPhenotypesuccess <- listMarkerLabels      <- listCellLabels  <- listFullMrkr <- listExpPhenotypeUpdate     <- listPhenotypes

# Default listPhenotypesuccess as a "No"
for ( q in markersToQuery ) {
    listPhenotypesuccess[[q]] <-   "No"
}

# Removes the flowCL_results directory if ResetArch is TRUE
if ( ResetArch == TRUE ) {
    unlink ( "flowCL_results", recursive = TRUE )
}

# Define directories for storing results.
save.dir             <- 'flowCL_results/'
save.dirResults      <- 'flowCL_results/results/'
save.dirParents      <- 'flowCL_results/parents/'
save.dirParentsQuery <- 'flowCL_results/parents_query/'
save.dirlacksPMP     <- 'flowCL_results/reverse_query/lacksPMP/'
save.dirhasPMP       <- 'flowCL_results/reverse_query/hasPMP/'
save.dirlowPMA       <- 'flowCL_results/reverse_query/lowPMA/'
save.dirhighPMA      <- 'flowCL_results/reverse_query/highPMA/'

# Create directories.
dir.create ( save.dir,             showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirResults,      showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirParents,      showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirParentsQuery, showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirlacksPMP,     showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirhasPMP,       showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirlowPMA,       showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirhighPMA,      showWarnings=FALSE, recursive=TRUE )


# Create a file for future quick matching of the short name with the ontology name of each marker.
fname <- paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
if( !file.exists ( fname ) ) {
    temp.table <- matrix (0, 2, 2)
#     temp.table[1] <- paste ( 'CD4' );    temp.table[2] <- paste ( 'CD4 molecule' )
    temp.table[1] <- paste ( 'CD3e' );    temp.table[2] <- paste ( 'CD3 epsilon' )
    temp.table[3] <- paste ( 'CD3z' );    temp.table[4] <- paste ( 'T cell receptor complex' )
    write.table ( t ( temp.table ), fname, sep = ",", col.names = FALSE, row.names = FALSE )
}

# Create null variables for memory allocation and easy additions of indices
rEG <- nAttrs <- eAttrs <- attrs <- list ( )

#################################################################################
# Main body starts here. This for-loop goes through each phenotype one at a time.
for ( q in markersToQuery ) {

    # Get starting time to keep track of each phenotype's running time.
    start <- Sys.time ( )

    #----------------------------------------------------------- Update phenotype names

    # Define a phenotype to test here.
    phenotype <- listPhenotypes[[q]]

    # Change the phenotype to have only hi and lo and not other variations
    phenotype <- gsub("bright", "hi", phenotype)
    phenotype <- gsub("bri",    "hi", phenotype)
    phenotype <- gsub("high",   "hi", phenotype)
    phenotype <- gsub("[+][+]", "hi", phenotype)
    phenotype <- gsub("dim",    "lo", phenotype)
    phenotype <- gsub("low",    "lo", phenotype)
    phenotype <- gsub("--",     "lo", phenotype)

    # Output for the user
    if ( Verbose == TRUE ) { cat("\nThe phenotype of interest is", phenotype, "\n" ) }

    phenotype <- gsub("CD3-", "CD3e-", phenotype)
    phenotype <- gsub("HLA-DR", "HLA.DR", phenotype)

    # The following breaks up the phenotype into single markers and sorts them
    # by their expression (only positive or negative expression implemented for now)
    # see supportingFunctions.R for details on how 'phenoParse' works.
    marker.list <- phenoParse ( phenotype )
    if ( is.null ( marker.list ) ) {
        return ()
    }
    # Change the . in HLA.DR to a - since the + and - signs are reserved for splitting the string.
    marker.list <- changeHLADR ( marker.list )
    # Creates another copy which will have the + and - signs. Used by treeDiagram and when searching files.
    marker.list.short <- marker.list

    # Put the +, -, hi and lo signs back on marker.list.short. Used by treeDiagram and when searching files.
    AddOn <- c("+", "-", "hi", "lo")
    for ( q3 in 1:length(marker.list)){
        if ( length ( marker.list.short[[q3]] ) != 0 ) {
            for ( q2 in 1:( length ( marker.list.short[[q3]] ) ) ) {
                marker.list.short[[q3]][q2] <- paste ( marker.list.short[[q3]][q2], AddOn[q3], sep = "" )
            }
        }
    }

    # Update the marker list with the full label names from the ontology.
    marker.list <- ontologyLabel ( marker.list, marker.list.short, Verbose = Verbose, save.dir = save.dir,
                                    que.hasProperLabel = que.hasProperLabel, que.hasProperSynonym = que.hasProperSynonym, prefix.info = prefix.info, endpoint = endpoint)

    # Make a list of the ontology names for each phenotype searched for, which will be exported to a table in .csv form

    listPhenotypeUpdate[[q]] <- phenoUnparse ( phenotype, marker.list )


    phenotype <- gsub("CD3e-", "CD3-", phenotype)
    phenotype <- gsub("HLA.DR", "HLA-DR", phenotype)
    marker.list.short[[2]][which (marker.list.short[[2]]== "CD3e-")] <- "CD3-"

    #----------------------------------------------------------- Update experiment phenotype names

    # Define a phenotype to test here.
    exp.phenotype <- listExpPhenotypes[[q]]

    exp.phenotype <- gsub("CD3-", "CD3e-", exp.phenotype)
    exp.phenotype <- gsub("HLA-DR", "HLA.DR", exp.phenotype)
    # Change the phenotype to have only hi and lo and not other variations
    exp.phenotype <- gsub("bright", ",", exp.phenotype)
    exp.phenotype <- gsub("bri",    ",", exp.phenotype)
    exp.phenotype <- gsub("high",   ",", exp.phenotype)
    exp.phenotype <- gsub("hi",     ",", exp.phenotype)
    exp.phenotype <- gsub("[+][+]", ",", exp.phenotype)
    exp.phenotype <- gsub("[+]",    ",", exp.phenotype)
    exp.phenotype <- gsub("dim",    ",", exp.phenotype)
    exp.phenotype <- gsub("low",    ",", exp.phenotype)
    exp.phenotype <- gsub("lo",     ",", exp.phenotype)
    exp.phenotype <- gsub("--",     ",", exp.phenotype)
    exp.phenotype <- gsub("-",      ",", exp.phenotype)
    exp.phenotype <- gsub(",$",      "", exp.phenotype)
    exp.phenotype <- gsub("HLA[.]DR", "HLA-DR", exp.phenotype)

    listExpPhenotypes[[q]] <- gsub("CD3e", "CD3", exp.phenotype)

    # Output for the user
    if ( Verbose == TRUE ) { cat("The experiment phenotype of interest is", gsub("CD3e", "CD3", exp.phenotype), "\n" ) }

    exp.marker.list <- strsplit ( x = exp.phenotype, split = "," )

    # Change the . in HLA.DR to a - since the + and - signs are reserved for splitting the string.
    exp.marker.list <- changeHLADR ( exp.marker.list )
    # Creates another copy which will have the + and - signs. Used by treeDiagram and when searching files.
    exp.marker.list.short <- exp.marker.list

    # Update the marker list with the full label names from the ontology.
    exp.marker.list <- ontologyLabel ( exp.marker.list, exp.marker.list.short, Verbose = Verbose, save.dir = save.dir,
                                    que.hasProperLabel = que.hasProperLabel, que.hasProperSynonym = que.hasProperSynonym, prefix.info = prefix.info, endpoint = endpoint )
    # Make a list of the ontology names for each phenotype searched for, which will be exported to a table in .csv form
    listExpPhenotypeUpdate[[q]]   <- paste ( unlist(exp.marker.list), collapse="\n")

    # remove CD3e- and replace with CD3-
    exp.phenotype <- gsub("CD3e", "CD3", exp.phenotype)
    exp.marker.list.short[[1]][which (exp.marker.list.short[[1]]== "CD3e")] <- "CD3"

    #----------------------------------------------------------- Query and save data

    # Test if the markers are a subset of the experimental markers
    returnTrue <- FALSE
    for ( j1 in 1 : length ( unlist ( marker.list ) ) ) {
        if ( length ( which ( unlist ( marker.list ) [j1] == unlist ( exp.marker.list ) ) ) < 1 ) {
            message(paste("The marker ", unlist ( marker.list.short ) [j1]," (", unlist ( marker.list ) [j1], ") in MarkerList needs to also be in ExpMrkrLst", sep=""))
            returnTrue <- TRUE
        }
    }
    if (returnTrue == TRUE)
        return()

    # Default skipQuery to TRUE. If it changes to FALSE then querying will have to be done.
    skipQuery <- TRUE
    # Create res as a NULL variable for easy additions of indices
    res <- NULL
    # Cycle through has, lacks, low and high markers.
    # For has, use a SPARQL query using the property 'has plasma membrane part' of each marker.
    # For the lacks markers, query 'lacks plasma membrane part'.
    # For the low markers, query 'has low plasma membrane amount'.
    # For the high markers, query 'has high plasma membrane amount'.
    for ( i in unlist ( marker.list.short ) ) {
        fname <- paste ( save.dirResults, 'results_', i, '.csv', sep = "" )
        if ( file.exists ( fname ) ) {
            tempCSV <- read.csv ( fname, as.is=TRUE )
            if ( length ( which ( tempCSV[,"Number.Of.Hits"] > 1 )) >= 1 ) { # Number Of Hits column
                warning ( "You are receiving more marker labels than markers that were input. Most likely the input markers are not defined correctly." )
            }
            # delete the last two columns to combine and redue the results of the last two columns with other markers
            tempCSV <- tempCSV[ - c ( ncol ( tempCSV ) , ncol ( tempCSV ) - 1 ) ]
            res <- rbind ( res, tempCSV )
        } else {
            if ( Verbose == TRUE ) { cat ( "At least one marker was not previously queried. Querying all.\n" ) }
            skipQuery <- FALSE
            break;
        }
    }

    # If all results files exist then no querying needs to be done
    if ( skipQuery == TRUE ) {

        clean.res <- tabulateResults ( res )

    } else {
        # Initialize result collector.
        res <- NULL
        for ( q3 in 1:length(marker.list)){
            for ( m in marker.list[[q3]] ) {
                if ( Verbose == TRUE ) { cat( "Locating marker", m, "\n" ) }
                # Get relevant information about the marker - the population names which
                # have "plasma membrane part" of this marker.
                fname <- paste ( save.dirResults, 'results_', marker.list.short[[q3]][which(marker.list[[q3]]==m)], '.csv', sep = "" )
                QueryFile <- list(que.hasPlasmaMembranePart,       que.lacksPlasmaMembranePart,       que.hasHighPlasmaMembraneAmount,       que.hasLowPlasmaMembraneAmount)
                if ( file.exists ( fname ) ) {
                    cur.res <- read.csv ( fname, as.is=TRUE )
                    cur.res <- cur.res[ - c ( ( ncol ( cur.res ) - 2 ) : ncol ( cur.res ) ) ]
                } else {
                    if ( q3 == 1 ) { cur.res <- queryMarker ( marker = m, query.file = QueryFile[[q3]], prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=TRUE ) }
                    if ( q3 == 2 ) { cur.res <- queryMarker ( marker = m, query.file = QueryFile[[q3]], prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=TRUE ) }
                    if ( q3 == 3 ) { cur.res <- queryMarker ( marker = m, query.file = QueryFile[[q3]], prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=TRUE ) }
                    if ( q3 == 4 ) { cur.res <- queryMarker ( marker = m, query.file = QueryFile[[q3]], prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=TRUE ) }
                    # temp fix for "\"has plasma membrane part\"@en" bug
                    # cur.res[,3] <- sapply(1:length(cur.res[,3]), function(x) {substr(cur.res[x,3], 2, nchar(cur.res[x,3])-4)}) # " and "@en removal

                    temp.loc <- which(colnames(cur.res) == "pl")
                    if (length(temp.loc) >= 1){
                        colnames(cur.res)[which(colnames(cur.res) == "pl")] <- "plabel"
                    }

                }
                if ( nrow ( cur.res ) == 0 ) {
                    if ( Verbose == TRUE ) {
                    if ( q3 == 1 ) {message( "No cell populations found which have plasma membrane part ", m, "!" )}
                    if ( q3 == 2 ) {message( "No cell populations found which lack plasma membrane part ", m, "!" )}
                    if ( q3 == 3 ) {message( "No cell populations found which have high plasma membrane amount ", m, "!" )}
                    if ( q3 == 4 ) {message( "No cell populations found which have low plasma membrane amount ", m, "!" )}
                    }

                }

                penalties <- rep ( 0, nrow ( cur.res ) )
                cur.res <- cbind ( cur.res, penalties )
                temp.clean.res <- tabulateResults ( cur.res )
                temp.name <- marker.list.short[[q3]][which ( m == marker.list[[q3]] )]
                fname <- paste ( save.dirResults, 'results_', temp.name, '.csv', sep = "" )
                # save the marker's data to skip the query next time
                write.csv ( temp.clean.res, fname, row.names = FALSE )
                res <- rbind ( res, cur.res )
                if ( length ( which ( temp.clean.res[,"Number Of Hits"]  > 1 )) >= 1 ) { # Number Of Hits column
                    warning ( "You are receiving more marker labels than markers that were input. Most likely the input markers are not defined correctly." )
                }
            }
        }

        # For each returned owl object ID, tabulate how many times the result was
        # returned. This is essentially the "hits" part of the score -- telling us
        # how many of the markers we have were matched to each population. From this,
        # the penalty tally will be subtracted to get a final overall score for each
        # population label.
        clean.res <- tabulateResults ( res )

    }

    # Save the results.
    fname <- paste ( save.dirResults, 'results_', phenotype, '.csv', sep = "" )
    write.csv ( clean.res, fname, row.names = FALSE ) # this is only redundant if the query is for only one marker
    if ( Verbose == TRUE ) { cat( 'Initial query results saved in\n  ', paste("[current directory]/", fname, sep=""), "\n" ) }
#     return(clean.res)
    # If there were no hits the code will skip to the next iteration instead of
    # proceeding and getting an error
    if ( nrow ( clean.res ) == 0 ) {
        if ( Verbose == TRUE ) {
            cat ( "Time elapsed:", timeOutput(start), "\n" )
            cat ( "Iterations at", which ( markersToQuery == q ) , "out of", length ( markersToQuery ), "\n" )
        }
        next
    }

    #----------------------------------------------------------- Refine clean.res
    # Refine list of parents by focusing on highest scores:
    scores <- as.numeric ( as.character ( clean.res[ , "Score"] ) )

    # If at least 3 perfect scores are found, only report them:
    if ( length ( which ( scores == length ( unlist ( marker.list ) ) ) ) >= 3 ) {
        cutoff.score <- length ( unlist ( marker.list ) )
    } else {
        # Typically there will be less than perfect scores (due to marker missing/
        # marker not declared in the population/non-typical phenotype)

        # Only show the matches with the most hits. If there are only a few that are off by 1 then show them as well.
        if ( length(unlist(marker.list)) <= 2 ) {
            cutoff.score <- as.integer ( clean.res[1,"Number Of Hits"] )
        } else {
            cutoff.score <- as.integer ( clean.res[1,"Number Of Hits"] ) - 1
        }
    }

    if ( Verbose == TRUE ) { cat ( "Using a cutoff score of", cutoff.score, "\n" ) }

    # Extract highly scored parents only:
    clean.res <- clean.res[which ( scores >= cutoff.score ), ]

        # Fixes small bug when there is only one result
    if ( length ( which ( scores >= cutoff.score ) ) == 1 ) {
        clean.res <- as.matrix ( clean.res )
        clean.res <- t ( clean.res )
    }

    #----------------------------------------------------------- Cell Label
    # Compile the new row of the .csv file by updating all the list functions.

    r <- updateLists ( clean.res=clean.res, length( which ( clean.res[ ,'Number Of Hits'] >= ( max ( as.numeric ( clean.res[ ,'Number Of Hits'] ) ) - 1 ) )), FALSE)
    listMarkerLabels[[q]] <- r[1] ; listCellLabels[[q]] <- r[2]
    listPhenotypeID[[q]]  <- r[3] ; listCellID[[q]]     <- r[4]

    temp.celllabel <- listCellLabels[[q]]
    temp.celllabel <- gsub ( "\n [+] more", "", temp.celllabel )
    temp.string <- list()
    temp.index <- gregexpr ( ")",  temp.celllabel )
    # Extract cell labels from text
    if ( ( temp.index[[1]][1] != - 1 ) ) {
        for ( q7 in 1:length ( temp.index[[1]] ) ) {
            if ( q7 < length ( temp.index[[1]] ) ) {
                if (q7 <= 8 )            {temp.string[q7] <- paste (substr ( temp.celllabel, temp.index[[1]][q7] + 2, temp.index[[1]][q7+1] - 3 ), sep = "" )}
                if (q7 >= 9 && q7 <= 98) {temp.string[q7] <- paste (substr ( temp.celllabel, temp.index[[1]][q7] + 2, temp.index[[1]][q7+1] - 4 ), sep = "" )}
                if (q7 >= 99 )           {temp.string[q7] <- paste (substr ( temp.celllabel, temp.index[[1]][q7] + 2, temp.index[[1]][q7+1] - 5 ), sep = "" )}
            }
            if ( q7 == length ( temp.index[[1]] ) ) {
                temp.string[q7] <- paste (substr ( temp.celllabel, temp.index[[1]][q7] + 2, nchar(temp.celllabel)), sep = "" )
            }
        }
    }
    # Store cell labels to return to the user
    finalResults[["Cell_Labels"]][[listPhenotypes[q]]] <- temp.string

    #----------------------------------------------------------- Marker_Groups and Markers

    postfix <- c("_lacksPMP.csv", "_hasPMP.csv", "_lowPMA.csv", "_highPMA.csv")
    query.file.list <- list(que.celllabel_lacksPMP, que.celllabel_hasPMP, que.celllabel_lowPMA, que.celllabel_highPMA)
    query.dir.list  <- c(save.dirlacksPMP, save.dirhasPMP, save.dirlowPMA, save.dirhighPMA)

    # Organize which markers are required, which are extra, which are in the experiment and not used,
    # and which additional ones would be required to get a perfect match.

    MarkerGroups <- MarkerGroupsFunc(
        temp.string=temp.string,
        query.dir.list=query.dir.list,
        postfix=postfix,
        save.dir=save.dir,
        query.file.list=query.file.list,
        prefix.info=prefix.info,
        Verbose=Verbose,
        marker.list.short=marker.list.short,
        exp.marker.list.short=exp.marker.list.short,
        AllMarkerGroups=FALSE,
        endpoint=endpoint
        )

    finalResults[["Marker_Groups"]][[listPhenotypes[q]]] <- MarkerGroups # List form
    # Display in a bracket form
    temp.list <- list()
    for ( q5 in 1 : length ( MarkerGroups ) ) {
        temp.list[[q5]] <- paste (
            "{", paste(MarkerGroups[[q5]][[1]], collapse=", "),
            "}", paste(MarkerGroups[[q5]][[2]], collapse=", "),
            "(", paste(MarkerGroups[[q5]][[3]], collapse=", "), ")",
            "[", paste(MarkerGroups[[q5]][[4]], collapse=", "), "]"
            )

        finalResults[["Markers"]][[listPhenotypes[q]]] <- temp.list # Bracket form
    }

    # Rank the results
    Ranking <- NULL
    KeepMarkerNoDouble <- rep(TRUE, length(MarkerGroups))
    for ( r1 in 1 : length(MarkerGroups)){ # remove + - hi lo
        tmp1 <- MarkerGroups[[r1]][[1]]
        tmp1 <- gsub("[+]$",  "", tmp1); tmp1 <- gsub("-$",  "", tmp1); tmp1 <- gsub("hi$", "", tmp1);   tmp1 <- gsub("lo$", "", tmp1)
        tmp2 <- MarkerGroups[[r1]][[3]]
        tmp2 <- gsub("[+]$",  "", tmp2); tmp2 <- gsub("-$",  "", tmp2); tmp2 <- gsub("hi$", "", tmp2);   tmp2 <- gsub("lo$", "", tmp2)

        tmp3 <- MarkerGroups[[r1]][[1]]
        tmp3 <- gsub("hi$", "+", tmp3);   tmp3 <- gsub("lo$", "+", tmp3)
        tmp4 <- MarkerGroups[[r1]][[3]]
        tmp4 <- gsub("hi$", "+", tmp4);   tmp4 <- gsub("lo$", "+", tmp4)

        DoubleMisUseMarker <- 0 # Calculates number of markers that were queried and matched but were queried with the wrong tag
        for ( w1 in 1 : length(tmp2) ) {
            if( length(which (tmp1 == tmp2[w1])) >= 1 ) { # equal without tag
                if( length(which (tmp3 == tmp4[w1])) == 0){ # not equal with tag
                    DoubleMisUseMarker <- DoubleMisUseMarker + 1
                    KeepMarkerNoDouble[r1] <- FALSE
                }
            }
        }

        PenaltiesOfMarkers <- length(MarkerGroups[[r1]][[1]]) # Calculates number of markers that were queried that are not part of the definition of the cell type
        TotalNumberOfMarkers <- length(MarkerGroups[[r1]][[2]]) + length(MarkerGroups[[r1]][[3]]) + length(MarkerGroups[[r1]][[4]]) # Number of markers to define the cell type
        # Scoring system: [ +1 for each hit - (4 * markers that were required but were given the wrong tag) - 2 * markers that were queried that were not required ] / all required to define the cell type
        Ranking[r1] <-  (length(MarkerGroups[[r1]][[2]]) - 4*DoubleMisUseMarker- 2*(PenaltiesOfMarkers - DoubleMisUseMarker)) / (TotalNumberOfMarkers)
    }

    if ( is.null(finalResults[["Ranking"]][[listPhenotypes[q]]]) && length(Ranking) == 1) { # fixes bug where only one numeric in the first slot would create a list of single numerics
        finalResults[["Ranking"]][[listPhenotypes[q]]] <- as.list(Ranking)
        finalResults[["Ranking"]][[listPhenotypes[q]]] <- finalResults[["Ranking"]][[listPhenotypes[q]]][[1]] # sloppy trick to set list up correctly
    } else {
        finalResults[["Ranking"]][[listPhenotypes[q]]] <- Ranking
    }

    #----------------------------------------------------------- Reorder Results

    SortOrder <- sort(finalResults[["Ranking"]][[listPhenotypes[q]]], decreasing=TRUE, index.return = TRUE )$ix

    # fixes a bug caused by removing all cell types because they all are
    if ( all ( KeepMarkerNoDouble == FALSE ) ) {
        warning ( paste0("All cell types of ", phenotype, " are using at least one marker with the wrong tag. ",
                    "These cases are deemed to be completely wrong and should be removed, however, ",
                    "they are the only ones left so flowCL will not remove them to avoid an error. ",
                    "Please check to see that you inputted correctly, and consider quierying something different.", collapse = "") )
        KeepMarkerNoDouble <- rep(TRUE, length(MarkerGroups))
    }

    tmpTrue <- which(KeepMarkerNoDouble == FALSE)

    if ( length(tmpTrue) >= 1 ) {
        for ( m1 in 1 : length(tmpTrue) ) {
            SortOrder <- SortOrder[-which ( SortOrder == tmpTrue[m1])]
        }
    }

    # Remembers if we need to add a "+ more" comment provided that we removed at least one cell type
    AddPlusMore <- FALSE
    if ( length(SortOrder) > MaxHitsPht ) {
        AddPlusMore <- TRUE # track if "+ more" needs to be added
        SortOrder <- SortOrder[1:MaxHitsPht] # only takes the top MaxHitsPht amount
    }

    finalResults[["Ranking"]][[listPhenotypes[q]]]        <- finalResults[["Ranking"]][[listPhenotypes[q]]][SortOrder]
    finalResults[["Markers"]][[listPhenotypes[q]]]        <- finalResults[["Markers"]][[listPhenotypes[q]]][SortOrder]
    finalResults[["Marker_Groups"]][[listPhenotypes[q]]]  <- finalResults[["Marker_Groups"]][[listPhenotypes[q]]][SortOrder]
    finalResults[["Cell_Labels"]][[listPhenotypes[q]]]     <- finalResults[["Cell_Labels"]][[listPhenotypes[q]]][SortOrder]

    if ( length(which(finalResults[["Ranking"]][[listPhenotypes[q]]] == 1 )) >= 1 ) {
        PerfectIndices <- which(finalResults[["Ranking"]][[listPhenotypes[q]]] == 1 ) #If perfect score, then only return that one.
        finalResults[["Ranking"]][[listPhenotypes[q]]]        <- finalResults[["Ranking"]][[listPhenotypes[q]]][PerfectIndices]
        finalResults[["Markers"]][[listPhenotypes[q]]]        <- finalResults[["Markers"]][[listPhenotypes[q]]][PerfectIndices]
        finalResults[["Marker_Groups"]][[listPhenotypes[q]]]  <- finalResults[["Marker_Groups"]][[listPhenotypes[q]]][PerfectIndices]
        finalResults[["Cell_Labels"]][[listPhenotypes[q]]]     <- finalResults[["Cell_Labels"]][[listPhenotypes[q]]][PerfectIndices]
        AddPlusMore <- FALSE
    }

    # change clean.res for updateLists
    indicesSort <- c()
    for ( g1 in 1 : length(finalResults[["Cell_Labels"]][[listPhenotypes[q]]]) ){
        indicesSort <- c(indicesSort, which(finalResults[["Cell_Labels"]][[listPhenotypes[q]]][g1] == clean.res[,"celllabel"]) )
    }
    if ( length(indicesSort == 1)  && nrow(clean.res) == 1 ) {
        clean.res <- clean.res
    } else {
        if (length(indicesSort) == 1 ) {
            clean.res <- t(as.matrix(clean.res[indicesSort,]))
        } else {
            clean.res <- clean.res[indicesSort,]
        }
    }

    # Compile the new row of the .csv file by updating all the list functions. This is the second time. listPhenotypes.csv uses these second case.
    r2 <- updateLists ( clean.res, MaxHitsPht, AddPlusMore)
    listMarkerLabels[[q]] <- r2[1] ; listCellLabels[[q]] <- r2[2]
    listPhenotypeID[[q]]  <- r2[3] ; listCellID[[q]]     <- r2[4]

    listMarkers[q] <- cleanMarkersList(finalResults[["Markers"]][[listPhenotypes[q]]], AddPlusMore)
    listRanking[q] <- cleanRankingList(finalResults[["Ranking"]][[listPhenotypes[q]]], AddPlusMore)

    #----------------------------------------------------------- Children, parents and treeDiagram

    # Added in to make the code faster when the user only wants to know if the markers
    # are in the ontology or not and does not want the tree diagrams
    if ( VisualizationSkip == FALSE ) {

        # Identify all parents of the matches.
        parent.res <- matrix ( nrow = 0, ncol = 5 )
        colnames ( parent.res ) <- c ( "x", "celllabel", "parent", "parentlabel", "score" )
        for ( i in 1:nrow ( clean.res ) ) {
            # Fixes a bug when trying to save cell type of "F4/80-negative adipose macrophage".
            # the "sub" replaces the "/" with a "_" to avoid R thinking that F4 is a folder.
            fname <- paste ( save.dirParentsQuery, sub ( "/","_", ( clean.res[i, "celllabel"] ) ), '.csv', sep = "" )
            if ( file.exists ( fname ) ) {  # read from file if there is one, instead of querying the same thing again
                res <- read.csv ( fname, as.is = TRUE )
            } else {
                res <- parentQuery ( child.label =  clean.res[i, "celllabel"], query.file = que.getParentClasses, prefix.info, endpoint=endpoint )
                write.csv ( res, fname, row.names = FALSE )
            }
            res <- cbind ( res, rep ( clean.res[i, "Score"], nrow ( res ) ) )
            colnames ( res ) [ ncol ( res ) ] <- "Score"
            parent.res <- rbind ( parent.res, res )
        }
        # Save parent information.
        fname <- paste ( save.dirParents, 'parent.res', phenotype, '.csv', sep = "" )
        write.csv ( parent.res, fname, row.names = FALSE )
        if ( Verbose == TRUE ) { cat ( 'Parent information saved in', fname, "\n" ) }

        # print(grep(x = parent.res[ , 'parentlabel'], pattern = "@en"))

        # January 2018 fix. For now the solution was to remove the @en labels. In the future it might make sense to uncomment the
        #  flowCL_query_data_getParentClasses parts in supportingFunctions.R
        if ( length(grep(x = parent.res[ , 'parentlabel'], pattern = "@en")) >= 1 ){
            parent.res <- parent.res[-grep(x = parent.res[ , 'parentlabel'], pattern = "@en"),]
        }

        summary <- table ( parent.res [ , 'parentlabel'] )

        scores2 <- sapply ( names ( summary ), function ( p ) mean ( as.numeric ( as.character ( parent.res [ which ( parent.res [ , 'parentlabel'] == p ), 'Score'] ) ) ) )
        # The higher up in the tree -- the more times a population is called a parent to
        # others -- the more certain we are that the label applies (e.g. 'cell' is usually
        # a parent to most population labels under investigation, and we are pretty sure
        # whatever phenotype we are working with is at least a cell!)
        # The higher the Score is, the more markers in our phenotype matched. Combining
        # these two measures gives an overall estimate of how specific and how reliable
        # the hits are.
        scored.summary <- summary / max ( summary ) * scores2 / max ( scores2 )
        s <- order ( scores2, scored.summary, decreasing = TRUE )
        # Create a list of the populations and their parents for visualization purposes.
        parent.analysis <- NULL
        for ( q1 in 1:length ( s ) ) {
            # Fixes a bug when trying to save cell type of "F4/80-negative adipose macrophage".
            # the "sub" replaces the "/" with a "_" to avoid R thinking that F4 is a folder.
            fname <- paste ( save.dirParentsQuery, sub ( "/", "_", names ( summary ) [ s[q1] ] ), '.csv', sep = "" )
            if ( file.exists ( fname ) ) { # read from file if there is one, instead of querying the same thing again
                res <- read.csv ( fname, as.is = TRUE )
            } else {
                res <- parentQuery ( child.label = names ( summary ) [ s[q1] ], query.file = que.getParentClasses, prefix.info, endpoint=endpoint )
                write.csv ( res, fname, row.names = FALSE )
            }
            parent.analysis [ length ( parent.analysis ) + 1 ] <- list ( res [ which ( is.element ( res [ , 'parentlabel'], names ( summary ) ) ), 'parentlabel'] )
        }
        names ( parent.analysis ) <- names ( summary ) [s]

        # Also create a list of the populations and their children for visualization purposes.
        child.analysis <- lapply ( names ( summary ) [s], function(x) {
            children <- c()
            for ( y in setdiff ( names ( summary ), x ) ) {
                if ( is.element ( x, parent.analysis[[y]] ) ) {
                    children <- c ( children, y )
                }
            }
            return ( children )
        })
        names ( child.analysis ) <- names ( summary ) [s]
        # Create and save a pdf file of a tree diagram

        temp.MarkerGroups <- finalResults[["Marker_Groups"]][[listPhenotypes[q]]]
        temp.MarkerGroups <- MarkerGroups

        for ( u1 in 1: length( temp.MarkerGroups ) ) {
            temp.MarkerGroups[[u1]] <- temp.MarkerGroups[[u1]][c(2,3,4)]
        }

        finalResults[[listPhenotypes[q]]] <- treeDiagram (
            child.analysis, clean.res, phenotype, OntolNamesTD, marker.list.short,
            marker.list, save.dir, listColours_flowCL = listColours_flowCL,
            MarkerGroups = temp.MarkerGroups, CellLabels = temp.string
            )

    }# end of visualization if statement

    #----------------------------------------------------------- The rest

    MarkerGroupsAmount <- NULL
    NonMatchMarkers <- NULL
    tmpMarkerGroup <- finalResults[["Marker_Groups"]][[listPhenotypes[q]]]
    for(p1 in 1: length(tmpMarkerGroup)){
        MarkerGroupsAmount[p1] <- length(tmpMarkerGroup[[p1]][[2]])
        NonMatchMarkers[p1] <- length(tmpMarkerGroup[[p1]][[1]]) + length(tmpMarkerGroup[[p1]][[3]]) + length(tmpMarkerGroup[[p1]][[4]])
    }
    temp.num1 <- which(MarkerGroupsAmount == length(unlist(marker.list)) )
    temp.num2 <- which(NonMatchMarkers == 0 )
    temp.num  <- intersect(temp.num1, temp.num2)

    # Check to see if all markers were hits (A perfect match).
    if ( length ( tmpMarkerGroup ) == 1 && length(temp.num) == 1 ) {
        listPhenotypesuccess[[q]] <- ( "Yes" )
    }
    if ( length ( tmpMarkerGroup ) > 1 && length(temp.num) > 1 ) {
        listPhenotypesuccess[[q]] <- paste( "Yes - ", length (temp.num) , " hits", sep = "" )
    }

    listPhenotypesuccess[[q]] <-  sub("5$", "5 or more", listPhenotypesuccess[[q]])

    # The time for one marker/phenotype to be queried
    if ( Verbose == TRUE ) { cat ( "Time elapsed:", timeOutput(start), "\n" ) }
    # Shown so the user knows how far the code has run
    if ( Verbose == TRUE ) { cat ( "Iterations at", which(markersToQuery==q) , "out of",length(markersToQuery), "\n" ) }
} # end of main body for-loop

#################################################################################

#----------------------------------------------------------- Save listPhenotypes.csv

# Creating a .csv file with the short name phenotypes and the ontology name phenotypes with a "Yes" or "No" indicating if it is in the ontology or not.
# Plus a list of the markers and the cell types with their ontology IDs for the cases with the maximum number of hits.
if ( length ( listPhenotypes ) == 1 ) {
    transpose.listP <- TRUE # if there is only one row
} else {
    transpose.listP <- FALSE # if there is two or more rows
}

listPhenotypes <- cbind (
    as.character(listPhenotypes),
    listPhenotypeUpdate,
    listExpPhenotypes,
    listExpPhenotypeUpdate,
    listPhenotypesuccess,
    listPhenotypeID,
    listMarkerLabels,
    listMarkers,
    listRanking,
    listCellID,
    listCellLabels
    )

colnames(listPhenotypes) <- c(
    "Short marker names",
    "Ontology marker names",
    "Experiment markers",
    "Ontology exper. names",
    "Successful Match?",
    "Marker ID",
    "Marker Label",
    "Marker Key",
    "Score (Out of 1)",
    "Cell ID",
    "Cell Label"
    )

# If there is only one marker queried the results will only have one row and this cases a small error. "transpose.listP" fixes this
if ( transpose.listP == FALSE ) { # if there is two or more rows
    listPhenotypes <- listPhenotypes[markersToQuery,]
}
# Save the results into a .csv file for outside of R viewing
fname <- paste ( save.dir, 'listPhenotypes.csv', sep = "" )
write.csv ( listPhenotypes, fname, row.names = FALSE )

#----------------------------------------------------------- Manipulate data for $Table

if ( transpose.listP == FALSE ) { # if there is two or more rows
    # Format the .csv file into a table format to show in R
    listPhenotypes[,2] <- gsub ( "\n", ", ", listPhenotypes[,2] )
    for ( q2 in 1:11 ) {
        for ( q1 in  1:MaxHitsPht ) {
            listPhenotypes[,q2] <- gsub ( paste ( "\n", q1, sep = "" ) , paste ( " ", q1, sep = "" ), listPhenotypes[,q2] )
}
        listPhenotypes[,q2] <- gsub ( "\n"  , ", ", listPhenotypes[,q2] )
        listPhenotypes[,q2] <- gsub ( ",  +", " ", listPhenotypes[,q2] )
    }
}

if ( transpose.listP == TRUE ) { # if there is only one row
    # Format the .csv file into a table format to show in R
    listPhenotypes[2] <- gsub ( "\n", ", ", listPhenotypes[2] )
    for ( q2 in 1:11 ) {
        for ( q1 in  1:MaxHitsPht ) {
            listPhenotypes[q2] <- gsub ( paste ( "\n", q1, sep = "" ) , paste ( " ", q1, sep = "" ), listPhenotypes[q2] )
        }
        listPhenotypes[q2] <- gsub ( "\n"  , ", ", listPhenotypes[q2] )
        listPhenotypes[q2] <- gsub ( ",  +", " ", listPhenotypes[q2] )
    }
}

SegmentTotal <- list()
# Create the output if []$Table is called for multiple phenotypes
if ( transpose.listP == FALSE ) { # if there is two or more rows
    for ( l1 in 1:length ( listPhenotypes[,1] ) ) {
        listTem <- t ( listPhenotypes[l1,] )
        columnof_listPhenotypes <- colnames ( listPhenotypes )
        segment <- strwrap ( listTem, 70 )
        segIndices <- c( 1, 2,  which(regexpr(paste("^", listExpPhenotypes[l1],"$", sep=""), segment) == 1 ),
                        which(regexpr(paste("^", listExpPhenotypes[l1],"$", sep=""), segment) == 1 ) + 1,
                        union ( which ( regexpr ( "No", segment ) == 1 ), which ( regexpr ( "Yes", segment ) == 1 ) ),
                        which ( regexpr("1)", segment ) == 1 ) )
#         segIndices <- c( 1,2, which(regexpr(exp.phenotype, segment) == 1 ), which(regexpr(exp.phenotype, segment) == 1 ) + 1, union ( which ( regexpr ( "No", segment ) == 1 ), which ( regexpr ( "Yes", segment ) == 1 ) ), which ( regexpr ( "1)", segment ) == 1 ) )
        temp <- rep( " ", length ( segment ) )

        temp[segIndices] <- columnof_listPhenotypes[1:length(segIndices)]
        segment <- t ( segment )
        colnames(segment) <- temp
        if ( l1 == 1 ) {
            SegmentTotal[segment[1]] <- list( t ( segment ) )
        } else {
            SegmentTotal[segment[1]] <- list( t ( segment ) )
        }
    }
    # Store results in a table
    finalResults["Table"] <- list(( SegmentTotal ))
}
# Create the output if []$Table is called for one phenotype
if ( transpose.listP == TRUE ) { # if there is only one row
    listTem <- t ( listPhenotypes )
    columnof_listPhenotypes <- colnames ( listPhenotypes )
    segment <- strwrap ( listTem, 70 )
    segIndices <- c( 1, 2,  which(regexpr(paste("^", listExpPhenotypes,"$", sep=""), segment) == 1 ),
                        which(regexpr(paste("^", listExpPhenotypes,"$", sep=""), segment) == 1 ) + 1,
                        union ( which ( regexpr ( "No", segment ) == 1 ), which ( regexpr ( "Yes", segment ) == 1 ) ),
                        which ( regexpr("1)", segment ) == 1 ) )
    temp <- rep( " ", length ( segment ) )

    temp[segIndices] <- columnof_listPhenotypes[1:length(segIndices)]

    segment <- t ( segment )
    colnames(segment) <- temp
    SegmentTotal <-  segment
    # Store results in a table
    finalResults["Table"] <- list ( t ( SegmentTotal ) )
}

#----------------------------------------------------------- End

# Removes the flowCL_results directory if KeepArch is FALSE
if ( KeepArch == FALSE ) {
    unlink ( "flowCL_results", recursive = TRUE )
}

# Output for the user
if ( Verbose == TRUE ) {
    cat ( "\nTotal time was: ", timeOutput ( initialTime ), "\n" )
    cat ( "Archive saved in \"[current directory]/", save.dir, "\"\n", sep = "" )
}

return ( finalResults )
}

