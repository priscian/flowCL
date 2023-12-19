#################################################################################
# query functions (flowCL: Semantic labelling of flow cytometric cell populations)
# Authors: Justin Meskas (jmeskas@bccrc.ca), Radina Droumeva (radina.droumeva@gmail.com )
#################################################################################
# Note:
#   The key regular expression is of the form: "CD19$|CD19[^0-9a-zA-Z][^ abes]".
#   This selects entries where the match is either of the exact marker passed
#   followed by the end of line, or followed by anything other than a number
#   or letter or followed by a space which is followed by an 'a', 'b', 'e' or 's'.
#   This will exclude alpha, beta, epsilon and single, however, it will exclude
#   any other words starting with a, b, e or s (A potential bug).
#   Research into proper use of Regular Expression" yielded a \b command
#   (i.e. \bbeta\b). Unfortunately this did not work.
#   The following scenarios will no longer occur:
#   Case 1: Looking for CD19 returns CD193, etc.
#   Case 2: Looking for Ly6 returns Ly6g, etc.
#   Case 3: Looking for CD8 returns CD8 alpha chain
#################################################################################

#################################################################################
# function used by ontologyLabel
ontolExceptions1 <- function ( marker.list, q1 ) {

    updatedException <- FALSE
#     if ( marker.list[q1] == "CD8" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "IgD" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "CD3" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "HLA-DR" ) { updatedException <- TRUE }
    return ( updatedException )
}
#################################################################################
# function used by ontologyLabel
ontolExceptions2 <- function ( marker.list, q1, Verbose, q3 ) {

#     if ( marker.list[q1] == "CD8" ) {
#         if ( Verbose == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"T cell receptor co-receptor CD8\"\n" ) }
#         marker.list[q1] <-"T cell receptor co-receptor CD8"
#     }
    if(marker.list[q1]=="IgD"){
        if ( Verbose == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"IgD immunoglobulin complex\"\n" ) }
        marker.list[q1] <-"IgD immunoglobulin complex"
    }
    if(marker.list[q1]=="CD3" && q3 == 2 ){ # negative marker
        if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[q1]), "has been FORCED to update to \"CD3 epsilon\"\n" ) }
        marker.list[q1] <-"CD3 epsilon"
    }
    if(marker.list[q1]=="CD3" && q3 != 2 ){
        if ( Verbose == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"alpha-beta T cell receptor complex\"\n" ) }
        marker.list[q1] <-"alpha-beta T cell receptor complex"
    }
    if(marker.list[q1]=="HLA-DR"){
        if ( Verbose == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"MHC class II histocompatibility antigen alpha chain HLA-DRA\"\n" ) }
        marker.list[q1] <-"MHC class II histocompatibility antigen alpha chain HLA-DRA"
    }
    return(marker.list)
}
#################################################################################
# Searches the ontology to find the correct label for each marker
ontologyLabel <- function ( marker.list, marker.list.short, Verbose="", save.dir="", que.hasProperLabel="", que.hasProperSynonym="", prefix.info="", endpoint="") {

    for ( q3 in 1:length(marker.list)){
        if ( length ( marker.list[[q3]] ) != 0 ) {
            for ( q1 in 1:length ( marker.list[[q3]] ) ) {
                # skips query if there is a file named "markers_ShortName_OntologyName" with all the marker ontology labels
                fname <-  paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
                if(file.exists(fname)){
                    markers_ShortName_OntologyName <- read.csv ( fname , header = FALSE )
                    markers_ShortName_OntologyName <- as.matrix ( markers_ShortName_OntologyName )

                    temp.marker.short <- as.character ( marker.list.short[[q3]][q1] )

                    temp.marker.short <- tempMarkerShort(temp.marker.short)

                    temp.location <- which ( markers_ShortName_OntologyName[,1] == temp.marker.short )
                    if ( length ( temp.location ) == 1 ) {
                        if ( Verbose == TRUE ) {  cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]), "is called", markers_ShortName_OntologyName[temp.location,2], "\n" ) }
                        marker.list[[q3]][q1] <- markers_ShortName_OntologyName[temp.location,2]
                        next
                    }
                }
                # Change the short name markers in marker.list to the marker labels in the ontology.

                # This is only needed for the ones that the code has trouble finding.
                # Hopefully with an updated ontology these next 8 lines can be deleted.

                updatedException <- ontolExceptions1 ( marker.list[[q3]], q1 )
                marker.list[[q3]] <- ontolExceptions2 ( marker.list[[q3]], q1, Verbose, q3)

                if ( updatedException == TRUE ) {
                    # Update the markers_ShortName_OntologyName.csv file
                    UpdateMarkers_ShortName_OntologyName(save.dir=save.dir, marker.list.short=marker.list.short, marker.list=marker.list, q3=q3, q1=q1)
                    next
                }

                #------------------------------------------------------------- Query exact match in synonym
                temp.marker <- marker.list[[q3]][q1]

                res <- queryMarker ( marker = temp.marker, query.file = que.hasProperSynonym, prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=FALSE )

                if ( nrow ( res ) == 1 ) {
                    temp.marker <- res[1,'label']
                    if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]), "has been updated to", temp.marker, "\n" ) }
                    marker.list[[q3]][q1] <- temp.marker
                    # Update the markers_ShortName_OntologyName.csv file
                    UpdateMarkers_ShortName_OntologyName(save.dir=save.dir, marker.list.short=marker.list.short, marker.list=marker.list, q3=q3, q1=q1)
                    next
                }

                #------------------------------------------------------------- Query match in Label
                temp.marker <- marker.list[[q3]][q1]

                res <- queryMarker ( marker = temp.marker, query.file = que.hasProperLabel, prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=FALSE )

                if ( nrow ( res ) == 1 ) {
                    temp.marker <- res[1,'label']
                    if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]), "has been updated to", temp.marker, "\n" ) }
                    marker.list[[q3]][q1] <- temp.marker

                    # Update the markers_ShortName_OntologyName.csv file
                    UpdateMarkers_ShortName_OntologyName(save.dir=save.dir, marker.list.short=marker.list.short, marker.list=marker.list, q3=q3, q1=q1)
                }
                if ( nrow ( res ) >= 2 ) {
                    if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]),"has not been changed, multiple possible markers in label\n" ) }
                }
                if ( nrow ( res ) == 0 | nrow ( res ) >= 2 ) {
                    temp.marker <- marker.list[[q3]][q1]
                    #------------------------------------------------------------- Query match in synonym


                    res <- queryMarker ( marker = temp.marker, query.file = que.hasProperSynonym, prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=FALSE )

                    # Small loop to check if the query is giving multiple label names. In this case the marker will not be changed
                    if ( nrow ( res ) >= 1 ) {
                        temp = res[1,'label']
                        for ( q2 in 1:nrow ( res ) ) {
                            if ( temp == res[q2,'label'] ) {
                                temp = res[q2,2]
                                sameLabels <- TRUE;
                            } else {
                                sameLabels <- FALSE;
                                break
                            }
                        }
                        if ( sameLabels == TRUE ) { # label exists and there is only one label
                            temp.marker <- res[1,'label']
                            if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]), "has been updated to", temp.marker, "\n" ) }
                            marker.list[[q3]][q1] <-temp.marker

                            # Update the markers_ShortName_OntologyName.csv file
                            UpdateMarkers_ShortName_OntologyName(save.dir=save.dir, marker.list.short=marker.list.short, marker.list=marker.list, q3=q3, q1=q1)
                        }
                        if ( sameLabels == FALSE ) { # label exists however there is two or more labels
                            if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]), "has not been changed, multiple possible markers in synonyms\n" ) }
                        }
                    }
                    if ( nrow ( res ) == 0 ) { # label does not exist
                        if ( Verbose == TRUE ) { cat ( "Marker", gsub("CD3e", "CD3", marker.list[[q3]][q1]), "could not be found\n" ) }
                    }
                }
            } # end of for-loop
        } # end of if statement
    }
    return ( marker.list )
}

#################################################################################
# Similar to queryMarker, but for finding the parents for a specific cell population
# by cell label.
# NOTE: I tried matching by ID and for some reason I could do it on the HE group
# but not here for some reason, that's why I'm matching by label...
parentQuery <- function ( child.label = "common myeloid progenitor",
                        query.file = "getParentClasses.txt", prefix.info="" , endpoint="") {

    # Concatenate the query preceded by all prefix information as a single string
    # to be passed to SPARQL.
    query <- paste ( c ( prefix.info, query.file ), collapse="\n" )

    # child.label <- gsub ( "\"",  "", child.label)
    # child.label <- gsub ( "@en", "", child.label)

    # Add "^" to beginning and "$" to end to find an exact match for the label.
    child.label <- paste ( "^", child.label, "$", sep = "", collapse = "" )
    query <- gsub ( "\\$label", child.label, query )
    # Execute query

    res <- SPARQL ( url = endpoint, query )$results

    return ( res )
}

#################################################################################
# Makes a query with SPARQL to the CL Ontology, and returns the results.
queryMarker <- function ( marker = NULL, query.file = "getMatchingSynonyms.txt",
                        celllabel = NULL, prefix.info="", Verbose="", endpoint="", exactMatch=FALSE) {

    # TO DO: improve input check
    # For now, a NULL, NA or length 0 phenotype will have empty results
    #  matrix returned.
#     if ( !is.character ( marker ) || length ( marker ) == 0 ) {
#         warning ( "No marker found." )
#         res <- matrix ( ncol = 2, nrow = 0 )
#         colnames ( res ) <- c ( "ID", "Synonym Match" )
#         return ( res )
#     }

    # Concatenate the query preceded by all prefix information as a single string
    # to be passed to SPARQL.
    query <- paste ( c ( prefix.info, query.file ), collapse="\n" )

    # Prepare marker for query by ensuring the marker is either followed by the
    # end of the line or it has a symbol other than a letter, a number or a space with
    # 'a', 'b', 'e' or 's' after it, as that may indicate a different marker.
    if ( !is.null ( marker ) ) {
        if ( exactMatch == TRUE ) {
            marker <- paste ( "^", marker, "$", collapse="", sep = "" )
            query <- gsub ( "\\$marker", marker, query )
        } else {
            marker <- paste ( marker, "$|", marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
    #         marker <- paste ("^", marker, "$|", marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
            query <- gsub ( "\\$marker", marker, query )
        }
    }

    if ( !is.null ( celllabel ) ) {
        celllabel <- paste ( "^", celllabel, "$", sep = "", collapse = "" )
        query <- gsub ( "\\$celllabel", celllabel, query )
    }

    # Execute query
    res <- SPARQL ( url = endpoint, query )$results

    browser()
    return ( res )
}



