#################################################################################
# supporting functions (flowCL: Semantic labelling of flow cytometric cell populations)
# Authors: Justin Meskas (jmeskas@bccrc.ca), Radina Droumeva (radina.droumeva@gmail.com )
#################################################################################


#################################################################################
# A test to see if the user has input "cd" instead of "CD"
cdTest <- function ( MarkerList ) {
    ShouldBreak <- FALSE
    for ( p1 in 1:length(MarkerList)){
        if ( regexpr ( pattern="cd", MarkerList[p1] )[[1]] >= 1) {
            warning ( paste("You have input \"cd\" in ", MarkerList[p1], ". flowCL is case sensitive. CD should be capitalized.", sep="") )
            ShouldBreak <- TRUE
        }
        if ( regexpr ( pattern="cD", MarkerList[p1] )[[1]] >= 1) {
            warning ( paste("You have input \"cD\" in ", MarkerList[p1], ". flowCL is case sensitive. CD should be capitalized.", sep="") )
            ShouldBreak <- TRUE
        }
        if ( regexpr ( pattern="Cd", MarkerList[p1] )[[1]] >= 1) {
            warning ( paste("You have input \"Cd\" in ", MarkerList[p1], ". flowCL is case sensitive. CD should be capitalized.", sep="") )
            ShouldBreak <- TRUE
        }
    }
    return (ShouldBreak)
}
#################################################################################
# Change the "." in HLA.DR to a - since the + and - signs are reserved for splitting the string
# Consider to change this to allow for any marker with a period to change to a dash
changeHLADR <- function ( marker.list ) {

    for ( q2 in 1:length(marker.list)){
        if ( length ( marker.list[[q2]] ) >= 1 ) {
            for ( q1 in 1:length ( marker.list[[q2]] ) ) {
                if ( marker.list[[q2]][q1] == "HLA.DR" ) {
                    marker.list[[q2]][q1] <- "HLA-DR"
                }
            }
        }
    }
    return ( marker.list )
}

#################################################################################
# Prepares the marker list to be put into listPhenotype and $Table
cleanMarkersList <- function ( temp, AddPlusMore ) {

    fullname <- c()
    for ( y1 in 1: length(temp)){
        fullname <- paste(fullname , y1, ") ", temp[y1] , "\n", sep="")
    }
    if ( AddPlusMore == FALSE ) {
        fullname <- substr ( fullname, 1, nchar ( fullname ) - 1 )
    } else {
        fullname <- paste ( fullname,  " + more", sep="" )
    }
    return ( fullname )
}

#################################################################################
# Prepares the ranking list to be put into listPhenotype and $Table
cleanRankingList <- function ( temp, AddPlusMore ) {

    fullname <- c()
    for ( y1 in 1: length(temp)){
        fullname <- paste(fullname , y1, ") ", round(temp[y1], digits = 3) , "\n", sep="")
    }
    if ( AddPlusMore == FALSE ) {
        fullname <- substr ( fullname, 1, nchar ( fullname ) - 1 )
    } else {
        fullname <- paste ( fullname,  " + more", sep="" )
    }
    return ( fullname )
}

#################################################################################
# Organize which markers are required, which are extra, which are in the experiment and not used,
# and which additional ones would be required to get a perfect match.
# reverse query
MarkerGroupsFunc <- function ( temp.string, query.dir.list, postfix, save.dir, query.file.list, prefix.info, Verbose, marker.list.short, exp.marker.list.short, AllMarkerGroups, endpoint=endpoint ) {

    MarkerGroups <- list()
    for ( y1 in 1 : length(temp.string) ){
        temp.res.all <- c()
        for ( y2 in 1:4 ) {
            fname <- paste( query.dir.list[y2], sub("/", "_", temp.string[[y1]][1]), postfix[y2], sep="")
#             # Fixes a bug when trying to save cell type of "F4/80-negative adipose macrophage".
#             # the "sub" replaces the "/" with a "_" to avoid R thinking that F4 is a folder.
#             fname <- paste ( save.dirParentsQuery, sub ( "/","_", ( clean.res[i, "celllabel"] ) ), '.csv', sep = "" )
            if( !file.exists ( fname ) ) {
                if ( Verbose == TRUE && y2 == 1 ) { cat ( "Querying cell label:", temp.string[[y1]][1], "\n" ) }
                temp.res <- queryMarker ( celllabel = temp.string[[y1]][1], query.file = query.file.list[[y2]], prefix.info = prefix.info, Verbose=Verbose, endpoint=endpoint, exactMatch=TRUE)
                if ( nrow(temp.res) == 0 ){ # If nothing create a non-blank csv file
                    temp.res <- matrix ( ncol = 2, nrow = 0 )
                    colnames ( temp.res ) <- c ( "ID", "Synonym Match" )
                }

                temp.loc <- which(colnames(temp.res) == "pl")
                if (length(temp.loc) >= 1){
                    colnames(temp.res)[which(colnames(temp.res) == "pl")] <- "plabel"
                }

                # Create 4 csv files for each returned cell label (has, lacks, low, high)
                write.csv ( temp.res, fname, row.names = FALSE )
            }
            else {
                temp.res <- read.csv ( fname, as.is=TRUE )
            }
            temp.res.all <- rbind(temp.res.all, temp.res) # combine all 4 csv files
        }

        fname <- paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
        markers_ShortName_OntologyName <- read.csv ( fname , header = FALSE )
        markers_ShortName_OntologyName <- as.matrix ( markers_ShortName_OntologyName )
        for ( p1 in 1 : nrow(temp.res.all)) {
            markerIndices <- which ( markers_ShortName_OntologyName[,2] == temp.res.all[p1,"markerlabel"] )
            if ( length ( markerIndices ) == 0 ) { # If the Ontology marker name is not in the csv file then use the Ontology marker name instead of the short marker label
                temp.res.all[p1,4] <- temp.res.all[p1,"markerlabel"]
            } else {
                if( length( markerIndices ) == 1) { # Use the short marker label in csv file
                    temp.res.all[p1,4] <- markers_ShortName_OntologyName[markerIndices , 1]
                } else { # pick the one from markers_ShortName_OntologyName that has the lowest amount of characters
                    temp.res.all[p1,4] <- markers_ShortName_OntologyName[markerIndices[which ( min(nchar(markers_ShortName_OntologyName[markerIndices , 1])) == nchar(markers_ShortName_OntologyName[markerIndices , 1]))], 1]
                }
            }
        }
        # Add tags
        suppressWarnings(temp.res.all[,6] <- temp.res.all[,4])
        if( length(which(temp.res.all[,6]=="CD3z")) >=1 ) {
            temp.res.all[which(temp.res.all[,6]=="CD3z"), 6] <- "CD3"
        }
        if( length(which(temp.res.all[,6]=="CD3e")) >=1 ) {
            temp.res.all[which(temp.res.all[,6]=="CD3e"), 6] <- "CD3"
        }
        suppressWarnings(temp.res.all[which(temp.res.all[,3] == "lacks_plasma_membrane_part")     ,6] <- paste(temp.res.all[which(temp.res.all[,3] == "lacks_plasma_membrane_part"),     6], "-", sep="")  )

        # suppressWarnings(temp.res.all[which(temp.res.all[,3] == "\"has plasma membrane part\"@en"),6] <- paste(temp.res.all[which(temp.res.all[,3] == "\"has plasma membrane part\"@en"),6], "+", sep="")  )
        suppressWarnings(temp.res.all[which(temp.res.all[,3] == "has plasma membrane part"),       6] <- paste(temp.res.all[which(temp.res.all[,3] == "has plasma membrane part"),       6], "+", sep="")  )

        suppressWarnings(temp.res.all[which(temp.res.all[,3] == "has_low_plasma_membrane_amount") ,6] <- paste(temp.res.all[which(temp.res.all[,3] == "has_low_plasma_membrane_amount"), 6], "lo", sep="") )
        suppressWarnings(temp.res.all[which(temp.res.all[,3] == "has_high_plasma_membrane_amount"),6] <- paste(temp.res.all[which(temp.res.all[,3] == "has_high_plasma_membrane_amount"),6], "hi", sep="") )
        # remove all duplicates of CD3's
        if( length(which(temp.res.all[,6]=="CD3+")) > 1 ) {
            temp.res.all <- temp.res.all[-which(temp.res.all[,6]=="CD3+")[-1], ]
        }
        if( length(which(temp.res.all[,6]=="CD3hi")) > 1 ) {
            temp.res.all <- temp.res.all[-which(temp.res.all[,6]=="CD3hi")[-1], ]
        }
        if( length(which(temp.res.all[,6]=="CD3lo")) > 1 ) {
            temp.res.all <- temp.res.all[-which(temp.res.all[,6]=="CD3lo")[-1], ]
        }
        temp.list.tags <- setdiff(temp.res.all[,6], unlist( marker.list.short ) )
        temp.list.No.tags <- temp.list.tags

        # Remove tags
        temp.list.No.tags <- gsub("[+]$",  "", temp.list.No.tags); temp.list.No.tags <- gsub("-$",  "", temp.list.No.tags)
        temp.list.No.tags <- gsub("hi$", "", temp.list.No.tags);   temp.list.No.tags <- gsub("lo$", "", temp.list.No.tags)

        AdditionalMarkers <- setdiff(  temp.list.No.tags, unlist( exp.marker.list.short ) ) # markers that are part of the cell label marker list but not part of the expermental marker list
        ExpNeededMarkers  <- intersect(temp.list.No.tags, unlist( exp.marker.list.short ) ) # markers that are part of the cell label marker list and part of the expermental marker list

        # Put the tags back onto the marker labels
        if(length(AdditionalMarkers) > 0 && length(temp.list.tags) > 0 ){
            for ( w1 in 1 : length(AdditionalMarkers)){
                for ( w2 in 1 : length(temp.list.tags)){
                    if( grepl(AdditionalMarkers[w1], temp.list.tags[w2])){
                        AdditionalMarkers[w1] <- temp.list.tags[w2]
                    }
                }
            }
        }

        # Put the tags back onto the marker labels
        if(length(ExpNeededMarkers) > 0 && length(temp.list.tags) > 0 ){
            for ( w1 in 1 : length(ExpNeededMarkers)){
                for ( w2 in 1 : length(temp.list.tags)){
                    if( grepl(ExpNeededMarkers[w1], temp.list.tags[w2])){
                        ExpNeededMarkers[w1] <- temp.list.tags[w2]
                    }
                }
            }
        }

        NeededMarkers <- intersect(temp.res.all[,6], unlist( marker.list.short ) ) # markers that are part of the cell label marker list and are part of the input marker list
        UnneededMarkers <- setdiff( unlist( marker.list.short ), temp.res.all[,6] ) # markers that are part of the input marker list but not part of the cell label marker list

        # check to see if there are copies of the same marker. This happens when has_pmp is in the CL at the same time as low_PMA or high_PMA
        if ( length(ExpNeededMarkers) >= 1 ) {
            tmp1 <- NeededMarkers
            tmp1 <- gsub("[+]$",  "", tmp1); tmp1 <- gsub("-$",  "", tmp1); tmp1 <- gsub("hi$", "", tmp1);   tmp1 <- gsub("lo$", "", tmp1)
            tmp2 <- ExpNeededMarkers
            tmp2 <- gsub("[+]$",  "", tmp2); tmp2 <- gsub("-$",  "", tmp2); tmp2 <- gsub("hi$", "", tmp2);   tmp2 <- gsub("lo$", "", tmp2)
            RemoveIndices <- NULL

            # for problems with has_pmp in the CL at the same time as low_PMA or high_PMA
            for ( w1 in 1 : length(tmp2) ) {
                if( length(which (tmp1 == tmp2[w1]) >= 1 ) ) {
                    RemoveIndices <- c(RemoveIndices, w1)
                }
            }
            if ( !is.null(RemoveIndices) ) { # Remove extra marker restriction
                ExpNeededMarkers <- ExpNeededMarkers[-RemoveIndices]
            }
        }
        if ( AllMarkerGroups == FALSE ){
            MarkerGroups[[y1]] <- list(UnneededMarkers, NeededMarkers, ExpNeededMarkers, AdditionalMarkers)
        } else {
            MarkerGroups[[y1]] <- list(NeededMarkers, ExpNeededMarkers, AdditionalMarkers)
        }
    }
    return(MarkerGroups)
}
#################################################################################
# Removes tags and leaves the hyphen in HLA-DR
tempMarkerShort <- function (temp.marker.short) {

    temp.marker.short <- sub("HLA-DR", "HLA.DR", temp.marker.short);
    temp.marker.short <- sub("-", "", temp.marker.short);
    temp.marker.short <- sub("[+]", "", temp.marker.short);
    temp.marker.short <- sub("lo", "", temp.marker.short);
    temp.marker.short <- sub("hi", "", temp.marker.short);
    temp.marker.short <- sub("HLA[.]DR", "HLA-DR", temp.marker.short);
    return(temp.marker.short)
}

#################################################################################
# Removes all the non first generation children from the child.analysis and then
# produces a flow chart.

treeDiagram <- function ( child.analysis, clean.res, phenotype, OntolNamesTD, marker.list.short, marker.list, save.dir, listColours_flowCL = "" ,
                            MarkerGroups = "", CellLabels = "" ) {


    # Sort the child.analysis by starting with the cell population which has
    # the most children to the one with the least-- i.e. is the most likely to be the
    # root parent, such as 'cell' or 'native cell':
    child.lengths <- unlist ( lapply(child.analysis, length ) )
    sort.child <- child.analysis[order ( child.lengths, decreasing = TRUE )]
    labels <- names ( sort.child )
    arrows.colour <- arrows.label <- arrows.dashed.solid <- arrows.head <- NULL
    if ( OntolNamesTD == TRUE ) { marker.list.short <- marker.list }

    # A quadruple for loop (probably not the most elegant, but it is still fast ). Start with q8,
    # which selects the label with the most children first. q9 stores an element from the q8 label. q10
    # and q11 look through all other elements in the q8 label and detects if the stored
    # element (q9) is a parent of other elements in the q8 label. If this is true, then the element
    # q11 is removed from the q8 label.
    # In short, something that is not a first generation child will be removed from the list.
    for ( q8 in 1:length ( labels ) ) {
        deletePoints <- NULL
        for ( q9 in 1:length ( child.analysis[labels[q8]][[1]] ) ) {

            if ( is.null ( child.analysis[labels[q8]][[1]][q9] ) ) { next }

            temp.label <- child.analysis[labels[q8]][[1]][q9]

            for ( q10 in 1:length ( child.analysis[temp.label][[1]] ) ) {
                for ( q11 in 1:length ( child.analysis[labels[q8]][[1]] ) ) {

                    if ( !is.null ( child.analysis[temp.label][[1]][q10] ) ) {
                        if ( ( child.analysis[temp.label][[1]][q10] == child.analysis[labels[q8]][[1]][q11] ) & ( q9!=q11 ) ) {
                            # stores points for removal
                            deletePoints <- c ( deletePoints, q11 )
                        }
                    }
                }
            }
        }
        if ( !is.null ( deletePoints ) ) {
            # removes all non first generation children
            child.analysis[labels[q8]][[1]] <- child.analysis[labels[q8]][[1]][-deletePoints]
        }
    }
    # Sets up Rgraphviz with all the node titles in "labels"
    rEG <- new ( "graphNEL", nodes = c ( labels, marker.list.short[[1]], marker.list.short[[2]], marker.list.short[[3]], marker.list.short[[4]]), edgemode="directed" )

    for ( q13 in 1:length ( labels ) ) {
        for ( q14 in 1:length ( child.analysis[labels[q13]][[1]] ) ) {
            if ( !is.null ( child.analysis[labels[q13]][[1]][q14] ) ) {
                # Adds an arrow from parent to child
                rEG <- addEdge ( child.analysis[labels[q13]][[1]][q14], labels[q13], rEG, 1 )
                arrows.colour <- c ( arrows.colour, "black" )
                arrows.label  <- c ( arrows.label , paste ( child.analysis[labels[q13]][[1]][q14], "~", labels[q13], sep = "" ) )
                arrows.dashed.solid <- c ( arrows.dashed.solid, "solid" )
                arrows.head <- c ( arrows.head, "open" )
            }
        }
    }
    labels.colour <- labels

    # Load all predetermined colours
    listColours_flowCL <- as.matrix ( listColours_flowCL )

    # Calculate the number of markers in each subsection of marker.list.short and calculate the indices
    NumMarkers <- c( length(marker.list.short[[1]]), length(marker.list.short[[2]]), length(marker.list.short[[3]]), length(marker.list.short[[4]]) )
    NumMarkersOrg <- NumMarkers
    NumMarkers[2] <- NumMarkers[1] + NumMarkers[2]
    NumMarkers[3] <- NumMarkers[2] + NumMarkers[3]
    NumMarkers[4] <- NumMarkers[3] + NumMarkers[4]
    NumMarkers    <- NumMarkers - NumMarkersOrg
    # Make lists of arrows start and end locations, colour and label
    count.markers <- 0
    for ( q15 in 1:length ( labels ) ) {
        for ( q16 in 1:length ( clean.res[,'celllabel'] ) ) {
            if ( clean.res[q16,'celllabel'] == ( labels[q15] ) ) {
                for ( q3 in 1:length(marker.list)){
                    if ( length ( marker.list[[q3]] ) != 0 ) {
                        for ( q17 in 1:length ( marker.list[[q3]] ) ) {
                            if ( grepl ( marker.list[[q3]][q17], clean.res[q16,'markerlabel'] ) == TRUE ) {
                                count.markers <- count.markers + 1
                                rEG <- addEdge ( marker.list.short[[q3]][q17], labels[q15], rEG, 1 )
                                arrows.colour <- c ( arrows.colour, listColours_flowCL[NumMarkers[q3]+q17] )
                                arrows.label  <- c ( arrows.label, paste ( marker.list.short[[q3]][q17],"~",labels[q15], sep = "" ) )
                                arrows.dashed.solid <- c ( arrows.dashed.solid, "dashed" )
                                arrows.head <- c ( arrows.head, "none" )
                            }
                        }
                    }
                }

                # make perfect matches green
                if ( length ( which(CellLabels == labels[q15]) ) != 0) {
                    if ( count.markers == length(unlist(MarkerGroups[[ as.numeric ( which(CellLabels == labels[q15]) )]] )) ) {
                        labels.colour[q15] <- "lightgreen"
                        count.markers <- 0
                    } else { # colour the partial matches a beige colour
                        labels.colour[q15] <- "bisque"
                        count.markers <- 0
                    }
                }

                break
            } else  { # colour all the non-important parents white
                labels.colour[q15] <- "white"
            }
        }
    }

    if ( length ( marker.list[[1]] ) != 0 ) { # colour positive markers sky blue
        for ( q19 in 1:( length ( marker.list[[1]] ) ) ) {
            labels.colour[length ( labels.colour ) + 1] <- "skyblue"
        }
    }
    if ( length ( marker.list[[2]] ) != 0 ) {
        for (q19 in 1:(length ( marker.list[[2]] ) ) ) { # colour negative markers a light red colour
            labels.colour[length ( labels.colour ) + 1] <- "lightcoral"
        }
    }
    if ( length ( marker.list[[3]] ) != 0 ) {
        for (q19 in 1:(length ( marker.list[[3]] ) ) ) { # colour negative markers a lighter blue colour
            labels.colour[length ( labels.colour ) + 1] <- "royalblue1"
        }
    }
    if ( length ( marker.list[[4]] ) != 0 ) {
        for (q19 in 1:(length ( marker.list[[4]] ) ) ) { # colour negative markers a darker blue colour
            labels.colour[length ( labels.colour ) + 1] <- "lightcyan3"
        }
    }
    nAttrs <- list ( )
    eAttrs <- list ( )

    nAttrs$fillcolor <- structure(c (labels.colour ), .Names = c (labels, marker.list.short[[1]], marker.list.short[[2]], marker.list.short[[3]], marker.list.short[[4]] ) )
    eAttrs$color <- structure(c (arrows.colour ), .Names = c (arrows.label ) )

#   # Even though it looks like arrows.dashed.solid and arrows.head are being used, they are not.
#   # This is how the vignette for Rgraphiz implements these, however it does not work currently.
#   eAttrs$style <- structure(c ( arrows.dashed.solid ), .Names = c ( arrows.label ) )
#   eAttrs$arrowhead <- structure(c ( arrows.head ), .Names = c ( arrows.label ) )

    attrs <- list( node = list ( shape="ellipse", fontsize = 14, fixedsize=FALSE ), graph = list ( rankdir = "BT" ) )
    # create a pdf flow chart
    child.file.name <- paste ( save.dir, "tree_", phenotype, ".pdf", sep = "" )

    pdf ( file = child.file.name, width = 10.5, height = 8 )
    suppressWarnings ( plot ( rEG, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs ) )
    dev.off ( )
    return ( list ( rEG, nAttrs, eAttrs, attrs ) )
}

#################################################################################
### Break up a phenotype into individual markers and 'signs'

phenoParse <- function ( phenotype )  {

    if ( !is.character ( phenotype ) ) {
        warning ( "Phenotype is not a valid string!" )
        try ( phenotype <- as.character ( phenotype ) )
    }

    # First split up the string based on +, -, hi or lo to get the markers
    markers <- unlist ( strsplit ( x = phenotype, split = "\\+|\\-|hi|lo" ) )

    for ( k in 1:length(markers)){
        if ( nchar(markers[k]) <= 2 ) {
            warning ( "Your input markers are most likely ill defined. Most markers have labels and synonyms that are at least 3 characters long." )
        }
    }

    # Next, split the original string based on the markers found above, leaving
    # only the signs (remove first result as there is no leading sign )
    signs <- unlist ( strsplit ( x = phenotype, split = paste ( markers, sep = "", collapse="|" ) ) )[-1]

    if (length(signs) != (length(markers))) {
        warning ( "The number of tags does not match the number of markers." )
        return ()
    }

    # Return a list of positive, negative, low and high markers
    res <- list (   `Positive` = markers[signs == "+"],  `Negative` = markers[signs == "-"] ,
                `High\\Bright` = markers[signs == "hi"], `Low\\Dim` = markers[signs == "lo"])
    return ( res )
}

#################################################################################
# Creates a string with the updated phenotypes to have better formatting for listPhenotypes.csv
phenoUnparse <- function ( phenotype, marker.list ) {

    # Find indices of where a certain pattern is found in the variable called phenotype
    posit.pos <- gregexpr ( pattern="[+]", phenotype )[[1]]
    posit.neg <- gregexpr ( pattern="[-]", phenotype )[[1]]
    posit.hi  <- gregexpr ( pattern="hi",  phenotype )[[1]]
    posit.lo  <- gregexpr ( pattern="lo",  phenotype )[[1]]

    if ( sum ( posit.pos, posit.neg, posit.hi, posit.lo ) == -4){
        warning ( "Phenotype does not contain a valid tag! (e.g. +, -, lo or hi)" )
        return ( phenotype )
    }

    # Combine all indices into one list
    allPosit <- sort(c(posit.pos, posit.neg, posit.hi ,posit.lo), decreasing=FALSE)
    # Remove all -1's caused by there being no markers in one or more of the sub-groups.
    if ( length ( which (allPosit == -1 ) ) >= 1 ){
        allPosit <- allPosit[-which(allPosit == -1)]
    }

    temp.string <- c ( )
    # Compile a string that maintains the same order from the inputted phenotype but makes it friendly for .csv files
    for ( k in 1:length(allPosit)){
        if ( length ( intersect ( allPosit[k], posit.pos ) ) >= 1 ) {
            temp.string <-  paste ( temp.string, marker.list[[1]][ which ( posit.pos == intersect ( allPosit[k], posit.pos))], "\n", sep = "" )
        }
        if ( length ( intersect ( allPosit[k], posit.neg ) ) >= 1 ) {
            temp.string <-  paste ( temp.string, marker.list[[2]][ which ( posit.neg == intersect ( allPosit[k], posit.neg))], "\n", sep = "" )
        }
        if ( length ( intersect ( allPosit[k], posit.hi  ) ) >= 1 ) {
            temp.string <-  paste ( temp.string, marker.list[[3]][ which ( posit.hi  == intersect ( allPosit[k], posit.hi ))], "\n", sep = "" )
        }
        if ( length ( intersect ( allPosit[k], posit.lo  ) ) >= 1 ) {
            temp.string <-  paste ( temp.string, marker.list[[4]][ which ( posit.lo  == intersect ( allPosit[k], posit.lo ))], "\n", sep = "" )
        }
    }

    # removes the \n from the last line for formatting reasons in the .csv file
    temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 )

    return ( temp.string )
}

#################################################################################
### Condense a results table to the unique entries only and tabulate repeated hits
#   TO DO: For now, this relies on the first column having the unique IDs of
#   the owl objects returned.
tabulateResults <- function ( res ) {

    number.of.hits <- table ( res[ , 1] )
    # fixed bug here, added as.vector()
    res <- cbind ( res, as.vector ( number.of.hits[res[, 1]] ) )

    colnames ( res )[ncol ( res )] <- "Number Of Hits"
    unique.ids <- unique ( res[, 1] )

    clean.res <- matrix ( "", nrow = length ( unique.ids ), ncol = ncol ( res ) + 1 )
    colnames ( clean.res ) <- c ( colnames ( res ), "Score" )
    rownames ( clean.res ) <- unique.ids
    for ( id in unique.ids ) {
        locate.entries <- which ( res[, 1] == id )
        clean.res[id, ] <- c ( id, res[locate.entries[1], 2],
                                apply(res[locate.entries, 3:( ncol ( res ) - 2 )], 2, paste,
                                collapse = "\n" ), sum ( res[locate.entries, 'penalties'] ),
                                res[locate.entries[1], ncol(res )], 0 )
        clean.res[id, "Score"] <- as.numeric ( clean.res[id, "Number Of Hits"] ) +
        sum(as.numeric(unlist(strsplit(clean.res[id, "penalties"], split="\n"))))
    }

    # Required because of a bug caused from only having 1 row (f[1,] notation was a problem )
    if ( nrow ( clean.res ) >= 2 ) {
        sort.scores <- sort ( as.numeric ( clean.res[, "Score"] ),
                                decreasing = TRUE, index.return = TRUE )$ix
        return ( clean.res[sort.scores, ] )
    } else  {
        return ( clean.res )
    }
}

#################################################################################
# Prints out the time since start_time. Used for optimizing code and for informing the user how long certain processes take.
# This function was copied from an online forum. Useful for keeping time.
timeOutput <- function ( start_time )  {
    start_time <- as.POSIXct ( start_time )
    dt <- difftime ( Sys.time ( ), start_time, units = "secs" )
    # Since you only want the H:M:S, we can ignore the date...
    # but you have to be careful about time-zone issues
    format ( .POSIXct ( dt, tz = "GMT" ), "%H:%M:%S" )
}
timeOutput ( Sys.Date ( ) )

#################################################################################
# Finds and stores information for display in a .csv file
updateLists <- function ( clean.res, MaxHitsPht, AddPlusMore ) {
    # Creates the lists MarkerLabels, CellLabels, PhenotypeID and CellID
    BreakTrue <- FALSE
    listMarkerLabels.temp <- listCellLabels.temp <- listPhenotypeID.temp <- listCellID.temp <- NULL

    if ( max ( as.numeric ( clean.res[ , 'Number Of Hits'] ) ) >= 1 ) {
        for ( q2 in 1:min ( length ( clean.res[ ,'Number Of Hits'] ), MaxHitsPht ) ) {

#             if ( q2 == 1 & ( clean.res[q2, 'Number Of Hits'] ) == max ( as.numeric ( clean.res[ , 'Number Of Hits'] ) ) ) {
            if ( q2 == 1 ) {
                listMarkerLabels.temp <- paste ( q2, ") ", ( clean.res[q2,'markerlabel'] ), "\n", sep = "" )
                listCellLabels.temp   <- paste ( q2, ") ", ( clean.res[q2,'celllabel'] ), "\n", sep = "" )
                # Look in the Label of the marker ID for "PR" then pulls out the PR and the numbers that follow
                temp.index <- gregexpr ( "PR", ( clean.res[q2,'marker'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 11 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n
                    }
                }
                # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
                temp.index <- gregexpr ( "GO", ( clean.res[q2,'marker'] ) )
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    if ( temp.string != "" ) {
                        temp.string <- paste ( temp.string, "\n", sep = "" )
                    }
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n
                    }
                }
                # Stores this string into the lists
                listPhenotypeID.temp  <- paste ( q2, ") ", temp.string, "\n", sep = "" )

                # Look in the Label of the Cell ID for "CL" then pulls out the CL and the numbers that follow
                temp.index <- gregexpr ( "CL", ( clean.res[q2,'x'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'x'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n
                    }
                }
                # Stores this string into the lists
                listCellID.temp  <- paste ( q2,") ", temp.string, "\n", sep = "" )
                next
            }
#             if ( q2 > 1 & ( clean.res[q2,'Number Of Hits'] ) == max ( clean.res[,'Number Of Hits'] ) ) {
            if ( q2 > 1 ) {
                listMarkerLabels.temp <- paste ( listMarkerLabels.temp , q2, ") ", ( clean.res[q2,'markerlabel'] ), "\n", sep = "" )
                listCellLabels.temp   <- paste ( listCellLabels.temp   , q2, ") ", ( clean.res[q2,'celllabel'] ), "\n", sep = "" )

                # Look in the Label of the marker ID for "PR" then pulls out the PR and the numbers that follow
                temp.index <- gregexpr ( "PR", ( clean.res[q2,'marker'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 11 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n
                    }
                }
                # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
                temp.index <- gregexpr ( "GO", ( clean.res[q2,'marker'] ) )
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    if ( temp.string != "" ) {
                        temp.string <- paste ( temp.string, "\n", sep = "" )
                    }
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n
                    }
                }
                # Stores this string into the lists
                listPhenotypeID.temp  <- paste ( listPhenotypeID.temp  , q2, ") ", temp.string,"\n", sep = "" )

                # Look in the Label of the Cell ID for "CL" then pulls out the CL and the numbers that follow
                temp.index <- gregexpr ( "CL", ( clean.res[q2,'x'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'x'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n
                    }
                }
                # Stores this string into the lists
                listCellID.temp  <- paste ( listCellID.temp  , q2, ") ", temp.string, "\n", sep = "" )
                next
            }
#             if ( ( clean.res[q2,'Number Of Hits'] ) != max ( as.numeric ( clean.res[ , 'Number Of Hits'] ) ) ) {
#             if ( ( clean.res[q2,'Number Of Hits'] ) != max ( as.numeric ( clean.res[ , 'Number Of Hits'] ) ) ) {
                BreakTrue <- TRUE
                break
#             }
        }# end of for loop
    }# end of if statement

    # If there is more than 5 elements to store, then the rest is cut off and a "+ more" is displayed
    if ( length ( as.numeric ( clean.res[ , 'Number Of Hits'] ) ) > MaxHitsPht & BreakTrue == FALSE ) {
        listMarkerLabels.temp <- paste ( listMarkerLabels.temp, "+ more" )
        listCellLabels.temp   <- paste ( listCellLabels.temp,   "+ more" )
        listPhenotypeID.temp  <- paste ( listPhenotypeID.temp,  "+ more" )
        listCellID.temp       <- paste ( listCellID.temp,       "+ more" )
    } else  {  # Removes the last \n from the string for better formatting into the .csv file
        listMarkerLabels.temp <- substr ( listMarkerLabels.temp, 1, nchar ( listMarkerLabels.temp ) - 1 )
        listCellLabels.temp   <- substr ( listCellLabels.temp,   1, nchar ( listCellLabels.temp )   - 1 )
        listPhenotypeID.temp  <- substr ( listPhenotypeID.temp,  1, nchar ( listPhenotypeID.temp )  - 1 )
        listCellID.temp       <- substr ( listCellID.temp,       1, nchar ( listCellID.temp )       - 1 )
    }

    if( AddPlusMore == TRUE ) {
        listMarkerLabels.temp <- paste ( listMarkerLabels.temp, "\n + more", sep="" )
        listCellLabels.temp   <- paste ( listCellLabels.temp,   "\n + more", sep="" )
        listPhenotypeID.temp  <- paste ( listPhenotypeID.temp,  "\n + more", sep="" )
        listCellID.temp       <- paste ( listCellID.temp,       "\n + more", sep="" )
    }
    return ( c ( listMarkerLabels.temp, listCellLabels.temp, listPhenotypeID.temp, listCellID.temp ) )
}
###############################################################3
# Update the markers_ShortName_OntologyName.csv file
UpdateMarkers_ShortName_OntologyName <- function(save.dir, marker.list.short, marker.list, q3, q1) {

    fname <-  paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
    temp.table <- read.csv ( fname, header = FALSE ) ;   temp.table <- as.matrix ( temp.table )
    temp.matrix <- matrix ( 0, length ( temp.table[,1] ) + 1, 2 )
    temp.matrix[1:length ( temp.table[,1] ) , ] <- temp.table
    temp.table <- temp.matrix
    temp.marker.short <- as.character ( marker.list.short[[q3]][q1] )

    temp.marker.short <- tempMarkerShort(temp.marker.short)

    temp.table[length ( temp.table[,1] ), 1] <- temp.marker.short
    temp.table[length ( temp.table[,1] ), 2] <- paste ( marker.list[[q3]][q1] )
    write.table(temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE)

}
###############################################################3
# List of phenotypes commonly used by HIPC
listPhenotypes_flowCL <- function ( ) {
return <- c(
"CD3+","CD4+","CD8+","CCR7+","CD45RA+","CD38+","HLA-DR+","CD127+","CD25+","CCR4+","CD45RO+","CXCR3+",
"CCR6+","CD19+","CD20+","CD27+","IgD+","CD24+","CD14+","CD11c+","CD123+","CD16+","CD56+",
"CD3-","CD4-","CD8-","CCR7-","CD45RA-","CD38-","HLA-DR-","CD127-","CD25-","CCR4-","CD45RO-","CXCR3-",
"CCR6-","CD19-","CD20-","CD27-","IgD-","CD24-","CD14-","CD11c-","CD123-","CD16-","CD56-",
"CD127lo","CD24hi","CD38hi","CD27hi","CD56hi","CD56lo",
"CD3+CD4+CD8-CD38+HLA-DR+",
"CD3+CD4+CD8-CCR7+CD45RA+",
"CD3+CD4+CD8-CCR7+CD45RA-",
"CD3+CD4+CD8-CCR7-CD45RA-",
"CD3+CD4+CD8-CCR7-CD45RA+",
"CD3+CD4-CD8+CD38+HLA-DR+",
"CD3+CD4-CD8+CCR7+CD45RA+",
"CD3+CD4-CD8+CCR7+CD45RA-",
"CD3+CD4-CD8+CCR7-CD45RA-",
"CD3+CD4-CD8+CCR7-CD45RA+",
"CD3+CD4+CD127loCD25+",
"CD3+CD4+CD127loCD25+CCR4+CD45RO-",
"CD3+CD4+CD127loCD25+CCR4+CD45RO+",
"CD3+CD4+CD127loCD25+CCR4+",
"CD3+CD4+CD127loCD25+CCR4+HLA-DR+",
"CD3-CD19+CD20+CD27-IgD+",
"CD3-CD19+CD20+CD27+IgD+",
"CD3-CD19+CD20+CD27+IgD-",
"CD3-CD19+CD20+CD24hiCD38hi",
"CD3-CD19+CD20-CD27hiCD38hi",
"CD3-CD19-CD20-CD14+",
"CD3-CD19-CD20-CD14+CD16+",
"CD3-CD19-CD20-CD14+CD16-",
"CD3-CD19-CD20-CD14-HLA-DR-CD16+CD56lo",
"CD3-CD19-CD20-CD14-HLA-DR-CD16-CD56hi",
"CD3-CD19-CD20-CD14-CD16-CD56-HLA-DR+",
"CD3-CD19-CD20-CD14-CD16-CD56-HLA-DR+CD11c-CD123+",
"CD3-CD19-CD20-CD14-CD16-CD56-HLA-DR+CD11c+CD123-",
"CD3+CD4+CD8-CD38+HLA-DR+",
"CD3+CD4+CD8-CXCR3+CCR6-",
"CD3+CD4+CD8-CXCR3-CCR6-",
"CD3+CD4+CD8-CXCR3-CCR6+",
"CD3+CD4-CD8+CD38+HLA-DR+",
"CD3+CD4-CD8+CXCR3+CCR6-",
"CD3+CD4-CD8+CXCR3-CCR6-",
"CD3+CD4-CD8+CXCR3-CCR6+")
}
####################################################
# List of colours used for the tree diagram
listColours_flowCL <- function ( ) {
return <- c(
"red", "blue", "forestgreen", "darkviolet", "gold2", "darkorange3", "aquamarine4", "aquamarine2", "khaki3", "deeppink", "tan4",
"darkmagenta", "coral2", "burlywood4", "hotpink", "lightblue2", "lightpink4", "mediumseagreen", "navyblue", "olivedrab1",
"orangered", "purple", "peru", "plum2", "tan3", "bisque1", "darkslategray3", "darkslategray1",
"gray23", "gray24", "gray25", "gray26", "gray27", "gray28", "gray29", "gray30", "gray31", "gray32", "gray33", "gray34",
"gray35", "gray36", "gray37", "gray38", "gray39", "gray40", "gray41", "gray42", "gray43", "gray44", "gray45", "gray46",
"gray47", "gray48", "gray49", "gray50", "gray51", "gray52", "gray53", "gray54", "gray55", "gray56", "gray57", "gray58",
"gray59", "gray60", "gray61", "gray62", "gray63", "gray64", "gray65", "gray66", "gray67", "gray68", "gray69", "gray70",
"gray71", "gray72", "gray73", "gray74", "gray75", "gray76", "gray77", "gray78")
}
##################################
# function for loading hasProperLabel data
flowCL_query_data_hasProperLabel <- function(){
return <- c("select distinct ?x ?label",
"{",
"?x a owl:Class.",
"?x rdfs:label ?label. ",
"FILTER regex(?label, \"$marker\", \"i\")",
"}")
}
##################################
# function for loading prefix.info data
flowCL_query_data_prefix.info <- function(){

return <- c("# Common prefix and abbreviation",
"PREFIX sesame: <http://www.openrdf.org/schema/sesame#>",
"prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>",
"prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>",
"prefix owl: <http://www.w3.org/2002/07/owl#>",
"prefix obo: <http://purl.obolibrary.org/obo/>",
"prefix oboinowl: <http://www.geneontology.org/formats/oboInOwl#>",
"prefix probe: <http://purl.obolibrary.org/obo/CL_0000903>",
"# Has plasma membrane part (has_pmp)",
"prefix has_pmp: <http://purl.obolibrary.org/obo/RO_0002104>",
"prefix lacks_pmp: <http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>",
"prefix has_high_pma: <http://purl.obolibrary.org/obo/cl#has_high_plasma_membrane_amount>",
"prefix has_low_pma: <http://purl.obolibrary.org/obo/cl#has_low_plasma_membrane_amount>",
"prefix definition: <http://purl.obolibrary.org/obo/IAO_0000115>")
}
##################################
# function for loading hasProperSynonym data
flowCL_query_data_hasProperSynonym <- function(){
return <- c("select distinct ?x ?label ?synonym ",
"where",
"{",
"?x a owl:Class.",
"?x rdfs:label ?label.",
"?x oboinowl:hasExactSynonym ?synonym. ",
"FILTER regex(?synonym, \"$marker\", \"i\")",
"}")
}
##################################
# function for loading getParentClasses data
flowCL_query_data_getParentClasses <- function(){
return <- c("# Find all parent classes of the cell type of interest. Note that for some reason,",
"# matching on ?x does not work, but matching on ?celllabel (x's label) does.",
"# Matching directly on ?x works on http://sparql.hegroup.org/sparql !",
# "select distinct ?x ?celllabel ?parent ?pl", # Jonathan's and Alan's fix to remove the @en by converting to a string
"select distinct ?x ?celllabel ?parent ?parentlabel",
"where",
"{",
# "  BIND (STR(?parentlabel) AS ?pl )", # Jonathan's and Alan's fix to remove the @en by converting to a string
"  ?parent a owl:Class.",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?parent.",
"  ?parent rdfs:label ?parentlabel.",
"  FILTER regex(?celllabel, \"$label\", \"i\")",
"}")
}
##################################
# function for loading hasPlasmaMembranePart data
flowCL_query_data_hasPlasmaMembranePart <- function(){
return <- c(
"select distinct ?x ?celllabel ?pl ?marker ?markerlabel",
# "select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  BIND (STR(?plabel) AS ?pl )", # Jonathan's and Alan's fix to remove the @en by converting to a string
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty has_pmp:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel.",
"  has_pmp: rdfs:label ?plabel.",
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",
"}")
}
##################################
# function for loading hasLowPlasmaMembraneAmount data
flowCL_query_data_hasLowPlasmaMembraneAmount <- function(){
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty has_low_pma:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel.  ",
"  has_low_pma: rdfs:label ?plabel.",
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",
"}")
}
##################################
# function for loading hasHighPlasmaMembraneAmount data
flowCL_query_data_hasHighPlasmaMembraneAmount <- function(){
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty has_high_pma:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel.  ",
"  has_high_pma: rdfs:label ?plabel.",
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",
"}")
}
##################################
# function for loading lacksPlasmaMembranePart data
flowCL_query_data_lacksPlasmaMembranePart <- function(){
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty lacks_pmp:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel. ",
"  lacks_pmp: rdfs:label ?plabel.",
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",
"}")
}

# ##################################
# # function for loading hasPMPsingle data
# flowCL_query_data_hasPMPsingle <- function(){
# return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
# # return <- c("select distinct ?x ?celllabel ?pl ?marker ?markerlabel",
# "where",
# "{",
# # "  BIND (STR(?plabel) AS ?pl )", # Jonathan's and Alan's fix to remove the @en by converting to a string
# "  ?x a owl:Class.",
# "  ?x rdfs:label ?celllabel.",
# "  ?x rdfs:subClassOf ?sub.",
# "  ?sub rdf:type owl:Restriction.",
# "  ?sub owl:onProperty has_pmp:.",
# "  ?sub owl:someValuesFrom ?marker.",
# "  ?marker rdfs:label ?markerlabel. ",
# "  has_pmp: rdfs:label ?plabel.",
# "  FILTER regex(?markerlabel, \"$marker\", \"i\")",
# "  FILTER regex(?celllabel, \"$celllabel\", \"i\")",
# "}")
# }
# ##################################
# # function for loading lacksPMPsingle data
# flowCL_query_data_lacksPMPsingle <- function(){
# return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
# "where",
# "{",
# "  ?x a owl:Class.",
# "  ?x rdfs:label ?celllabel.",
# "  ?x rdfs:subClassOf ?sub.",
# "  ?sub rdf:type owl:Restriction.",
# "  ?sub owl:onProperty lacks_pmp:.",
# "  ?sub owl:someValuesFrom ?marker.",
# "  ?marker rdfs:label ?markerlabel. ",
# "  lacks_pmp: rdfs:label ?plabel.",
# "  FILTER regex(?markerlabel, \"$marker\", \"i\")",
# "  FILTER regex(?celllabel, \"$celllabel\", \"i\")",
# "}")
# }

##################################
# function for loading the date
flowCL_query_date <- function(){
return <- c("PREFIX :<http://purl.obolibrary.org/obo/cl.owl#>",
"PREFIX geo-pos:<http://www.w3.org/2003/01/geo/wgs84_pos#>",
"PREFIX uberon:<http://purl.obolibrary.org/obo/uberon#>",
"PREFIX umbel-ac:<http://umbel.org/umbel/ac/>",
"PREFIX sw-vocab:<http://www.w3.org/2003/06/sw-vocab-status/ns#>",
"PREFIX ff:<http://factforge.net/>",
"PREFIX music-ont:<http://purl.org/ontology/mo/>",
"PREFIX dc-term:<http://purl.org/dc/terms/>",
"PREFIX om:<http://www.ontotext.com/owlim/>",
"PREFIX opencyc-en:<http://sw.opencyc.org/2008/06/10/concept/en/>",
"PREFIX factbook:<http://www.daml.org/2001/12/factbook/factbook-ont#>",
"PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>",
"PREFIX pext:<http://proton.semanticweb.org/protonext#>",
"PREFIX ot:<http://www.ontotext.com/>",
"PREFIX dc:<http://purl.org/dc/elements/1.1/>",
"PREFIX onto:<http://www.ontotext.com/>",
"PREFIX foaf:<http://xmlns.com/foaf/0.1/>",
"PREFIX yago:<http://mpii.de/yago/resource/>",
"PREFIX umbel:<http://umbel.org/umbel#>",
"PREFIX pkm:<http://proton.semanticweb.org/protonkm#>",
"PREFIX wordnet16:<http://xmlns.com/wordnet/1.6/>",
"PREFIX owl:<http://www.w3.org/2002/07/owl#>",
"PREFIX gr:<http://purl.org/goodrelations/v1#>",
"PREFIX wordnet:<http://www.w3.org/2006/03/wn/wn20/instances/>",
"PREFIX opencyc:<http://sw.opencyc.org/concept/>",
"PREFIX wordn-sc:<http://www.w3.org/2006/03/wn/wn20/schema/>",
"PREFIX nytimes:<http://data.nytimes.com/>",
"PREFIX dbp-prop:<http://dbpedia.org/property/>",
"PREFIX geonames:<http://sws.geonames.org/>",
"PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>",
"PREFIX dbpedia:<http://dbpedia.org/resource/>",
"PREFIX oasis:<http://psi.oasis-open.org/iso/639/#>",
"PREFIX geo-ont:<http://www.geonames.org/ontology#>",
"PREFIX umbel-en:<http://umbel.org/umbel/ne/wikipedia/>",
"PREFIX ubprop:<http://purl.obolibrary.org/obo/ubprop#>",
"PREFIX bbc-pont:<http://purl.org/ontology/po/>",
"PREFIX ptop:<http://proton.semanticweb.org/protontop#>",
"PREFIX lingvoj:<http://www.lingvoj.org/ontology#>",
"PREFIX fb:<http://rdf.freebase.com/ns/>",
"PREFIX dbtune:<http://dbtune.org/bbc/peel/work/>",
"PREFIX obo:<http://purl.obolibrary.org/obo/>",
"PREFIX psys:<http://proton.semanticweb.org/protonsys#>",
"PREFIX umbel-sc:<http://umbel.org/umbel/sc/>",
"PREFIX dbp-ont:<http://dbpedia.org/ontology/>",
"PREFIX xsd:<http://www.w3.org/2001/XMLSchema#>",
"PREFIX ub:<http://www.lehigh.edu/~zhp2/2004/0401/univ-bench.owl#>",
"PREFIX oboInOwl:<http://www.geneontology.org/formats/oboInOwl#>",
"PREFIX skos:<http://www.w3.org/2004/02/skos/core#>",

"select ?s  ?p",
"where",
"{",
"?s owl:versionIRI ?p.",
"}"
)
}
##################################
# function for querying celllabel's lacks PMP
flowCL_query_data_celllabel_lacksPMP <- function(){
return <- c(
"select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty lacks_pmp:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel. ",
"  lacks_pmp: rdfs:label ?plabel.",
"  FILTER regex(?celllabel, \"$celllabel\", \"i\")",
"}"
)
}
##################################
# function for querying celllabel's has PMP
flowCL_query_data_celllabel_hasPMP <- function(){
return <- c(
# "select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"select distinct ?x ?celllabel ?pl ?marker ?markerlabel",
"where",
"{",
"  BIND (STR(?plabel) AS ?pl )", # Jonathan's and Alan's fix to remove the @en by converting to a string
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty has_pmp:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel. ",
"  has_pmp: rdfs:label ?plabel.",
"  FILTER regex(?celllabel, \"$celllabel\", \"i\")",
"}"
)
}
##################################
# function for querying celllabel's low PMA
flowCL_query_data_celllabel_lowPMA <- function(){
return <- c(
"select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty has_low_pma:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel. ",
"  has_low_pma: rdfs:label ?plabel.",
"  FILTER regex(?celllabel, \"$celllabel\", \"i\")",
"}"
)
}
##################################
# function for querying celllabel's high PMA
flowCL_query_data_celllabel_highPMA <- function(){
return <- c(
"select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",
"{",
"  ?x a owl:Class.",
"  ?x rdfs:label ?celllabel.",
"  ?x rdfs:subClassOf ?sub.",
"  ?sub rdf:type owl:Restriction.",
"  ?sub owl:onProperty has_high_pma:.",
"  ?sub owl:someValuesFrom ?marker.",
"  ?marker rdfs:label ?markerlabel. ",
"  has_high_pma: rdfs:label ?plabel.",
"  FILTER regex(?celllabel, \"$celllabel\", \"i\")",
"}"
)
}
##################################
# unit test function to test if the server is connected and ready to be queried
test.flowCL.connection <- function()
{
    # require("RUnit")
    testsuite <- RUnit::defineTestSuite("flowCL.check",
                    dirs = system.file("unitTests", package="flowCL"),
                    testFileRegexp = "^test_.*\\.R$", testFuncRegexp = "^test.+")
    testResult <- RUnit::runTestSuite(testsuite)
    RUnit::printTextProtocol(testResult)
}
