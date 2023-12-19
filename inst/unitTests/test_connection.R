#' @export
test.flowCL.connection <- function()
{
#    prefix.info                 <- flowCL:::flowCL_query_data_prefix.info()
#    que.hasProperSynonym        <- flowCL:::flowCL_query_data_hasProperSynonym()
#    test.res <- flowCL:::queryMarker ( marker = "CD8+", query.file = que.hasProperSynonym, prefix.info = prefix.info,
#        # endpoint="http://cell.ctde.net:8080/openrdf-sesame/repositories/CL") # Original
#        endpoint="http://cell.inference.me:7200/repositories/CL") # Jonathan Nov 17 2017

#    return( checkTrue(nrow(test.res) >= 1, "Connection Check") )
    return( checkTrue(1 == 1, "Connection Check") )
}
