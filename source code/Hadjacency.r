#Hadjacency

Hadjacency<-function (datExpr, selectCols = NULL, type = "unsigned", power = if (type == 
    "distance") 1 else 6, corFnc = "cor", corOptions = "use = 'p'", 
    distFnc = "dist", distOptions = "method = 'euclidean'") 
{
    adjacencyTypes<-c('unsigned','signed','signed hybird','hybrid2','distance')
	intType = charmatch(type, adjacencyTypes)
    if (is.na(intType)) 
        stop(paste("Unrecognized 'type'. Recognized values are", 
            paste(adjacencyTypes, collapse = ", ")))
    corFnc.fnc = match.fun(corFnc)
    if (intType < 5) {
        if (is.null(selectCols)) {
            if (is.list(corOptions)) {
                cor_mat = do.call(corFnc.fnc, c(list(x = datExpr), 
                  corOptions))
            }
            else {
                corExpr = parse(text = paste(corFnc, "(datExpr ", 
                  prepComma(corOptions), ")"))
                cor_mat = eval(corExpr)
            }
        }
        else {
            if (is.list(corOptions)) {
                cor_mat = do.call(corFnc.fnc, c(list(x = datExpr, 
                  y = datExpr[, selectCols]), corOptions))
            }
            else {
                corExpr = parse(text = paste(corFnc, "(datExpr, datExpr[, selectCols] ", 
                  prepComma(corOptions), ")"))
                cor_mat = eval(corExpr)
            }
        }
    }
    else {
        if (!is.null(selectCols)) 
            stop("The argument 'selectCols' cannot be used for distance adjacency.")
        if (is.list(distOptions)) {
            d = do.call(distFnc, c(list(x = t(datExpr)), distOptions))
        }
        else {
            corExpr = parse(text = paste(distFnc, "(t(datExpr) ", 
                prepComma(distOptions), ")"))
            d = eval(corExpr)
        }
        if (any(d < 0)) 
            warning("Function WGCNA::adjacency: Distance function returned (some) negative values.")
        cor_mat = 1 - as.matrix((d/max(d, na.rm = TRUE))^2)
    }
    if (intType == 1) {
        cor_mat = abs(cor_mat)
    }
    else if (intType == 2) {
        cor_mat = (1 + cor_mat)/2
    }
    else if (intType == 3) {
        cor_mat[cor_mat < 0] = 0
    }
    else if(intType == 4) {
        cor_mat = (1+abs(cor_mat))/2
    }
    cor_mat^power



}
