###
## helper functions for spatial analyses, such as getting a matrix of distances
## or a tidygraph of neighbors
###

dist_mat <- function(df, normalize=FALSE, cutoff=NULL) {
    require(geosphere)
    
    dm <- df %>%
        select(long, lat) %>%
        distm()
    
    colnames(dm) <- df$id
    rownames(dm) <- df$id

    if (normalize) dm <- dm / max(dm)
    if (!is.null(cutoff)) {
        dm[dm <= cutoff] <- 1
        dm[dm > cutoff] <- 0
        diag(dm) <- 0 # no self-loops
    }
    return(dm)
}
        
knn_mat <- function(dm, k=1) {
    diag(dm) <- max(dm) + 1 # so this def isnt the min
    ret = dm
    for (r in seq_len(nrow(dm))) {
        nc = 0
        while (nc < k) {
            ## for house/row r, find and mark smallest distance in ret, and "remove"
            ## from dm so can find next smallest
            ret[r, which.min(dm[r, ])] <- -1
            dm[r, which.min(dm[r, ])] <- max(dm) + 1
            nc <- nc + 1
        }
    }
    ## make adj matrix from marked pairs
    ret[ret >= 0] <- 0
    ret[ret < 0] <- 1
    return(ret)
}

make_cov_net <- function(df_org, method=c('knn', 'knn_rev', 'cutoff'), ...) {
    require(tidygraph)
    require(tidyverse)
    require(geosphere)

    ## get adj matrix according to specified method
    if (method == 'knn') { # "my cov effect my nns"
        dm <- dist_mat(df_org)
        am <- knn_mat(dm, list(...)$k)
    }
    if (method == 'knn_rev') { # "my cov effect those who have me as a nn"
        dm <- dist_mat(df_org)
        am <- t(knn_mat(dm, list(...)$k))
    }
    if (method == 'cutoff') {
        am <- dist_mat(df_org, cutoff=list(...)$h)
    }
    g <- as_tbl_graph(am, directed=TRUE) %>%
        mutate(
            long = df_org$long,
            lat = df_org$lat,
            inf = as.factor(df_org$infestation),
            !!!df_org[, 5:ncol(df_org)]
        ) %>%
        activate(edges)

    from_idx <- pull(g, from)
    to_idx <- pull(g, to)
    edge_attr <- df_org[from_idx, 5:ncol(df_org)-1]

    mutate(g, !!!edge_attr) %>%
        mutate(target_inf = as.factor(df_org$infestation[to_idx])) %>%
        activate(nodes)
}

add_neighbor_attr <- function(tbl_graph, attr, new_col, fun, ...) {
    attr_col <- enquo(attr)
    ## grab original attr for each node
    node_attr <- tbl_graph %N>%
        as_tibble %>%
        pull(!!attr_col)

    ## find neighbors of each node and apply fun
    tbl_graph %N>%
        mutate(
            !!new_col := map_dbl(
                  local_members(mindist=1, mode='in'),
                  ~fun(node_attr[.x]),
                  ...
              )
        )
}

