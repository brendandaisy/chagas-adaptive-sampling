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

dist_mat_sncar <- function(df, d_low = NULL) {
  require(geosphere)
  
  dm <- df %>%
    select(long, lat) %>%
    distm()
  
  colnames(dm) <- df$id
  rownames(dm) <- df$id
  diag(dm) <- 0
  
  if (is.null(d_low)) {
    return(dm / mean(dm))
  }
  
  dm / d_low
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

make_cov_net <- function(df_org, method=c('knn', 'knn_rev', 'cutoff', 'decay'), ...) {
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
  if (method == 'decay') {
    am <- exp(-.05 * dist_mat(df_org))
    diag(am) <- 0
    am <- am / max(am)
  }
  g <- as_tbl_graph(am, directed=TRUE) %>%
    mutate(
      long = df_org$long,
      lat = df_org$lat,
      !!!df_org[, 5:(ncol(df_org))]
    )
  
  # from_idx <- pull(g, from)
  # to_idx <- pull(g, to)
  # edge_attr <- df_org[from_idx, 5:ncol(df_org)]
  # 
  # mutate(g, !!!edge_attr) %>%
  #   mutate(target_inf = as.factor(df_org$infestation[to_idx])) %>%
  #   activate(nodes)
}

sum_factor <- function(f) {
  as.double(table(f))
}

pull_attr <- function(tbl_graph, attr, E=FALSE) {
  attr_col <- enquo(attr)
  if (E) {
    tbl_graph %E>%
      as_tibble %>%
      pull(!!attr_col)
  }
  else {
    tbl_graph %N>%
      as_tibble %>%
      pull(!!attr_col)
  }
}

## For summing/counting across edges
count_edge_attr <- function(tbl_graph, attr) {
  ## grab original attr for each node
  node_attr <- pull_attr(tbl_graph, attr)
  
  tbl_graph %N>%
    mutate(nb = local_members(mindist=1, mode='in')) %>%
    as_tibble %>%
    select(source_inf = infestation, nb) %>%
    pmap_dfr(~tibble(target_inf = .x, cov = node_attr[.y]))
}

count_neighbor_attr <- function(tbl_graph, attr) {
  nb_df <- tbl_graph %N>%
    mutate(
      nb_prof = map_local(
        mindist=1,
        mode='in',
        .f = function(neighborhood, node, graph) {
          pull_attr(neighborhood, attr)
        }
      )
    ) %>%
    as_tibble %>%
    select(name, infestation, nb_prof) %>%
    unnest(nb_prof)
  
  if (is.factor(nb_df$nb_prof)) {
    nb_df %<>%
      add_count(name, nb_prof) %>%
      distinct %>%
      pivot_wider(names_from=nb_prof, values_from=n) %>%
      mutate_at(2:ncol(.), ~ replace_na(.x, replace = 0))
    
    return(nb_df)
  }

  nb_df %>%
    group_by(name, infestation) %>%
    summarise(nb_prof = sum(nb_prof)) %>%
    ungroup
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

### For factors where you want each level done
add_neighbor_attr2 <- function(tbl_graph, attr, new_col, fun, ...) {
  ## grab original attr for each node
  node_attrs <- tbl_graph %N>%
    as_tibble %>%
    select(contains(attr))
  
  new_cols <- str_c(as.character(quote(fun)), colnames(node_attrs))
  
  ## find neighbors of each node and apply fun    
  tbl_graph %N>%
    bind_cols(
      map2_dfr(
        local_members(mindist=1, mode='in'),
        1:ncol(node_attr),
        ~tibble(!!enquo(new_cols[.y]) := fun(node_attrs[.x, .y])),
        ...
      )
    )
}
