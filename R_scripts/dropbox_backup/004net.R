
comp <- function(a, b) {
  if (a == b) { return(0) }
  if (a > b) { return(1) }
  return(2)
}

tree.add <- function(tree, value) {
  new.node <- function(value) {
    node <- list(high = NULL, low = NULL)
    attr(node, "value") <- value
    return(node)
  }
  if (is.null(tree)) { return(new.node(value)) }
  sub.tree <- tree
  path <- numeric(0)
  while (TRUE) {
    cr <- comp(value, attr(sub.tree, "value"))
    if (cr == 0) {
      return(tree)
    }
    path <- c(path, cr)
    sub.tree <- sub.tree[[cr]]
    if (is.null(sub.tree)) {
      tree[[path]] <- new.node(value)
      return(tree)
    }
  }
}



