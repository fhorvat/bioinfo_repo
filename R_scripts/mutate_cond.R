library(dplyr)

# mutate in rows for which cond is TRUE
mutate_cond <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  
  # Initialize any new variables as new_init
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data %>% filter(condition) %>% mutate(...)
  .data
}
