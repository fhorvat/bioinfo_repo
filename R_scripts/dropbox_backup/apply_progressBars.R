# by 
by_pb <- function(TheTitle=NULL,Data,INDICES,FUN, ...)
{
    env <- environment()
    pb_Total <- length(unique(INDICES))
    counter <- 0
    # pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    pb <- winProgressBar(title=TheTitle, min = 0, max = pb_Total, width = 300)
    
    wrapper <- function(...){
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        # setTxtProgressBar(get("pb", envir=env), curVal +1)
        setWinProgressBar(get("pb", envir= env),curVal +1)
        FUN(...)
    }
    res <- by(Data,INDICES,wrapper,...)
    close(pb)
    res
}

# STATUS: WORKING, but only tested once or twice,
# tested with most ?apply examples
# ISSUES/TODO: MARGIN argument cannot take a
# vector like 1:2 that is more than one numeric
apply_pb <- function(X, MARGIN, FUN, ...)
{
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir= env)
        setTxtProgressBar(get("pb", envir= env),
                          curVal +1)
        FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
}

## NOT RUN:
# apply_pb(anscombe, 2, sd, na.rm=TRUE)

## large dataset
# df <- data.frame(rnorm(30000), rnorm(30000))
# apply_pb(df, 1, sd)

###############################################################

lapply_pb <- function(X, FUN, ...)
{
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   

    # wrapper around FUN
    wrapper <- function(...){
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        setTxtProgressBar(get("pb", envir=env), curVal +1)
        FUN(...)
    }
    res <- lapply(X, wrapper, ...)
    close(pb)
    res
}

## NOT RUN:
# l <- sapply(1:20000, function(x) list(rnorm(1000)))
# lapply_pb(l, mean)

###############################################################

sapply_pb <- function(X, FUN, ...)
{
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    
    wrapper <- function(...){
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        setTxtProgressBar(get("pb", envir=env), curVal +1)
        FUN(...)
    }
    res <- sapply(X, wrapper, ...)
    close(pb)
    res
}

## NOT RUN:
# l <- sapply(1:20000, function(x) list(rnorm(1000))
# sapply_pb(l, mean)
