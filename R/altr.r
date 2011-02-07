altr <- function(object, bu=burst(object)[1])
{
    names(object) <- unlist(lapply(object, function(x) attr(x, "burst")))
    attach(object[[bu]])
}

dltr <- function(object, bu=burst(object)[1])
{
    names(object) <- unlist(lapply(object, function(x) attr(x, "burst")))
    detach(object[[bu]])
}
