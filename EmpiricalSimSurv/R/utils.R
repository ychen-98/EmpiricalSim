# Null-coalescing operator: return x if non-NULL, else y
`%||%` <- function(x, y) if (!is.null(x)) x else y
