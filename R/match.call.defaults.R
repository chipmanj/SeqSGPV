#' @title Get list of inputs
#'
#' @export
match.call.defaults <- function(...) {

  # Return inputs of call, including default inputs
  # https://stackoverflow.com/questions/14397364/match-call-with-default-arguments
  # Answer from Neal Fultz

  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))

  for(i in setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )


  match.call(sys.function(sys.parent()), call)
}
