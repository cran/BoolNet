
# Split up <expression> into symbols.
# Returns a vector of symbols.
scan <- function(expression)
{
  # add extra whitespace to brackets and operators
  expression <- gsub("("," ( ", expression, fixed=TRUE)
  expression <- gsub(")"," ) ", expression, fixed=TRUE)
  expression <- gsub("!"," ! ", expression, fixed=TRUE)  
  expression <- gsub("&"," & ", expression, fixed=TRUE)
  expression <- gsub("|"," | ", expression, fixed=TRUE)    

  # strip multiple whitespace characters
  expression <- gsub("[ ]+", " ", expression)
  expression <- gsub("^[ ]+", "", expression)
  expression <- gsub("[ ]+$", "", expression)
  
  # split up at whitespace positions
  res <- strsplit(expression, " ", fixed=TRUE)[[1]]
  return(res)
}

# Parse a Boolean function in <expression>,
# and build a corresponding parse tree.
parseBooleanFunction <- function(expression)
{
  
  # internal function to step forward one symbol
  advanceSymbol <- function()
  {
    pos <<- pos + 1
    if (pos > length(symbols))
      return(NA)
    else
      return(symbols[pos])
  }
  
  # internal function to parse a sub-expression in brackets
  parseExpression <- function()
  {
    operators <- c()
    children <- c()
    symb <- advanceSymbol()
    
    while (TRUE)
    {
      if (symb == "(")
        # a new sub-expression
        children[[length(children)+1]] <- parseExpression()
      else
      if (symb == "!")
        # a negation
        children[[length(children)+1]] <- parseNegation()
      else
      if (symb %in% c(")","&","|"))
        # something unexpected happened
        stop(paste("Unexpected symbol:",symb))
      else
        # an atom (variable or constant)
        children[[length(children)+1]] <- list(type="atom",name=symb,
                                               negated=FALSE)
      
      symb <- advanceSymbol()  
      
      if (is.na(symb) || symb == ")")
        # end of expression was reached      
        break
        
      if (!symb %in% c("&","|"))
        stop("Operator expected!")
      
      operators <- c(operators,symb)
      
      symb <- advanceSymbol()
    }
    
    if (length(children) == 1)
      # an operator with only one child
      return(children[[1]])
    else
    if (length(unique(operators)) == 1)
      # an n-ary operator
      return(list(type="operator",operator=operators[1],
                  negated=FALSE, operands=children))
    else
    # a mixture of multiple operators => ensure correct precedence
    {
      i <- 1
      startPos <- NA
      operators <- c(operators, "|")
      operands <- list()

      while (i <= length(operators))
      # identify AND operators and move them to a subtree
      {
        if (operators[i] == "&" && is.na(startPos))
          # start of a conjunction
          startPos <- i
        if (operators[i] == "|" && !is.na(startPos))
        # end of a conjunction => create subtree
        {
          subOp <- list(type="operator", operator="&",
                        negated=FALSE, operands=children[startPos:i])
          operands[[length(operands)+1]] <- subOp
          startPos <- NA
        }
        else
        if (is.na(startPos))
          operands[[length(operands)+1]] <- children[[i]]
        i <- i + 1  
      }
      return(list(type="operator",operator="|",
                  negated=FALSE, operands=operands))
    }
  }
  
  # Internal function to parse a negation
  parseNegation <- function()
  {
    symb <- advanceSymbol()
    if (symb == "(")
    {
      res <- parseExpression()
    }
    else
    if (symb %in% c(")","&","|"))
      stop(paste("Unexpected symbol:",symb))
    else
      res <- list(type="atom",name=symb)
    res$negated = TRUE  
    return(res)
  }

  symbols <- scan(expression)
  
  pos <- 0
  
  return(parseExpression())
}

# Currently unused - 
# Regenerate the expression string 
# from the parse tree <tree>
stringFromParseTree <- function(tree)
{
  res <- switch(tree$type,
    operator = 
    {
      paste({if (tree$negated) "!" else ""},
            "(",
            paste(sapply(tree$operands,stringFromParseTree), 
                  collapse=paste(" ",tree$operator," ",sep="")),
            ")", sep="")
    },    
    atom = paste({if (tree$negated) "!" else ""},
                 tree$name,
                 sep=""))
  return(res)  
}
