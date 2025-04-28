# Luke Keele's utils
# Create a data frame with squares and/or interactions
make.sq_inter = function(df, is.square=T, is.inter=T, keep.marginal=T) {
  
  if(keep.marginal)
  {
    df.appended = df
  }
  
  if(is.inter)
  {
    df.inter = matrix(rep(0, nrow(df)* choose(ncol(df),2)), nrow =nrow(df))
    
    # introduce a pseudo outcome and run a linear model
    y.pseudo = df[,1]
    df.with.outcome = data.frame(y.pseudo = y.pseudo, df)
    
    df.inter = data.frame(model.matrix(y.pseudo~.*.-1, data = df.with.outcome))
    names.inter = colnames(df.inter)
    
    if(keep.marginal)
    {
      df.appended = df.inter
    }
    if(!keep.marginal)
    {
      df.appended = as.data.frame(df.inter[,-(1:ncol(df))])
      colnames(df.appended) = names.inter[-(1:ncol(df))]
    }
    
  }
  
  if(is.square)
  {
    df.square  = df 
    for(j in 1:ncol(df))
    {
      df.square[,j] = (df[,j])^2
      colnames(df.square)[j] = paste(colnames(df)[j], "^2")
    }
    
    if(is.inter)
    {
      df.appended = data.frame(df.appended, df.square)
    }
    if(!is.inter)
    {
      if(keep.marginal)
      {
        df.appended = data.frame(df.appended, df.square)
      }
      if(!keep.marginal)
      {
        df.appended = data.frame(df.square)
      }
    }
  }
  return(as.data.frame(df.appended))
}