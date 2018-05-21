





aa=as.data.frame(id=b)
bb=as.data.frame(cbind(id = seq(1:4),t(h)))
g = join(aa, bb)


a = NaN
a[is.nan(a)] = 0



system.time(discrete_convolution(aaa,b))

discrete_convolution <- function(VectorSmall,VectorBig){
  lengthSmall = length(VectorSmall)
  lengthBig = length(VectorBig)
  Tem = matrix(0.0,(lengthBig + lengthSmall -1), lengthSmall)
  for (i in 1:lengthSmall) {
    Tem[i:(lengthBig-1+i),i] = VectorBig * VectorSmall[i]
  }
  Result = rowSums(Tem)
  return(Result)
}

b = rep(a,3)

outer(a,b,FUN = fctlength)
table(unlist(a))

a=c(1,2,3)
b = rep(a,3)







