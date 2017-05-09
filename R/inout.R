# Is point in polygon?
#ray-casting algorithm True if in false if out
inout=function(s,vertices){
  count=0
  for(i in 1:(nrow(vertices)-1)){
    if(intersect(s,vertices[i,],vertices[i+1,])){
      count=count+1
    }
  }
  if(count %% 2 != 0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

intersect=function(s,vertex1,vertex2){
  if(s[2]==vertex1[2]|s[2]==vertex2[2]){
    s[2]=s[2]+ 0.000001
  }
  #vertex1 must be below vertex2
  if(vertex1[2]>vertex2[2]){
    swap=vertex1
    vertex1=vertex2
    vertex2=swap
  }
  if(s[2]<vertex1[2]|s[2]>vertex2[2]){
    out=FALSE
  }else if(s[1] > max(vertex1[1], vertex2[1]) ){
    out=FALSE
  }else{
    if(s[1] < min(vertex1[1], vertex2[1])){
      out=TRUE
    }else{
      if(vertex1[1]!=vertex2[1]){
        m_red=(vertex2[2]-vertex1[2])/(vertex2[1]-vertex1[1])
      }else{
        m_red=Inf
      }
      if(vertex1[1]!=s[1]){
        m_blue=(s[2]-vertex1[2])/(s[1]-vertex1[1])
      }else{
        m_blue=Inf
      }
      if(m_blue>=m_red){
        out=TRUE
      }else{
        out=FALSE
      }
    }
  }
  return(out)
}
