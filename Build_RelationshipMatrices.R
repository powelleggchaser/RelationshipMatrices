#Chapter 13

#### Exercise 1 #####
 
#Build A

# Own funtion to build Numerator Relationship Matrices (A matrix)

#Build A and A inverse through recursive method, and complete inversion
#Inverse not feasible with large scale data
BuildA<-function(pedigree,inverse=TRUE){
  
  A = matrix(nrow=length(pedigree$ID),ncol=length(pedigree$ID))
  for (i in 1:nrow(pedigree)){
    i_dam = as.numeric(as.character(pedigree$dam[i]))
    i_sire = as.numeric(as.character(pedigree$sire[i]))
    j = as.numeric(as.character(pedigree$ID[pedigree$ID!=i]))
    if (i_sire==0&i_dam==0){
      A[i,i]=1
      for(each in j){
        A[i,j[each]]=0; A[j[each],i]=0
      }
    }
    else if (i_sire!=0&i_dam==0){
      A[i,i]=1
      for(each in j){
        A[i,j[each]]=0.5*(A[j[each],i_sire]);A[j[each],i]=A[i,j[each]]
      }
    }
    else if (i_sire==0&i_dam!=0){
      A[i,i]=1
      for(each in j){
        A[i,j[each]]=(0.5*(A[j[each],i_dam]));A[j[each],i]=A[i,j[each]]
      }
    }
    else if (i_sire!=0&i_dam!=0){
      A[i,i]= 1 + (0.5*(A[i_sire,i_dam]))
      for(each in j){
        A[i,j[each]]=(0.5*((A[j[each],i_sire])+(A[j[each],i_dam]))); A[j[each],i]=A[i,j[each]]
      }
    }
  }
  if ((inverse==TRUE)|(inverse==T)){
    if (det(A)==0){
      print("Matrix is Singular. Non-Invertible")
    }
    else{
      Ainv<-solve(A)
      Ainv[is.na(Ainv)]=0
      return(Ainv)
    }
  }
  else{
    A[is.na(A)]=0
    return(A)
  }
}

BuildAinv_NoInb<-function(pedigree,Tinv,D_matrix,method="T"){
  if (method=="T"){
    TinvP=t(Tinv)
    D_inv=solve(D_matrix)
    Ainv=TinvP%*%Dinv%*%Tinv
    return(Ainv)
  }
  if (method=="NoT"){
    M = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
    for (i in 1:nrow(pedigree)){
      i_dam = as.numeric(as.character(pedigree$dam[i]))
      i_sire = as.numeric(as.character(pedigree$sire[i]))
      if (i_sire==0&i_dam==0){
        M[i,i_sire]=0
        M[i,i_dam]=0
      }
      if (i_sire!=0&i_dam==0){
        M[i,i_sire]=0.5
        M[i,i_dam]=0
      }
      if (i_sire==0&i_dam!=0){
        M[i,i_sire]=0
        M[i,i_dam]=0.5
      }
      if (i_sire!=0&i_dam!=0){
        M[i,i_sire]=0.5
        M[i,i_dam]=0.5
      }
    }
    I = diag(length(pedigree$ID))
    first = I-t(M)
    Dinv=solve(D_matrix)
    second = (I-M)
    Ainv = first%*%Dinv%*%second
    return(Ainv)
  }
}

#Build M matrix - Traces flow of genes from parents to offspring (e.g. 0.5 in each of the prent progent off diagonals)
#M_squared = traces flow of grandparent genes to individual 
#M power t - traces flow of individuals t generations ago to individual
#n_generation>1 may not work
BuildM <- function(pedigree,n_generations=1){
  M = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
  for (i in 1:nrow(pedigree)){
    i_dam = as.numeric(as.character(pedigree$dam[i]))
    i_sire = as.numeric(as.character(pedigree$sire[i]))
    if (i_sire==0&i_dam==0){
      M[i,i_sire]=0
      M[i,i_dam]=0
    }
    if (i_sire!=0&i_dam==0){
      M[i,i_sire]=0.5
      M[i,i_dam]=0
    }
    if (i_sire==0&i_dam!=0){
      M[i,i_sire]=0
      M[i,i_dam]=0.5
    }
    if (i_sire!=0&i_dam!=0){
      M[i,i_sire]=0.5
      M[i,i_dam]=0.5
    }
  }
  if (n_generations>1){
    M_ancestors = M^n_generations
    return(M_ancestors)
  }
  else {
    return(M)
  }
}

#Use A, T & D matrices to find A inverse 
BuildT<-function(pedigree,inverse=TRUE){
  
  T_mat = matrix(nrow=length(pedigree$ID),ncol=length(pedigree$ID))
  for (i in 1:nrow(pedigree)){
    i_dam = as.numeric(as.character(pedigree$dam[i]))
    i_sire = as.numeric(as.character(pedigree$sire[i]))
    j = as.numeric(as.character(pedigree$ID[pedigree$ID!=i]))
    T_mat[i,i]=1
    if (i_sire==0&i_dam==0){
      for(each in j){
        T_mat[i,j[each]]=0; T_mat[j[each],i]=0
      }
    }
    else if (i_sire!=0&i_dam==0){
      for(each in j){
        T_mat[i,j[each]]=0.5*(T_mat[i_sire,j[each]]); T_mat[j[each],i]=0
      }
    }
    else if (i_sire==0&i_dam!=0){
      for(each in j){
        T_mat[i,j[each]]=(0.5*(T_mat[i_dam,j[each]])); T_mat[j[each],i]=0
      }
    }
    else if (i_sire!=0&i_dam!=0){
      for(each in j){
        T_mat[i,j[each]]=(0.5*((T_mat[i_sire,j[each]])+(T_mat[i_dam,j[each]]))); T_mat[j[each],i]=0
      }
    }
  }
  if ((inverse==TRUE)|(inverse==T)){
    if (det(T_mat)==0){
      print("Matrix is Singular. Non-Invertible")
    }
    else{
      Tinv<-solve(T_mat)
      Tinv[is.na(Tinv)]=0
      return(Tinv)
    }
  }
  if (inverse=="IM"){
    M = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
    for (i in 1:nrow(pedigree)){
      i_dam = as.numeric(as.character(pedigree$dam[i]))
      i_sire = as.numeric(as.character(pedigree$sire[i]))
      if (i_sire==0&i_dam==0){
        M[i,i_sire]=0
        M[i,i_dam]=0
      }
      if (i_sire!=0&i_dam==0){
        M[i,i_sire]=0.5
        M[i,i_dam]=0
      }
      if (i_sire==0&i_dam!=0){
        M[i,i_sire]=0
        M[i,i_dam]=0.5
      }
      if (i_sire!=0&i_dam!=0){
        M[i,i_sire]=0.5
        M[i,i_dam]=0.5
      }
    }
    I = diag(length(pedigree$ID))
    Tinv = I - M
    Tinv[is.na(Tinv)]=0
    return(Tinv)
  }
  else{
    T_mat[is.na(T_mat)]=0
    return(T_mat)
    
  }
}
BuildD<-function(pedigree,A_matrix,inverse=TRUE){
  
  D = matrix(nrow=length(pedigree$ID),ncol=length(pedigree$ID))
  for (i in 1:nrow(pedigree)){
    i_dam = as.numeric(as.character(pedigree$dam[i]))
    i_sire = as.numeric(as.character(pedigree$sire[i]))
    i_paternal_grandsire = as.numeric(as.character(pedigree$sire[i_sire]))
    if (length(i_paternal_grandsire)==0){
      i_paternal_grandsire=0
      F_sire = 0
      F_dam = 0
    }
    i_paternal_granddam = as.numeric(as.character(pedigree$dam[i_sire]))
    if (length(i_paternal_granddam)==0){
      i_paternal_granddam=0
      F_sire = 0
      F_dam = 0
    }
    i_maternal_grandsire = as.numeric(as.character(pedigree$sire[i_dam]))
     if (length(i_maternal_grandsire)==0){
      i_maternal_grandsire=0
      F_sire = 0
      F_dam = 0
    }
    i_maternal_granddam = as.numeric(as.character(pedigree$dam[i_dam]))
    if (length(i_maternal_granddam)==0){
      i_maternal_granddam=0
      F_sire = 0
      F_dam = 0
    }
    if ((i_maternal_granddam!=0)&(i_maternal_grandsire!=0)){
      F_dam = (0.5*A_matrix[i_maternal_grandsire,i_maternal_granddam])
    }
    if ((i_paternal_granddam!=0)&(i_paternal_grandsire!=0)){
      F_sire = (0.5*A_matrix[i_paternal_grandsire,i_paternal_granddam])
    }
    if (i_sire==0&i_dam==0){
      D[i,i]=1
    }
    if (i_sire!=0&i_dam==0){
      D[i,i]=0.75 - (0.25*(F_sire))
    }
    if (i_sire==0&i_dam!=0){
      D[i,i]=0.75 - (0.25*(F_dam))
    }
    if (i_sire!=0&i_dam!=0){
      D[i,i]=0.5 - (0.25*(F_sire+F_dam))
    }
  }
  if ((inverse==TRUE)|(inverse==T)){
    Dinv = 1/D
    Dinv[is.na(Dinv)]=0
    return(Dinv)
  }
  else{
    D[is.na(D)]=0
    return(D)
  }
}

############### -------------- TEST FUNCTIONS -------------- ###############

#ID =c(1:10)
#dam = c(rep(0,5),2,4,4,7,7)
#sire = c(rep(0,5),1,3,5,6,8)
#sex = c("M","F","M","F","M","M","F","M","M","F")
#weight = c(rep(NA,5),400,300,312,405,298)
#pedigree = data.frame(cbind(ID,sire,dam,sex,weight))

ID = c(1:6)
sire = c(rep(0,2),1,1,4,5)
dam = c(rep(0,2),2,0,3,2)
pedigree = data.frame(cbind(ID,sire,dam))

#Build A
NRM = BuildA(pedigree,inverse=F)
#Buld inverse of A
Ainv = BuildA(pedigree,inverse=T)

#Build M
M = BuildM(pedigree,n_generations = 1)

#Build T
T_matrix = BuildT(pedigree,inverse=F)
#Buld inverse of T
Tinv = BuildT(pedigree,inverse=T)
Tinverse = BuildT(pedigree,inverse="IM")

#Build D
D = BuildD(pedigree,NRM,inverse = FALSE)
#Buld inverse of D
Dinv = BuildD(pedigree,NRM,inverse=T)

# Build inverse A using A=TDT' : With T
NI_Ainv=BuildAinv_NoInb(pedigree,T_matrix,D,method="T")

#Build inverse A using (I-M')Dinv(I-M) : Without T
NI_Ainv_NoT=BuildAinv_NoInb(pedigree,T_matrix,D,method="NoT")


BuildM(pedigree,n_generations = 8)

