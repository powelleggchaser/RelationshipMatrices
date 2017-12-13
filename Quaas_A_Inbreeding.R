#Quaas A matrix with Inbreeding
#Requires n*(n+1)/2 Operations; computational time ~ n squared; n = size of dataset(individuals)
#Previously calcualted L matrix needs stored if new data will be added
#Columnwise Calculation of L 

QuaasA = function(pedigree){

  L = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
  Alpha = rep(0,nrow(L))
  Ainv=matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
  #Dinv = rep(0,length(pedigree$ID))
  i = c(1:nrow(pedigree))
  #i = as.numeric(c((j+1):(max(j))))
  for (each in i){
    j = i[i!=each]
    i_sire = pedigree$sire[each]
    i_dam = pedigree$dam[each]
    sire_range = c(1:i_sire)
    dam_range = c(1:i_dam)
    L_sire = rep(0,i_sire)
    L_dam = rep(0,i_dam)
    if ((i_sire==0)&(i_dam==0)){
      L[each,each] = sqrt((1 -(0.25*(0 + 0))))
      Alpha[each]=1/(L[each,each]^2)
      Ainv[each,each] = Ainv[each,each] + Alpha[each]
    }
    if ((i_sire!=0)&(i_dam==0)){
      for (ind in sire_range){
        L_sire[ind] = (L[i_sire,ind]^2) 
      }
      L_sire = sum(L_sire)
      L[each,each] = sqrt((1 - (0.25*(L_sire))))
      Alpha[each]=1/(L[each,each]^2)
      j = i[i!=each]
      for (individual in j){
        L[each,individual] = 0.5*(L[i_sire,individual]+0)
      }
      Ainv[each,each] = Ainv[each,each] + Alpha[each]
      Ainv[i_sire,each] = Ainv[i_sire,each] - (Alpha[each]/2)
      Ainv[each,i_sire] = Ainv[each,i_sire] - (Alpha[each]/2)
      Ainv[i_sire,i_sire] = Ainv[i_sire,i_sire] + (Alpha[each]/4)
    }
    if ((i_sire==0)&(i_dam!=0)){
      for (indiv in dam_range){
        L_dam[indiv] = (L[i_dam,indiv]^2) 
      }
      L_dam = sum(L_dam)
      L[each,each] = sqrt((1 - (0.25*(L_dam))))
      Alpha[each]=1/(L[each,each]^2)
      j = i[i!=each]
      for (individual in j){
        L[each,individual] = 0.5*(L[0+i_dam,individual])
      }
      Ainv[each,each] = Ainv[each,each] + Alpha[each]
      Ainv[i_dam,each] = Ainv[i_dam,each] - (Alpha[each]/2)
      Ainv[each,i_dam] = Ainv[each,i_dam] - (Alpha[each]/2)
      Ainv[i_dam,i_dam] = Ainv[i_dam,i_dam] + (Alpha[each]/4)
    }
    if ((i_sire!=0)&(i_dam!=0)){
      for (ind in sire_range){
        L_sire[ind] = (L[i_sire,ind]^2) 
      }
      L_sire = sum(L_sire)
      for (indiv in dam_range){
        L_dam[indiv] = (L[i_dam,indiv]^2) 
      }
      L_dam = sum(L_dam)
      L[each,each] = sqrt((1 - (0.25*(L_sire + L_dam))))
      Alpha[each]=1/(L[each,each]^2)
      j = i[i!=each]
      for (individual in j){
        L[each,individual] = 0.5*(L[i_sire,individual]+L[i_dam,individual])
      }
      Ainv[each,each] = Ainv[each,each] + Alpha[each]
      Ainv[i_sire,each] = Ainv[i_sire,each] - (Alpha[each]/2)
      Ainv[each,i_sire] = Ainv[each,i_sire] - (Alpha[each]/2)
      Ainv[i_dam,each] = Ainv[i_dam,each] - (Alpha[each]/2)
      Ainv[each,i_dam] = Ainv[each,i_dam] - (Alpha[each]/2)
      Ainv[i_sire,i_sire] = Ainv[i_sire,i_sire] + (Alpha[each]/4)
      Ainv[i_sire,i_dam] = Ainv[i_sire,i_dam] + (Alpha[each]/4)
      Ainv[i_dam,i_sire] = Ainv[i_dam,i_sire] + (Alpha[each]/4)
      Ainv[i_dam,i_dam] = Ainv[i_dam,i_dam] + (Alpha[each]/4)
    }
    
  }
  return(Ainv)
}

#TEST

ID = c(1:6)
sire = c(rep(0,2),1,1,4,5)
dam = c(rep(0,2),2,0,3,2)
pedigree = data.frame(cbind(ID,sire,dam))

Ainv_QuaasA = QuaasA(pedigree)


