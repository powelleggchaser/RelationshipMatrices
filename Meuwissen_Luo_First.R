# Meuwissen and Luo Algorithm 1st Algorithm
# Calculating Inbreeding Based on the L matrix
# Independent Row-wise calculation of L - suitable for updating with new data

#Neither Work

L = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
D = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
A = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
F = rep(0,nrow(L))
#F[1] = -1
i = c(1:nrow(L))
for (each in i){
	ancestors = each
	i_sire = pedigree$sire[each]
  i_dam = pedigree$dam[each]
  sire_range = c(1:i_sire)
  dam_range = c(1:i_dam)
  L_sire = rep(0,i_sire)
  L_dam = rep(0,i_dam)
	L[each,each] = 1
	if ((i_sire!=0)&(i_dam!=0)){
	  ancestors = c(ancestors,i_sire,i_dam)
	  for (person in ancestors){
	    L[each,i_sire] = L[each,i_sire] + (0.5*(L[each,each]))
	    L[each,i_dam] = L[each,i_dam] + (0.5*(L[each,each]))
	    for (ind in sire_range){
	      L_sire[ind] = (L[i_sire,ind]^2) 
	    }
	    L_sire = sum(L_sire)
	    for (indiv in dam_range){
	      L_dam[indiv] = (L[i_dam,indiv]^2) 
	    }
	    L_dam = sum(L_dam)
	    D[each,each] = (0.5 - (0.25*(L_sire + L_dam)))
	    A[each,each] = A[each,each] + ((L[each,j]^2)*(D[j,j]))
	    ancestors = ancestors[ancestors!=j] 
	  }
	}
	if ((i_sire!=0)&(i_dam==0)){
	  for (person in ancestors){
	    ancestors = c(ancestors,i_sire)
	    L[each,i_sire] = L[each,i_sire] + (0.5*(L[each,each]))
	    for (ind in sire_range){
	      L_sire[ind] = (L[i_sire,ind]^2) 
	    }
	    L_sire = sum(L_sire)
	    A[each,each] = A[each,each] + ((L[each,j]^2)*(D[j,j]))
	    D[each,each] = (0.75 - (0.25*(L_sire)))
	    ancestors = ancestors[ancestors!=j]
	  }
	}
	if ((i_sire==0)&(i_dam!=0)){
	  for(person in ancestors){
	    ancestors = c(ancestors,i_dam)
	    L[each,i_dam] = L[each,i_dam] + (0.5*(L[each,each]))
	    for (indiv in dam_range){
	      L_dam[indiv] = (L[i_dam,indiv]^2) 
	    }
	    L_dam = sum(L_dam)
	    D[each,each] = (0.75 - (0.25*(L_dam)))
	    A[each,each] = A[each,each] + ((L[each,j]^2)*(D[j,j]))
	    ancestors = ancestors[ancestors!=j]
	  }
	}
	if ((i_sire==0)&(i_dam=0)){
	  for (person in ancestors){
	    D[each,each] = 1
	    A[each,each] = A[each,each] + ((L[each,j]^2)*(D[j,j]))
	    ancestors = ancestors[ancestors!=j] 
	  }
	}
	F[each] = A[each,each] - 1
}

# Quaas Updated Meuwissen and Luo Algorithm 1st Algorithm


L = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
D = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
A = matrix(0,nrow=length(pedigree$ID),ncol=length(pedigree$ID))
F = rep(-1,nrow(L))
i = c(1:nrow(L))

for (each in i){
  F[each]=0
  male_ancestors = vector()
  female_ancestors = vector()
  i_sire = pedigree$sire[each]
  if (i_sire!=0){ male_ancestors = c(male_ancestors,i_sire); L[i_sire,i_sire]=1}
  i_dam = pedigree$dam[each]
  if (i_dam!=0){ female_ancestors=c(female_ancestors,i_dam); L[i_dam,i_dam]=1}
  if ((length(male_ancestors)>0)&(length(female_ancestors)>0)){
    j = max(male_ancestors); k = max(female_ancestors)
    j_sire = pedigree$sire[j]; j_dam = pedigree$dam[j]
    if (j > k){
      if (j_sire!=0){ male_ancestors= c(male_ancestors,j_sire); L[i_sire,j_sire]=L[i_sire,j_sire]+(0.5*(L[i_sire,j]))}
      if (j_dam!=0){ male_ancestors=c(male_ancestors,j_dam); L[i_sire,j_dam]=L[i_sire,j_dam]+(0.5*(L[i_sire,j]))}
      male_ancestors = male_ancestors[male_ancestors!=j]
    }
    else if (k > j){
      k_sire = pedigree$sire[k]; k_dam = pedigree$dam[k]
      if (k_sire!=0){ female_ancestors=c(female_ancestors,k_sire); L[i_dam,k_sire]=L[i_dam,k_sire]+(0.5*(L[i_dam,k]))}
      if (k_dam!=0){ male_ancestors=c(female_ancestors,k_dam); L[i_dam,i=k_dam]=L[i_dam,k_dam]+(0.5*(L[i_dam,k]))}
      female_ancestors = female_ancestors[female_ancestors!=k]
    }
    else if (j == k){
      if (j_sire!=0){ 
        male_ancestors = c(male_ancestors,j_sire); L[i_sire,j_sire]=L[i_sire,j_sire]+(0.5*(L[i_sire,j]))
        female_ancestors=c(female_ancestors,j_sire); L[i_dam,j_sire]=L[i_dam,j_sire]+(0.5*(L[i_dam,j]))
      }
      if (j_dam!=0){ 
        male_ancestors = c(male_ancestors,j_dam); L[i_sire,j_dam]=L[i_sire,j_dam]+(0.5*(L[i_sire,j]))
        female_ancestors=c(female_ancestors,j_dam); L[i_dam,j_dam]=L[i_dam,j_dam]+(0.5*(L[i_dam,j]))
      }
    }
    F[each] = F[each] + ((L[i_sire,j])*(L[i_dam,j])*(0.5)*(D[j,j]))
  }
}

