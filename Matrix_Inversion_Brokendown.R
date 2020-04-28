######## Simple Script to demonstrate the process of Matrix Inversion ########

## Create a 4x4 Matrix
#X = matrix(c(3,1,1,1 ,1,1,0,0, 1,0,1,0, 1,0,0,1),byrow = T,nrow=4)
X = matrix(c(4,2,2,2,1,1,2,1,1),byrow = T,nrow=3)
#function to find inverse
solve(X) #produces singlular matrix

#### ---- LETS BUILD INVERSE MANUALLY ----- #####

ncol = ncol(X)
nrow = nrow(X)

#### ---- Build matrix of Minors ---- ####

# Calculating Minors involves working through the matrix element-wise, and removing values from the
# row and column of each matrix. The remaining values are then solved for their determinant. D = (a11.a22 - a21.a12). 
# The determinant is then stored in the corresponding location of the orinigal matrix element.
# After working through each element of the matrix, a new matrix, The Matrix of Minors will have been produced.
# Easieset to learn with a 3x3 matrix, initially.
m_minors = matrix(0,ncol=ncol,nrow=nrow)
for (r in (1:nrow)){
 for (c in (1:ncol)){
   m = X[-c(r),-c(c)]
   m_minors[r,c] = (m[1,1]%*%m[2,2])-(m[1,2]%*%m[2,1])
 }
}

# Build a matrix of cofactors. A matrix of alternate +s and -s is multiplied with the Matrix of Minors 
#cof = matrix(c(1,-1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1),byrow=T,nrow=4)
cof = matrix(c(1,-1,1,-1,1,-1,1,-1,1),byrow=T,nrow=3)
m_cofactors = cof*m_minors

#Buld the Adjugate matrix by transposing the matrix of cofactors
m_adjugate = t(m_cofactors)

#multiple top row of original matrix (X) by the top row of the adjugte matrix to generate determinant
det_m = X[1,]*m_adjugate[1,]
determinant = sum(det_m)

# Build the Inverse of X matrix by multiplying the adjugate matrix by (1/determinant)
X_inv = (1/determinant)*m_adjugate

