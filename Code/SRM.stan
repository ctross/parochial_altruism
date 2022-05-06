
data{
 int N;
 int P;
 int M;

 real Q[2];

 int Giving[N,N];
 int Taking[N,N];
 int Reducing[N,N];
 int Transfer[N,N];
 int Friends[N,N];

 matrix[N,N] SameEthnicity;
 matrix[N,N] SameSex;
 matrix[N,N] Other;
 matrix[N,N] Relatedness;
 matrix[N,N] Marriage;

 vector[N] Age;
 vector[N] Male;
 vector[N] CantWork;
 vector[N] Grip;
 vector[N] Sad;
 vector[N] NoFood;
 vector[N] GoodsValues;
 vector[N] Indigenous;
 vector[N] NotThere;
 vector[N] MissingFocal;
}

parameters {
 vector[P] BG; 
 vector[P] BL; 
 vector[P] BP; 
 vector[P] BT; 
 vector[P] BF; 

 vector[2] FA_raw[N, M]; 
 vector<lower=0>[2] FA_SD[M];
 cholesky_factor_corr[2] FA_L_chol[M];

 matrix[N,N] D_raw[M];
 real<lower=0> D_SD[M];
 cholesky_factor_corr[2] D_L_chol[M];
}

model{
//######################################################## Local storage
 real FocalFactors;
 vector[N] TargetFactors;
 row_vector[N] DyadFactors;
 vector[N] Theta;

 vector[N] F[M];
 vector[N] A[M];

 vector[2] FA[N, M];  
 matrix[N,N] D[M];
 vector[2] scrap;

//######################################################## Priors
 BG ~ normal(0,1.5);
 
 BL ~ normal(0,1.5);

 BP ~ normal(0,1.5);
 
 BT ~ normal(0,1.5);

 BF ~ normal(0,1.5);

 for(i in 1:N){
 for(m in 1:M){
  FA_raw[i,m] ~ normal(0,1);
 }}

 for(m in 1:M){
  to_vector(D_raw[m]) ~ normal(0,1);
 }

for(m in 1:M){
 FA_SD[m] ~ exponential(1.5);
 FA_L_chol[m] ~ lkj_corr_cholesky(2);
}

 for(m in 1:M){
 D_SD[m] ~ exponential(1.5);
 D_L_chol[m] ~ lkj_corr_cholesky(2);
}

//######################################################## Transformed Priors


 for(i in 1:N){
  for(m in 1:M){
   FA[i,m] = FA_SD[m] .* (FA_L_chol[m]*FA_raw[i, m]);
  }}

 for(m in 1:M){
 for(i in 1:N){
  F[m,i] = FA[i,m,1];
  A[m,i] = FA[i,m,2];
  }}

for(m in 1:M){
   for(i in 1:(N-1)){
   for(j in (i+1):N){  
      scrap[1] = D_raw[m,i,j];
      scrap[2] = D_raw[m,j,i];
      scrap = rep_vector(D_SD[m],2) .* (D_L_chol[m]*scrap);
    D[m,i,j] = scrap[1];           
    D[m,j,i] = scrap[2];                       
    }}
    
   for(i in 1:N){
    D[m,i,i] = -99;                                
    }
  }


//######################################################## Model Allocation Data
 for(i in 1:N){
  if(MissingFocal[i]==0){
  FocalFactors  = Q[1]*(BG[2]*Age[i] + BG[3]*Male[i] + BG[5]*CantWork[i] + BG[6]*Grip[i] + BG[7]*Sad[i] + BG[8]*NoFood[i] + BG[9]*GoodsValues[i]) + 
                  Q[2]*(BG[4]*Indigenous[i]) + 
                  BG[1] + F[1,i];  
  
  TargetFactors = Q[1]*(BG[11]*Age + BG[12]*Male + BG[13]*CantWork + BG[14]*Grip + BG[15]*Sad + BG[16]*NoFood + BG[17]*GoodsValues + BG[18]*NotThere) +
                  A[1];  
    
  DyadFactors   = Q[1]*(BG[19]*Relatedness[i] + BG[20]*Marriage[i] + BG[21]*SameSex[i]) +
                  Q[2]*(BG[22]*SameEthnicity[i]*(1-Indigenous[i]) + BG[23]*SameEthnicity[i]*Indigenous[i]) + 
                  D[1][i];

  Giving[i] ~ multinomial(softmax((FocalFactors + TargetFactors + to_vector(DyadFactors)) .* to_vector(Other[i])));
  }
 }  
 
//######################################################## Model Taking Data
 for(i in 1:N){
     if(MissingFocal[i]==0){
  FocalFactors  = Q[1]*(BL[2]*Age[i] + BL[3]*Male[i] + BL[5]*CantWork[i] + BL[6]*Grip[i] + BL[7]*Sad[i] + BL[8]*NoFood[i] + BL[9]*GoodsValues[i]) + 
                  Q[2]*(BL[4]*Indigenous[i]) + 
                  BL[1] + F[2,i];  
  
  TargetFactors = Q[1]*(BL[11]*Age + BL[12]*Male + BL[13]*CantWork + BL[14]*Grip + BL[15]*Sad + BL[16]*NoFood + BL[17]*GoodsValues + BL[18]*NotThere) +
                  A[2];  
    
  DyadFactors   = Q[1]*(BL[19]*Relatedness[i] + BL[20]*Marriage[i] + BL[21]*SameSex[i]) +
                  Q[2]*(BL[22]*SameEthnicity[i]*(1-Indigenous[i]) + BL[23]*SameEthnicity[i]*Indigenous[i]) + 
                  D[2][i];

  Theta = softmax((FocalFactors + TargetFactors + to_vector(DyadFactors)) .* to_vector(Other[i]));
  Taking[i] ~ multinomial(Theta);  
  }
 }   
 
//######################################################## Model Punishment Data
 for(i in 1:N){
     if(MissingFocal[i]==0){
  FocalFactors  = Q[1]*(BP[2]*Age[i] + BP[3]*Male[i] + BP[5]*CantWork[i] + BP[6]*Grip[i] + BP[7]*Sad[i] + BP[8]*NoFood[i] + BP[9]*GoodsValues[i]) + 
                  Q[2]*(BP[4]*Indigenous[i]) + 
                  BP[1] + F[3,i];  
  
  TargetFactors = Q[1]*(BP[11]*Age + BP[12]*Male + BP[13]*CantWork + BP[14]*Grip + BP[15]*Sad + BP[16]*NoFood + BP[17]*GoodsValues + BP[18]*NotThere) +
                  A[3];  
    
  DyadFactors   = Q[1]*(BP[19]*Relatedness[i] + BP[20]*Marriage[i] + BP[21]*SameSex[i]) +
                  Q[2]*(BP[22]*SameEthnicity[i]*(1-Indigenous[i]) + BP[23]*SameEthnicity[i]*Indigenous[i]) + 
                  D[3][i];

  Reducing[i] ~ multinomial(softmax((FocalFactors + TargetFactors + to_vector(DyadFactors)) .* to_vector(Other[i])));
 }     
 }
 
//######################################################## Model Transfer Data
 for(i in 1:N){
      if(MissingFocal[i]==0){
  FocalFactors  = Q[1]*(BT[2]*Age[i] + BT[3]*Male[i] + BT[5]*CantWork[i] + BT[6]*Grip[i] + BT[7]*Sad[i] + BT[8]*NoFood[i] + BT[9]*GoodsValues[i]) + 
                  Q[2]*(BT[4]*Indigenous[i]) + 
                  BT[1] + F[4,i];  
  
  TargetFactors = Q[1]*(BT[11]*Age + BT[12]*Male + BT[13]*CantWork + BT[14]*Grip + BT[15]*Sad + BT[16]*NoFood + BT[17]*GoodsValues + BT[18]*NotThere) + 
                  A[4];  
    
  DyadFactors   = Q[1]*(BT[19]*Relatedness[i] + BT[20]*Marriage[i] + BT[21]*SameSex[i]) +
                  Q[2]*(BT[22]*SameEthnicity[i]*(1-Indigenous[i]) + BT[23]*SameEthnicity[i]*Indigenous[i]) + 
                  D[4][i];
  
  Theta = softmax((FocalFactors + TargetFactors + to_vector(DyadFactors)) .* to_vector(Other[i]));
  Transfer[i] ~ multinomial(Theta);
   }
 }   

//######################################################## Model Friendship Data
 for(i in 1:N){
      if(MissingFocal[i]==0){
  FocalFactors  = Q[1]*(BF[2]*Age[i] + BF[3]*Male[i] + BF[5]*CantWork[i] + BF[6]*Grip[i] + BF[7]*Sad[i] + BF[8]*NoFood[i] + BF[9]*GoodsValues[i]) + 
                  Q[2]*(BF[4]*Indigenous[i]) + 
                  BF[1] + F[5,i];  
  
  TargetFactors = Q[1]*(BF[11]*Age + BF[12]*Male + BF[13]*CantWork + BF[14]*Grip + BF[15]*Sad + BF[16]*NoFood + BF[17]*GoodsValues + BF[18]*NotThere) +
                  A[5];  
    
  DyadFactors   = Q[1]*(BF[19]*Relatedness[i] + BF[20]*Marriage[i] + BF[21]*SameSex[i]) +
                  Q[2]*(BF[22]*SameEthnicity[i]*(1-Indigenous[i]) + BF[23]*SameEthnicity[i]*Indigenous[i]) + 
                  D[5][i];
  
  Theta = softmax((FocalFactors + TargetFactors + to_vector(DyadFactors)) .* to_vector(Other[i]));
  Friends[i] ~ multinomial(Theta);
   }
 } 

}
                                                     
generated quantities{
 matrix[2,2] FA_corr[M];  
 matrix[2,2] D_corr[M];
 real FA_Rho[M];
 real D_Rho[M];

 for(m in 1:M){
  FA_corr[m] = tcrossprod(FA_L_chol[m]); 
  D_corr[m] = tcrossprod(D_L_chol[m]); 
  FA_Rho[m] = FA_corr[m,1,2]; 
  D_Rho[m] = D_corr[m,1,2];
 }
}   
                 
               
 