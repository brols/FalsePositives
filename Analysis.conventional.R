setwd("C:\\Dropbox\\R\\False_Positives_FGSP\\Results\\")
source("C:\\Dropbox\\R\\False_Positives_FGSP\\DataManip.R")
library ("R2jags")
library("R2WinBUGS")
library("BRugs")
library ("mailR")
# data is in a list called "dat"
# Y and W are 3D arrays, nsite x nvisit x nyear 
# date has been centered and standardized
######################################################
## Conventional occupancy model, data remains split into Y and W matrices
## equal p11, omit r11
######################################################
sink( "Dynamic_Multiple-Methods_FGSP_conventional.txt") 
cat("
## SINGLE-METHOD p11 equal ##
model{

## PRIORS
for (m in 1:M){
  psi1[m] ~ dunif(0, 1)} #m
  p11.b1 ~ dunif(-10, 10)
  p11.b3 ~ dunif(-10, 10)
  phi.mu ~ dunif(-10, 10)
  gamma.mu ~ dunif(-10, 10)
  p11.mu ~ dunif(-10, 10)
  phi.tau <-phi.sigma*phi.sigma
  gamma.tau <-gamma.sigma*gamma.sigma
  p11.tau <-p11.sigma*p11.sigma
  phi.sigma ~ dgamma(.1, .1)
  gamma.sigma ~ dgamma(.1, .1)
  p11.sigma ~ dgamma(.1, .1)

## LIKELIHOOD
## Ecological Model
  for (i in 1:nsite){
    z[i,1]~dbern(psi1[agg[i]])
    for(j in 2:nyear){
      muZ[i,j]<-z[i,j-1]*phi[i,j-1] + (1-z[i,j-1])*gamma[i,j-1]
      z[i,j] ~ dbern(muZ[i,j])
        }#j year
      } # i sites 
  for (i in 1:nsite){
    for (j in 2:nyear){
      logit(phi[i,j-1])<- phi.a[agg[i],j-1]
      logit(gamma[i,j-1])<- gamma.a[agg[i],j-1]
    }}
  for (m in 1:M){
    for (j in 2:nyear){
      phi.a[m,j-1] ~ dnorm (phi.mu, phi.tau)
      gamma.a[m,j-1] ~ dnorm (gamma.mu, gamma.tau)
    }}

## Observation model
  for (t in 1:nvisit){
    for (j in 1:nyear){
      for(i in 1:nsite){
        logit(p11[i,t,j])<- p11.a[agg[i],j] + p11.b1*date[i,t,j] + p11.b3*hr[i,t,j] 
          } #i 
        } #j 
      }#t
  for (m in 1:M){
    for (j in 1:nyear){ 
      p11.a[m,j] ~ dnorm (p11.mu, p11.tau)
}}

  for (j in 1:nyear){
    for(i in 1:nsite){
      for (t in 1:nvisit){
      Y[i,t,j] ~ dbern(p[i,t,j])
      p[i,t,j] <- z[i,j]*p11[i,t,j]     
      r[i,t,j] <- z[i,j]*p11[i,t,j]      
        } # t uncertain visits
      for (k in 1:nvisit){
        W[i,k,j] ~ dbern(r[i,k,j])
          } # k certain visits
        }# j year
      } # i sites

# DERIVED PARAMETERS
# Within aggregations
for (m in 1:M){
for (j in 1:nyear){
  p11.2[m,j]<-ilogit(p11.a[m,j] + p11.b1*mean(date[,,]) +  p11.b3*mean(hr[,,]))
  } # j 
  for (j in 1:nyear-1){
    gamma.2[m,j]<-ilogit(gamma.a[m,j])
    phi.2[m,j]<-ilogit(phi.a[m,j])
    }} #j #m
for (m in 1:M){
  psi[m,1] <- psi1[m]
  } #m

n.occ[1]<-sum(z[1:nsite,1])

for (j in 2:nyear){
  n.occ[j] <- sum(z[1:nsite,j])
  for(m in 1:M){
    psi[m,j] <- psi[m,j-1]*phi.2[m,j-1] + (1-psi[m,j-1])*gamma.2[m,j-1]
    turnover[m,j-1] <- (1-psi[m,j-1]) * gamma.2[m,j-1]/psi[m,j]
    growthr[m,j-1] <- psi[m,j]/psi[m,j-1]
}} #j #m

p11.mu2<-ilogit(p11.mu)
gamma.mu2<-ilogit(gamma.mu)
phi.mu2<-ilogit(phi.mu)
p11.sigma2<- ilogit(p11.sigma)
gamma.sigma2<- ilogit(gamma.sigma)
phi.sigma2<- ilogit(phi.sigma)

} # end model
    ",fill=TRUE)
sink()

n.chains=5
n.thin=5
n.iter=5000
n.burnin=floor(n.iter/2)
params<-c("p11.b1",  "p11.b3",
          "phi.mu2", "phi.sigma2", "gamma.mu2", "gamma.sigma2",
          "p11.mu2", "p11.sigma2", 
          "p11.2", "gamma.2", "phi.2"
          ,"psi", "n.occ", "growthr", "turnover")

inits <- function()list (list(psi1 = runif(3, 0.01, 0.99), 
                              z = G                      
)) 
 
st<-system.time(
out.j.conv<-jags(data=dat, 
              inits=c(inits(),inits(),inits(),inits(),inits()),
              model.file="Dynamic_Multiple-Methods_FGSP_conventional.txt",
              n.iter=n.iter, n.chains=n.chains, n.burnin=n.burnin, 
              parameters.to.save=params))
save.image("conv_modeloutput_20_Nov.RData")
send.mail(from="brianrolek@gmail.com", to="3347043143@messaging.sprintpcs.com", subject="Model finished!", 
         body=paste("Your conventional model has finished processing.", date(), st[3]),
         smtp= list(host.name= "smtp.gmail.com", port=465, 
         user.name="brianrolek", passwd="######", ssl=T),
         authenticate=T, send=T)

