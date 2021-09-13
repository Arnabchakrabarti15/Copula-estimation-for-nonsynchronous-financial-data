library(copula)
N=500
compare = c(0,0)
for(i in 1:100){
rho_val = runif(1, 0,1)
uncorrected_est = 0
for(j in 1:10){
T1 = 330
#set.seed(230)
l = sort(runif(N, min = 0, max = T1))
#set.seed(334)
ll = sort(runif(N, min = 0, max = T1))
marker <- c(rep(0,length(l)),rep(1,length(ll)))
loci <- cbind(c(l,ll), marker)
locate <- loci[order(loci[,1]),]
location <- locate[,1]
n=length(location)
dt=diff(c(0,location))
#mean(dt)
#*************

normal <- normalCopula(param = rho_val, dim = 2)
multivariate_dist <- mvdc(copula = normal,
              margins = c("norm", "norm"),
              paramMargins = list(list(mean = 0, sd=1),
                                  list(mean = 0, sd=1)) )
x000 <- rMvdc(length(dt), multivariate_dist)
normal.cop <- normalCopula(dim=2)
fit.cop<- fitCopula(normal.cop,pobs(x000),method="ml")
rho <- coef(fit.cop)
#print(rho)


x000 <- x000*sqrt(dt)
normal.cop <- normalCopula(dim=2)
fit.cop<- fitCopula(normal.cop,pobs(x000),method="ml")
rho <- coef(fit.cop)
#print(rho)


Xt=cbind(cumsum(x000[,1]),cumsum(x000[,2]), locate)
#head(Xt)
#plot(Xt[,3],Xt[,2],type='l',col="Red", xlim=c(0,12), ylim=c(-4,4))
#lines(Xt[,3],Xt[,1],type='l',col="Green")
St1<-Xt[which(Xt[,4]==0),c(1,3)]
St2<-Xt[which(Xt[,4]==1),c(2,3)]
n1 <- length(St1[,2])
n2 <- length(St2[,2])
#plot(St1[,2], St1[,1], ylim=range(c(St1[,1],St2[,1])), xlab="Time", ylab="return")
#lines(St2[,2], St2[,1], col="red")

#*********************
# SYNCHRONIZATION #
n = n1
t = cbind(Xt[which(Xt[,4]==0),3], Xt[which(Xt[,4]==1),3])
y = cbind(Xt[which(Xt[,4]==0),1], Xt[which(Xt[,4]==1),2])
t2 = sort(c(t[,1],t[,2]))
T = max(t[1,1],t[1,2])
if(T==t[1,2]){
	d = t[,1]-T
	k = which.min(abs(d[d<=0]))
	new = c(y[k,1], y[1,2],t[k,1], t[1,2])
	loki = T
}else{
	d=t[,2]-T
	k = which.min(abs(d[d<=0]))
     	new = c(y[1,1], y[k,2], t[1,1], t[k,2])
	loki = T}

while(T <min(t[n,1], t[n,2]))
     {
	df = t[,1]-T
	df2= t[,2]-T
      kf = which.min(abs(df[df<=0]))
	kf2 = which.min(abs(df2[df2<=0]))
	T = max(t[kf+1,1],t[kf2+1,2])
	if(T==t[kf2+1,2]){
		d = t[,1]-T
		k = which.min(abs(d[d<=0]))
		new = rbind(new, c(y[k,1], y[kf2+1,2],t[k,1], t[kf2+1,2]))
            loki = c(loki, T)
	}else{
		d=t[,2]-T
		k = which.min(abs(d[d<=0]))
      	new = rbind(new, c(y[kf+1,1], y[k,2], t[kf+1,1], t[k,2]))
		loki = c(loki, T)}
      }
mathu2 = cbind(diff(new[,1]), diff(new[,2]))
normal.cop <- normalCopula(dim=2)
fit.cop2<- fitCopula(normal.cop,pobs(mathu2),method="ml")
rho_cap2 <- coef(fit.cop2)
#print(rho_cap2)
uncorrected_est = c(uncorrected_est, rho_cap2)
}
unc.est = mean(uncorrected_est[-1])
compare = rbind(compare, c(rho_val, unc.est))
}
xxx = compare[-1,2]
yyy = compare[-1,1]
op = par(cex=1.2)
plot(xxx, yyy, pch=20, xlab= "Estimated (uncorrected) correlation", 
	ylab = "True correlation")
abline(lm(yyy~xxx), col= "red", lwd=2)
summary(lm(yyy~xxx))

save_data = data.frame(xxx,yyy)
write.csv(save_data, "gaussian_fig_8_right.csv")

