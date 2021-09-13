library(copula)
N=5000
est = 0
uncorrected_est = 0
R=0
for(j in 1:100){
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
rho_val = 0.8
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


new22 = 0
over = 0
di1 = 0
di2 = 0 
for(i in 2:length(new[,1])){
	if(new[i,4]<new[i,3] && new[i-1,4]<new[i-1,3]){
		overlap = new[i,4] - new[i-1,3]
		dino1 = new[i,4] - new[i-1,4]
		dino2 = new[i,3] - new[i-1,3]
		}else if(new[i,4]<new[i,3] && new[i-1,4]>=new[i-1,3]){
			overlap = new[i,4] - new[i-1,4]
			dino1 = new[i,4] - new[i-1,4]
			dino2 = new[i,3] - new[i-1,3]
			}else if(new[i,4]>=new[i,3] && new[i-1,4]<new[i-1,3]){
				overlap = new[i,3] - new[i-1,3]
				dino1 = new[i,4] - new[i-1,4]
				dino2 = new[i,3] - new[i-1,3]
				}else{
					overlap = new[i,3] - new[i-1,4]
					dino1 = new[i,4] - new[i-1,4]
					dino2 = new[i,3] - new[i-1,3]
					}
	#print(overlap)
	#print(dino)
	#const = sqrt(new21[i,4]-new21[i-1,4])*sqrt(new21[i,3]-new21[i-1,3])/overlap	
	new22 = rbind(new22, overlap/sqrt(dino1*dino2))	#new23 = rbind(new23, c(new21[i,1], new21[i,2]))
	over = c(over,overlap)
	di1 = c(di1,dino1)
	di2 = c(di2,dino2)
}
new22 = new22[-1,]
over = over[-1]
di1 = di1[-1]
di2 = di2[-1]
coco = mean(over)/(sqrt(mean(di1))*sqrt(mean(di2)))
#-.4*coco
#-.4*mean(new22)
#print(rho_cap2)
unbiased_estimate = rho_cap2/coco
est = c(est, unbiased_estimate)
uncorrected_est = c(uncorrected_est, rho_cap2)

# SYNCHRONIZATION #
n = n1
t = cbind(Xt[which(Xt[,4]==0),3], Xt[which(Xt[,4]==1),3])
y = cbind(Xt[which(Xt[,4]==0),1], Xt[which(Xt[,4]==1),2])
k = N/2
T2 = min(max(t[,1]),max(t[,2]))-max(min(t[,1]),min(t[,2]))
tttt = seq(from = max(min(t[,1]),min(t[,2])), to = min(max(t[,1]),max(t[,2])), by = T2/k)
tttt = tttt[-1]
pair = c(0,0,0)
for(i in 1:length(tttt)){
	d = t[,1]-tttt[i]
	m = which.min(abs(d[d<=0]))
	t[m,1]
	d2 = t[,2]-tttt[i]
	m2 = which.min(abs(d2[d2<=0]))
	if(t[m,1]!=0 && t[m2,2]!=0){
		pair<-rbind(pair, c(y[m,1],y[m2,2], tttt[i]))
		}	
}    
pair = pair[-1,]
#**************************
dif_var1 <- diff(pair[,1])
dif_var2 <- diff(pair[,2])
mat <- cbind(dif_var1, dif_var2)
normal.cop <- normalCopula(dim=2)
fit.cop<- fitCopula(normal.cop,pobs(mat),method="ml")
# Coefficients
rho <- coef(fit.cop)
R = c(R, rho)
}
mean(est[-1])
sd(est[-1])
mean(uncorrected_est[-1])
sd(uncorrected_est[-1])
mean(R[-1])
sd(R[-1])
op = par(cex=1.2)
boxplot(est[-1],uncorrected_est[-1], R[-1], 
	names=c("Corrected", "Refresh time", "Previous tick"), 
	col = c("skyblue", "indianred1", "peachpuff"), border = "brown")
abline(h = rho_val, col="blue")

save_data = data.frame(est[-1],uncorrected_est[-1], R[-1])
write.csv(save_data, "boxplot2.csv")
