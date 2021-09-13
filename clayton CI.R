library(copula)
N=500
points.ci = c(0,0)
for(i in 1:100){
rho_ken = 0+i/125
rho_val = 2*rho_ken/(1-rho_ken)
uncorrected_rho=0
for(j in 1:1000){
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

cp <- claytonCopula(param = c(rho_val), dim = 2) #normalCopula(param = rho_val, dim = 2)
multivariate_dist <- mvdc(copula = cp,
              margins = c("norm", "norm"),
              paramMargins = list(list(mean = 0, sd=1),
                                  list(mean = 0, sd=1)) )
x000 <- rMvdc(length(dt), multivariate_dist)


x000 <- x000*sqrt(dt)

Xt=cbind(cumsum(x000[,1]),cumsum(x000[,2]), locate)
St1<-Xt[which(Xt[,4]==0),c(1,3)]
St2<-Xt[which(Xt[,4]==1),c(2,3)]
n1 <- length(St1[,2])
n2 <- length(St2[,2])

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
asynch.Kendall = cor(mathu2[,1], mathu2[,2], method = "kendall")

#uncorrected_est = c(uncorrected_est, rho_cap2)
uncorrected_rho = c(uncorrected_rho, asynch.Kendall)
}
uncorrected_rho = uncorrected_rho[-1]
#rho_est.val = 2*uncorrected_rho/(1-uncorrected_rho)
#quantile(rho_est.val, c(.025, .975))
quan = quantile(uncorrected_rho, c(.025, .975))
points.ci = rbind(points.ci, c(rho_ken, quan[1]), c(rho_ken, quan[2]))
}
points.ci = points.ci[-1,]
## Figure with confidence interval
op <- par(cex = 1.2)
plot(points.ci[,2], points.ci[,1], xlab = "Estimated tau (uncorrected)", ylab="true value of tau")
for(i in 1:100){
points(points.ci[(2*i-1):(2*i),2], points.ci[(2*i-1):(2*i),1], pch=20)
lines(points.ci[(2*i-1):(2*i),2], points.ci[(2*i-1):(2*i),1])
}

# save the outcome for later use: 
save_data = data.frame(points.ci[,2],points.ci[,1])
write.csv(save_data, "clayton_CI_fig_9_left.csv")


## same picture changing the axis
#plot( points.ci[,1],points.ci[,2], ylab = "uncorrected tau", xlab="true tau")
#for(i in 1:30){
#points( points.ci[(2*i-1):(2*i),1], points.ci[(2*i-1):(2*i),2],pch=20)
#lines( points.ci[(2*i-1):(2*i),1],points.ci[(2*i-1):(2*i),2], col="Red")
#}

## Band picture
#plot(points.ci[,1], points.ci[,2])
#points(points.ci[,1], points.ci[,3])
#points.ci.0 = rep(30,0)
#points.ci.1 = rep(30,0)
#points.ci.2 = rep(30,0)
#for(i in 1:30){
#points.ci.0[i]=points.ci[2*i,1]
#points.ci.1[i]=points.ci[2*i-1,2]
#points.ci.2[i]=points.ci[2*i,2]
#} 
#plot(points.ci.2, points.ci.0, xlim=c(-.05,.5), xlab="uncorrected tau",
#		ylab="true tau", pch=20)
#lines(points.ci.2, points.ci.0, col="blue")
#points(points.ci.1, points.ci.0, pch=20)
#lines(points.ci.1, points.ci.0, col="blue")

