data = read.csv("dddd4.csv")# Facebook
timeDatech<- as.character(data[,1])
hour = as.numeric(substring(timeDatech,12,13))
plot(hour)
hour = (hour)*3600
minute = as.numeric(substring(timeDatech,15,16))
minute = minute*60
second = as.numeric(substring(timeDatech,18,23))
time = hour+minute+second
newdata4 = data.frame(time,data[,2])

data2 = read.csv("dddd6.csv")# Apple
timeDatech2<- as.character(data2[,1])
hour2 = as.numeric(substring(timeDatech2,12,13))
plot(hour2)
hour2 = (hour2)*3600
minute2 = as.numeric(substring(timeDatech2,15,16))
minute2 = minute2*60
second2 = as.numeric(substring(timeDatech2,18,23))
time2 = hour2+minute2+second2
newdata6 = data.frame(time2,data2[,2])

#xxx = cbind( c(head(time), head(time2)) ,c(rep(1,6),rep(2,6)))
#plot(xxx[,1], xxx[,2], ylim=c(0,3),yaxt='n', pch=19)

op = par(cex=1.2)
plot(head(time), rep(1,6), ylim=c(.8,2), pch=19, col="red" ,ylab="", yaxt='n', xlab="Time")
points(head(time2), rep(1.5,6), pch= 17, col = "blue", ,yaxt='n')
legend(0.5,2, legend=c("Facebook","Apple"), 
		col=c("red", "blue"), pch=c(19,17))

par(mar = c(5, 5, 5, 5))
op = par(cex=1.2)
plot(time[1:17], rep(1,17), ylim=c(.7,1.8),xlim=c(0,6.5), pch=19, col="red" ,ylab="", yaxt='n', xlab="Time")
points(time2[1:10], rep(1.5,10), pch= 17, col = "blue", ,yaxt='n')
abline(v=time[1:17], col="gray", lty=2)
abline(v=time2[1:10], col="gray", lty=2)
axis(2, c(1,1.5), las=2,  labels=c("Facebook","Apple"))



