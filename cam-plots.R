argv <- commandArgs(trailingOnly=T)

t<-read.csv(argv)
  t<-na.omit(t)
site <- substr(argv,1,(nchar(argv)-4))
years <- max(unique(as.numeric(t$year)))-1
 
x<-rep(1:365, years)
  x<-x[180:length(x)]
    x<-c(x, 1:abs(length(x)-length(t$ag.mdBiomass)))
     
i<-seq(from=1,to=length(t$ag.mdBiomass), by=100)

if(!file.exists("graphs")) dir.create("graphs")
png(height=800, width=1200, filename=file.path('graphs',paste(site,".png",sep="")))
par(mfrow=c(3,1))

plot(t$ag.mdBiomass, xaxt='n', main=paste("site no:",site), xlab="day of year", ylab="median AgB(g)/individual/m2", type="l", col="white")
  grid(lwd=1.2)
    lines(t$ag.mdBiomass, xaxt='n', xlab="day of year", ylab="median biomass(g)/individual/m2", type="l", col="#0066CC", lwd=2)
      axis(1, at=i, labels=x[i])

plot(t$p.size, xaxt='n', main=paste("site no:",site), xlab="day of year", ylab="# of individuals / m2", type="l", col="white")
  grid(lwd=1.2)
    lines(t$p.size, xaxt='n', xlab="day of year", ylab="# of individuals / m2)", type="l", col="#009933", lwd=2)
      axis(1, at=i, labels=x[i])
      
plot(t$sb.size, xaxt='n', main=paste("site no:",site), xlab="day of year", ylab="# of seeds / m2", type="l", col="white")
  grid(lwd=1.2)
    lines(t$sb.size, xaxt='n', xlab="day of year", ylab="median biomass of living plants (g)", type="l", col="#FF9933", lwd=2)
      axis(1, at=i, labels=x[i])
      
graphics.off()
