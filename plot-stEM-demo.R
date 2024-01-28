# 03/17/2021

# draw stochastic EM demo plots

par(mar = c(2,2,2,2))

y = rnorm(150)
y[1:30] = y[1:30] + seq(from=5,to=0,length.out=30)


# taking the last iteration result
plot(y,type = 'l', yaxt='n')
abline(h = 0, col='gray')
points(x=150,y=y[150],cex=2,
       col='red',lwd=2)

# averaging the last several iters
plot(y,type = 'l', yaxt='n')
abline(h = 0, col='gray')
polygon(x=c(120,120,150,150), y=c(-2,2.5,2.5,-2), 
        border='red',lwd=2)

# multiple chains
y1 = rnorm(150,mean = 2)
y1[1:50] = y1[1:50] + seq(from=-4,to=0,length.out=50)

plot(y,type = 'l', yaxt='n', ylim = c(-5,7))
lines(y1, col='blue')
abline(h = 1, col='gray')
points(x=c(150,150),y=c(y[150],y1[150]),cex=2,
       col='red',lwd=2)

