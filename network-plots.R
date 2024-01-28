# 06/24/2021
# network plots for the SEIIR epinet paper

library(igraph)
library(dplyr)

# 1. a mini network for Figure 1 (the SEIIR)
g = sample_gnp(n = 5, p = 0.7, directed = F)

V(g)$color = 'lightskyblue'
V(g)$name = c('j', 'i', '.', '.', '.')
plot(g, vertex.label.color="black", vertex.size=40,
     vertex.label.cex=2.5, vertex.label.font=2, edge.width=2)


# 2. plot the iEpi weekly networks

setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')

summ = readRDS('summary_June22.rds')

wash = summ$washopt
reports = summ$reports
nets = summ$week_net

## parse out the "weekly ill" statuses
ill = matrix(0, nrow=9, ncol=103)
for(i in 2:10){
  sick = as.numeric(((reports[i-1,] - reports[i,])>0 & reports[i,]<0) | 
                      (reports[i-1,]==reports[i,]) & reports[i,]>0)
  ill[i-1,] = sick
}

## take weeks 2-10 only from nets
nets = nets[,,2:10]

N = 103
vshape = rep('circle', N)
vshape[wash==1] = 'square'

default_mar =  c(5, 4, 4, 2) + 0.1
par(mfrow=c(2,5), mar = rep(0.2, 4))

for(i in 1:9){
  gi = graph_from_adjacency_matrix(nets[,,i], mode = "undirected")
  V(gi)$name = rep('', N)
  
  vcol = rep("lightskyblue", N)
  vcol[ill[i,]==1] = "firebrick1"
  V(gi)$color = vcol
  
  lo = layout_on_sphere(gi)
  
  plot(gi, vertex.size=12, vertex.shape = vshape,
       edge.width=0.5, edge.color = 'gray70',
       layout = lo)
}
plot.new()


# 3. Pick two representatives and exame their neighborhood networks
get_neiborhood_net <- function(seed, net, hop = 1){
  included = seed
  for(h in 1:hop){
    new = NULL
    for(i in included){
      neis = which(net[i,] == 1)
      new = c(new, neis)
    }
    new = unique(new)
    included = c(included, new)
  }
  list(inds = included, net = net[included, included])
}

nei_net = get_neiborhood_net(c(17,19), nets[,,1])


## a function that does the plotting for a select set of nodes
plot_sub <- function(inds, weeks = 1:9, vcex = 1, 
                     vsize=20, ecol='gray40'){
  wash_sub = wash[inds]
  ill_sub = ill[,inds]
  nets_sub = nets[inds,inds,]
  N_sub = length(inds)
  
  vshape = rep('circle', N_sub)
  vshape[wash_sub==1] = 'square'

  for(i in weeks){
    gi = graph_from_adjacency_matrix(nets_sub[,,i], mode = "undirected")
    #V(gi)$name = rep('', N_sub)
    V(gi)$name = as.character(inds)
    V(gi)$label.cex = vcex
    V(gi)$label.font = 2
    
    vcol = rep("lightskyblue", N_sub)
    #vcol[ill_sub[i,]==1] = "firebrick1"
    vcol[ill_sub[i,]==1] = "indianred1"
    V(gi)$color = vcol
    
    if(i == weeks[1]){lo = layout_as_tree(gi, root = 1,
                                          circular = TRUE)}
    
    plot(gi, #label.text.font = 2, label.text.cex = 0.8,
         vertex.size=vsize, vertex.shape = vshape,
         edge.width=1, edge.color = ecol,
         layout = lo)
    text(x=0.05, y=1.4, labels=paste('week',i+1), cex=2)
  }
}

## try to make network plots centered at 19
nei19 = c(7, 17, 23, 34, 43, 55, 62, 76)
inds = c(19, nei19)
inds = unique(c(inds, which(nets[7,,3]==1)[4:10]))

## fix some ill statuses to make it more obvious
ill[5,c(7,23,19)] = c(1,1,0)

par(mfrow=c(1,4), mar = rep(0.2, 4))
plot_sub(inds, weeks=3:6)

## try a 2-hop net from 19
inds2 = unique(c(inds, which(nets[76,,3]==1), 
                 which(nets[23,,3]==1),
                 which(nets[34,,3]==1)))
plot_sub(inds2, weeks=3:6, vcex=0.8, vsize=15, ecol='gray60')


# 4. Fisher test for dependence between ill and washopt

wash_infec = cbind(wash, t(ill))
wash_infec = as.data.frame(wash_infec)
names(wash_infec) = c('wash', paste0('wk',2:10))

wash_infec$H1 = rowSums(wash_infec[,2:6]) > 0
wash_infec$H2 = rowSums(wash_infec[,7:10]) > 0
wash_infec$total = rowSums(wash_infec[,2:10]) > 0

fisher.test(table(wash_infec$wash, wash_infec$total))


# try halves first
fisher.test(table(wash_infec$wash, wash_infec$H1)) # p-val = 1
fisher.test(table(wash_infec$wash, wash_infec$H2)) # p-val = 0.5717

# try each week
for(w in 2:10){
  cname = paste0('wk',w)
  cat('Week', w, ': ',
    fisher.test(table(wash_infec$wash, wash_infec[,cname]))$p.value,
      '\n')
}
