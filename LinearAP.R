
# ReplicateSim <- function(init, theta., K){
#   init_g <- init
#   for (i in 1:K){
#     sample_adj <- SimLinearAP(init = init_g, theta = theta., n = i*1000)
#     sample_g <- graph_from_adjacency_matrix(sample_adj$adj, mode = "directed", weighted = NULL, diag = TRUE)
#     init_g <- sample_g
#   }
#   return(sample_adj)
# }
# 
# g <- ReplicateSim(init = init_g, theta. = params, K = 4)
# sample_g <- graph_from_adjacency_matrix(g$adj, mode = "directed", weighted = NULL, diag = TRUE)


NodeSample <- function(edge, j, nedge, nnode, delta){
  W <- runif(1, 0, nedge+nnode*delta)
  if (W <= nedge){ v <- edge[ceiling(W), j]}
  else {v <- ceiling((W - nedge)/delta)}
  stopifnot(v <= nnode)
  return(v)
}

SimLinearAP <- function(init, theta, n) {
  init_adj <- get.adjacency(init)
  
  adj <- matrix(init_adj, nrow = dim(init_adj)[1], ncol = dim(init_adj)[1])
  edge_list <- get.edgelist(init)
  
  # main code 
  t <- sum(adj) #number of edges
  stopifnot(t == dim(edge_list)[1])
  N <- dim(adj)[1] #number of nodes
  while (t < n){
    scheme <- runif(1, 0, 1)
    if (scheme < theta[1]){  
      v_1 <- N + 1        #Assign the newest node index to v_1
      v_2 <- NodeSample(edge_list, 2, t, N, theta[3])
      N <- N + 1        #Update Node count
      #update adjacency matrix
      adj <- cbind(adj, 0); adj <- rbind(adj, 0)
      }
    else if ((theta[1] <= scheme) & (scheme < theta[1] + theta[2] )){
      v_1 <- NodeSample(edge_list, 1, t, N, theta[4])
      v_2 <- NodeSample(edge_list, 2, t, N, theta[3])
    }
    else if (scheme >= theta[1]+theta[2]){
      v_1 <- NodeSample(edge_list, 1, t, N, theta[4])
      v_2 <- N + 1 
      N <- N + 1
      #update adjacency matrix
      adj <- cbind(adj, 0); adj <- rbind(adj, 0)
    }
    t <- t + 1
    edge_list <- rbind(edge_list, c(v_1 , v_2))
    adj[v_1, v_2] <- adj[v_1, v_2] + 1
    
    #Make sure node and edge count add up
    stopifnot(t == dim(edge_list)[1])
    stopifnot(N == dim(adj)[1])
  }
  return(list(adj = adj, edge_list = edge_list))
}

## Estimation & Inference 
library(rootSolve); library(pracma)
EstLinearAP <- function(adj, r, normalize) {
  func <- function(x, degree) c(1:max(degree))/(c(1:max(degree)) + rep(x, max(degree)))
  t <- sum(adj)
  
  #Beta 
  B <- 1 - (dim(adj)[1]/t)
  #Delta_in
  in_degree <- colSums(adj)
  N.0 <- sum(in_degree == 0)
  N.i <- rep(0, max(in_degree))
  for (j in 1:max(in_degree)){
    N.i[j] <- sum(in_degree > j)
  }
  fun_delta_in <- function(delta_in) {
    sum(N.i*func(x = delta_in, degree = in_degree))*(1 + delta_in*(1 - B))/t - (((N.0/t)+ B)/(1 - ((N.0/t))*(delta_in/(
     1+(1-B)*delta_in))))
  }
  sol_delta_in <- uniroot(fun_delta_in, c(0.0001, r))$root

  #Alpha
  A <- ((N.0/t)+B)/(1- (N.0/t)*(sol_delta_in/(1+(1-B)*sol_delta_in)))-B
  
  #Delta_out
  out_degree <- rowSums(adj)
  N.0_out <- sum(out_degree==0)
  N.i_out <- rep(0, max(out_degree))
  for (j in 1:max(out_degree)){
    N.i_out[j] <- sum(out_degree > j)
  }
  fun_delta_out <- function(delta_out) {
    sum(N.i_out*func(x = delta_out, degree = out_degree))*(1 + delta_out*(1 - B))/t - (((N.0_out/t)+ B)/(1 - ((N.0_out/t))*(delta_out/(
      1+(1-B)*delta_out))))
  }
  sol_delta_out <- uniroot(fun_delta_out, c(0.0001, r), extendInt = "yes")$root
  
  #Gamma
  G <- ((N.0_out/t)+B)/(1- (N.0_out/t)*(sol_delta_out/(1+(1-B)*sol_delta_out)))-B
  
  if (normalize == 1){
    A_0 <- A
    A <- (A_0*(1-B))/(A_0+G)
    G <- (G*(1-B))/(A_0 + G)
    
    fun_in <- function(d_in){
      (1/t)*(sum(in_degree > 0)/d_in + sum(N.i/(c(1:max(in_degree))+rep(d_in, max(in_degree))))) - (1-A-B)/d_in -
        (A+B)*(1-B)/(1+(1-B)*d_in)
    }
    sol_delta_in <- uniroot(fun_in, c(0.0001, r))$root
    
    fun_out <- function(d_out){
      (1/t)*(sum(out_degree > 0)/d_out + sum(N.i_out/(c(1:max(out_degree))+rep(d_out, max(out_degree))))) - (1-G-B)/d_out -
        (G+B)*(1-B)/(1+(1-B)*d_out)
    }
    sol_delta_out <- uniroot(fun_out, c(0.0001, r))$root
  }
  
  return(list(Beta = B, delta_in = sol_delta_in, Alpha = A, delta_out = sol_delta_out, Gamma = G))
}

est <- EstLinearAP(adj = sample_adj$adj, r = 3, normalize = 1)



library(igraph)
set.seed(134)

init_g <- sample_pa(20, directed = TRUE)
params <- c(0.4, 0.2, 1, 1)

sample_adj <- SimLinearAP(init = init_g, theta = params, n = 1000)
sample_g <- graph_from_adjacency_matrix(sample_adj$adj, mode = "directed", weighted = NULL, diag = TRUE)

set.seed(120)
l <- layout.fruchterman.reingold(sample_g)
par(mar = c(1,1,1,1))

plot(sample_g, layout=l, vertex.size = 0.4*degree(sample_g, mode = "total"), 
     vertex.label = NA, vertex.shape = c("circle"), 
     edge.width = 0.5, edge.arrow.size=0.05, edge.arrow.width=0.5,
     vertex.color=c("blue"), vertex.border = c("blue"))


# in_degree <- colSums(sample_adj$adj)
# out_degree <- rowSums(sample_adj$adj)
# 
save(sample_g,l, file="/Users/MacUser/Desktop/Network/Final Project/SampleSim.RData")



#### french
source("/Users/MacUser/Desktop/Network/Final Project/Functions.R")

french <- read.delim("/Users/MacUser/Desktop/Network/Final Project/french.txt", sep = " ", header = FALSE)
french <- french[, -3]
french <- graph.data.frame(french, directed = TRUE)
french_adj <- matrix(as_adjacency_matrix(french), nrow = 8325, ncol = 8325)
french_theta <- EstLinearAP(adj = french_adj, normalize = 0)

read.table(file = "/Users/MacUser/Desktop/Network/Final Project/out.patentcite",
           sep = '\t', header = FALSE)

dd.sample_g <- degree.distribution(french, mode= "in")
d <- (0:(max(degree(french, mode = "in"))-1))
ind <- (dd.sample_g!=0)
plot(d[ind], dd.sample_g[ind], log = "xy", col = "blue", 
     xlab = "Log-Degree", ylab = "Log-Frequency",  
     main = "Log-Log Degree Distribution")


---
  title: "Fitting a Random Linear Preferential Attachment Model for Directed Graphs"
number_sections: True
output: pdf_document
subtitle: 'Jimmy Ting-Yuan Kuo'
header-includes:
  - \setlength\parindent{24pt}
bibliography: MA703.bib
---

