
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
EstLinearAP <- function(adj, r , normalize = 1) {
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
  sol_delta_in <- uniroot(fun_delta_in, c(0.0001, r), extendInt = "yes")$root
  
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
    sol_delta_in <- uniroot(fun_in, c(-r, r), extendInt = "yes")$root
    
    fun_out <- function(d_out){
      (1/t)*(sum(out_degree > 0)/d_out + sum(N.i_out/(c(1:max(out_degree))+rep(d_out, max(out_degree))))) - (1-G-B)/d_out -
        (G+B)*(1-B)/(1+(1-B)*d_out)
    }
    sol_delta_out <- uniroot(fun_out, c(-r, r), extendInt = "yes")$root
  }
  
  return(list(Beta = B, delta_in = sol_delta_in, Alpha = A, delta_out = sol_delta_out, Gamma = G))
}


PlotDegree <- function(graph1, graph2, direction){
  dd1 <- degree.distribution(graph1, mode = direction)
  d1 <- (0:(max(degree(graph1, mode = direction))-1))
  ind1 <- (dd1!=0)
  
  dd2 <- degree.distribution(graph2, mode = direction)
  d2 <- (0:(max(degree(graph2, mode = direction))-1))
  ind2 <- (dd2!=0)
  
  if (direction == "in"){
    plot(d1[ind1], dd1[ind1], log = "xy", col = "blue", ylim=c(min(dd1[ind1], dd2[ind2]), 
       max(dd1[ind1], dd2[ind2]) ), cex = 0.1, cex.lab = 0.1,
         xlab = "Log-Degree", ylab = "Log-Frequency",  
         main = "Log-Log In-Degree", cex.main = 0.7, cex.axis = 0.5)
    points(d2[ind2], dd2[ind2], log = "xy", col = "red", cex = 0.1, cex.lab = 0.1)
   }
  if (direction == "out"){
    plot(d1[ind1], dd1[ind1], log = "xy", col = "blue", ylim=c(min(dd1[ind1], dd2[ind2]), 
        max(dd1[ind1], dd2[ind2]) ), cex = 0.1, cex.lab = 0.1,
         xlab = "Log-Degree", ylab = "Log-Frequency",  
         main = "Log-Log Out-Degree", cex.main = 0.7, cex.axis = 0.5)    
    points(d2[ind2], dd2[ind2], log = "xy", col = "red", cex = 0.1, cex.lab = 0.1)
    
  }
  if (direction == "total"){
    plot(d1[ind1], dd1[ind1], log = "xy", col = "blue", ylim=c(min(dd1[ind1], dd2[ind2]), 
       max(dd1[ind1], dd2[ind2]) ), cex = 0.1, cex.lab = 0.1,
         xlab = "Log-Degree", ylab = "Log-Frequency",  
         main = "Log-Log Total-Degree", cex.main = 0.7, cex.axis = 0.5)    
    points(d2[ind2], dd2[ind2], log = "xy", col = "red", cex = 0.1, cex.lab = 0.1)
  }
  
  legend("topright", legend=c("Simulated", "Real"), col = c("blue", "red"), pch=16, cex=0.5, bty = "n")
}

Partition <- function(matrix, K){
    rowIdx <- seq_len(nrow(matrix))    
    lapply(split(rowIdx, cut(rowIdx, pretty(rowIdx, K))), function(x) matrix[x, ])
}


sort_plot <- function(g,cluster= cluster_infomap, mode="directed"){
    library(igraph)
    library(Matrix)
    if(grepl("atrix",class(g))){
        ig <- graph_from_adjacency_matrix(as.matrix(g)>0,mode=mode)
    }else{
        ig <- g
    }
    n <- length(V(ig))
    if(is.null(cluster)){
        ord <- seq(n)
    } else {
        com <- membership(cluster(ig))
        deg <- rowSums(as.matrix(g[]))
        ord <- order(max(deg)*com+deg)
        deg2 <- colSums(as.matrix(g[]))
        ord2 <- order(max(deg2)*com+ deg2)
    }
    print(image(g[][ord,ord2],
        xlab="",ylab="",sub=""))
    ord
}

CalcExp <- function(theta){
	k_in <- 1 + (1+  theta[3]*(1-theta[2]))/(theta[1] + theta[2])
	k_out <- 1 + (1+  theta[4]*(1-theta[2]))/(1-theta[1])
	return(c(k_in, k_out))
	
}

PowerFit <- function(g, direction){
	m <- displ$new(fun.zero.omit(degree(g, mode = direction))) 
	min <- estimate_xmin(m)
	m$setXmin(min)
	estimate_pars(m)
	bs_p = bootstrap_p(m)
	bs_sd <- sd(bs_p$bootstraps[,3])
	return(c(m$pars, m$xmin, m$pars - 1.96*bs_sd, m$pars + 1.96*bs_sd, bs_p$p))
}

par_dynamics <- function(partition){
  K <- length(partition)
  X <- rep(0, 5)
  for (i in 1:K){
    graph <- graph.data.frame(partition[i], directed=TRUE)
    graph_adj <- matrix(as_adjacency_matrix(graph), nrow = length(V(graph)), 
                        ncol = length(V(graph)))
    est <- EstLinearAP(adj = graph_adj, r = 10)
    bit_theta <- c(est$Alpha, est$Beta, est$delta_in, est$delta_out)
    X <- rbind(X, est)
  }
  X <- X[-1, ]
  return(X)
}

gini_dynamics <- function(partition, direction){
  K <- length(partition)
  X <- rep(0, K)
  for (i in 1:K){
    graph <- graph.data.frame(partition[i], directed = TRUE)
    X[i] <- gini.index(degree(graph, mode = direction))$statistic
  }
  return(X)
}

dynamics <- function(edge, K){
  N <- dim(edge)[1] %/% K
  M <- K + 1
  X <- rep(0, 2)
  for (j in 1:M){
    if (j < M){
      S <- as.integer(j*N); print(S)
      graph <- graph.data.frame(edge[1:S, ], directed = TRUE)
      graph_adj <- matrix(as_adjacency_matrix(graph), nrow = length(V(graph)), 
                          ncol = length(V(graph)))
      bit_est <- EstLinearAP(graph_adj, r = 15)
      t <- c(bit_est$Alpha, bit_est$Beta, bit_est$delta_in, bit_est$delta_out)
      est <- c(S, assortativity_degree(graph, directed = "TRUE"), transitivity(graph, "average")
               ,gini.index(degree(graph, mode = "in"))$statistic, gini.index(degree(graph, mode = "out"))$statistic, bit_est, CalcExp(theta = t))
      X <- rbind(X, est)
    }
    else {
      graph <- graph.data.frame(edge, directed = TRUE)
      graph_adj <- matrix(as_adjacency_matrix(graph), nrow = length(V(graph)), 
                          ncol = length(V(graph)))
      bit_est <- EstLinearAP(graph_adj, r = 15)
      t <- c(bit_est$Alpha, bit_est$Beta, bit_est$delta_in, bit_est$delta_out)
      est <- c(dim(edge)[1], assortativity_degree(graph, directed = "TRUE"), transitivity(graph, "average"), gini.index(degree(graph, mode = "in"))$statistic,
               gini.index(degree(graph, mode = "out"))$statistic, bit_est, CalcExp(theta = t))
      X <- rbind(X, est)      
    }
  }
  X <- X[-1, ]
  return(X)
}


part_dynamics <- function(edge, K){
  N <- dim(edge)[1] %/% K
  M <- K + 1
  X <- rep(0, 2)
  Q <- 1
  for (j in 1:M){
    if (j < M){
      S <- as.integer(j*N); print(S)
      graph <- graph.data.frame(edge[Q:S, ], directed = TRUE)
      graph_adj <- matrix(as_adjacency_matrix(graph), nrow = length(V(graph)), 
                          ncol = length(V(graph)))
      bit_est <- EstLinearAP(graph_adj, r = 6)
      t <- c(bit_est$Alpha, bit_est$Beta, bit_est$delta_in, bit_est$delta_out)
      est <- c(S, assortativity_degree(graph, directed = "TRUE"), transitivity(graph, "average")
               ,gini.index(degree(graph, mode = "in"))$statistic, gini.index(degree(graph, mode = "out"))$statistic, bit_est, CalcExp(theta = t))
      X <- rbind(X, est)
      Q <- S + 1
    }
    else {
      graph <- graph.data.frame(edge[Q:dim(edge)[1], ], directed = TRUE)
      graph_adj <- matrix(as_adjacency_matrix(graph), nrow = length(V(graph)), 
                          ncol = length(V(graph)))
      bit_est <- EstLinearAP(graph_adj, r = 6)
      t <- c(bit_est$Alpha, bit_est$Beta, bit_est$delta_in, bit_est$delta_out)
      est <- c(dim(edge)[1], assortativity_degree(graph, directed = "TRUE"), transitivity(graph, "average"), gini.index(degree(graph, mode = "in"))$statistic,
               gini.index(degree(graph, mode = "out"))$statistic, bit_est, CalcExp(theta = t))
      X <- rbind(X, est)      
    }
  }
  X <- X[-1, ]
  return(X)
}



if (FALSE){
GetAge <- function(V, edge, direction){
	m <- matrix(0, )
	for (k in 1:length(V)){
		for (i in 1:2){
			for (j in 1:dim(edge)[1]){
				if (i == 1){
					v <- V[k]
					
				}
				if (i == 2){
					v <- V[k]
				}
				
				
			}
		}
	}
}
}