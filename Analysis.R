
source("/Users/MacUser/Desktop/Network/Final Project/LinearAP.R")
source("/Users/MacUser/Desktop/Network/Final Project/Functions.R")

library(igraph)
set.seed(134)

init_g <- sample_pa(20, directed = TRUE)
params <- c(0.4, 0.2, 1, 1)

sample_adj <- SimLinearAP(init = init_g, theta = params, n = 1000)
sample_g <- graph_from_adjacency_matrix(sample_adj$adj, mode = "directed", weighted = NULL, diag = TRUE)

dd.sample_g <- degree.distribution(sample_g, mode= "in")
d <- (0:(max(degree(sample_g, mode = "in"))-1))
ind <- (dd.sample_g!=0)
plot(d[ind], dd.sample_g[ind], log = "xy", col = "blue", 
     xlab = "Log-Degree", ylab = "Log-Frequency",  
     main = "Log-Log Degree Distribution")

in_degree <- colSums(sample_adj$adj)
out_degree <- rowSums(sample_adj$adj)

plot(in_degree, out_degree)

#Bitcoin 

library(igraph);
bit_edge <- read.csv("/Users/MacUser/Desktop/Network/Final Project/soc-sign-bitcoinotc.csv", header = FALSE)
bit_edge <- bit_edge[,-4]
bit_edge <- bit_edge[, -3]
bit_graph <- graph.data.frame(bit_edge, directed=TRUE)
bit_adj <- matrix(as_adjacency_matrix(bit_graph), nrow = 5881, ncol = 5881)

bit_est <- EstLinearAP(adj = bit_adj, r = 10, normalize = 1)
bit_theta <- c(bit_est$Alpha, bit_est$Beta, bit_est$delta_in, bit_est$delta_out)

sim_bit <- SimLinearAP(init = init_g, theta = bit_theta, n = 35592)
sim_bit_g <- graph_from_adjacency_matrix(sim_bit$adj, mode = "directed", weighted = NULL, diag = TRUE)
V(sim_bit_g)$name <- V(sim_bit_g)

par(mfrow=c(1,2)) 
par(mar = c(3,2,1,2))
PlotDegree(graph1 = sim_bit_g, graph2 = bit_graph, direction = "in")
PlotDegree(graph1 = sim_bit_g, graph2 = bit_graph, direction = "out")

in_degree_real <- colSums(bit_adj)
out_degree_real <- rowSums(bit_adj)

in_degree_sim <- colSums(sim_bit$adj)
out_degree_sim <- rowSums(sim_bit$adj)

par(mfrow=c(1,2)) 
par(mar = c(2,2,2,2))
plot(in_degree_real, out_degree_real, log = "xy", col = "red", cex = 0.1, cex.lab = 0.1, main = 
       "Joint In-Out Degree of Real Network (Log Scale)", xlab = "Log(in-degree)", ylab = "Log(out-degree)")
plot(in_degree_sim, out_degree_sim, log = "xy", col ="blue", cex = 0.1, cex.lab = 0.1,
     main = "Joint In-Out Degree of Simulated Netowork (Log Scale)", xlab = "Log(in-degree)", ylab = "Log(out-degree)" )


est_sim_bit <- EstLinearAP(adj = sim_bit$adj, r = 10, normalize = 1)
est_sim_bit <- c(est_sim_bit$Alpha, est_sim_bit$Beta, est_sim_bit$delta_in, est_sim_bit$delta_out)
CalcExp(est_sim_bit)

save(bit_graph, bit_adj, sim_bit, sim_bit_g, file="/Users/MacUser/Desktop/Network/Final Project/Bit.RData")

###############
#Evolution of Joint degree distribution 

test <- Partition(matrix = bit_edge, K = 4)
test0 <- test[1]
test1 <- rbind(test$`(0,1e+04]`, test$`(1e+04,2e+04]`)
test2 <- rbind(test1, test$`(2e+04,3e+04]`)

test0 <- graph.data.frame(test0, directed=TRUE)
test0 <- matrix(as_adjacency_matrix(test0), nrow = dim(as_adjacency_matrix(test0))[1], 
                ncol = dim(as_adjacency_matrix(test0))[1])

test1 <- graph.data.frame(test1, directed=TRUE)
test1 <- matrix(as_adjacency_matrix(test1), nrow = dim(as_adjacency_matrix(test1))[1],
                ncol = dim(as_adjacency_matrix(test1))[1])

test2 <- graph.data.frame(test2, directed=TRUE)
test2 <- matrix(as_adjacency_matrix(test2), nrow = dim(as_adjacency_matrix(test2))[1],
                ncol = dim(as_adjacency_matrix(test2))[1])


real0_in <- colSums(test0); real0_out <- rowSums(test0)
real1_in <- colSums(test1); real1_out <- rowSums(test1)
real2_in <- colSums(test2); real2_out <- rowSums(test2)

par(mfrow=c(1,2)) 
par(mar = c(2,2,2,2))
plot(in_degree_real, out_degree_real, log = "xy", pch = 16, col = "red", cex = 0.2, cex.lab = 0.2, main = 
       "Joint In-Out Degree of Real Network (Log Scale)", xlab = "Log(in-degree)", ylab = "Log(out-degree)")
points(real0_in, real0_out,log = "xy",   pch = 10, cex = 0.2, cex.lab = 0.2, col = "gray25")
#points(real1_in, real1_out, log = "xy", pch = 4, cex = 0.2, cex.lab = 0.2, col = "red")
#points(real2_in, real2_out,log = "xy",  pch = 6, cex = 0.2, cex.lab = 0.2, col = "gray69")
legend("topleft", legend=c("First 10,000 Edges", "All 35,592 Edges"), col = c("gray25", "red"), pch=c(10, 16), cex=0.8, bty = "n")


test_sim <- Partition(matrix = sim_bit$edge_list, K = 4)
test0_sim <- graph.data.frame(test_sim[1], directed=TRUE)
test0_sim <- matrix(as_adjacency_matrix(test0_sim), nrow = dim(as_adjacency_matrix(test0_sim))[1], 
                ncol = dim(as_adjacency_matrix(test0_sim))[1])

sim0_in <- colSums(test0_sim); sim0_out <- rowSums(test0_sim)

par(mar = c(2,2,2,2))
plot(in_degree_sim, out_degree_sim, log = "xy", pch = 16, col = "blue", cex = 0.2, cex.lab = 0.2, main = 
       "Joint In-Out Degree of Simulated Network (Log Scale)", xlab = "Log(in-degree)", ylab = "Log(out-degree)")
points(sim0_in, sim0_out,log = "xy",   pch = 10, cex = 0.2, cex.lab = 0.2, col = "gray25")
#points(real1_in, real1_out, log = "xy", pch = 4, cex = 0.2, cex.lab = 0.2, col = "red")
#points(real2_in, real2_out,log = "xy",  pch = 6, cex = 0.2, cex.lab = 0.2, col = "gray69")
legend("topleft", legend=c("First 10,000 Edges", "All 35,592 Edges"), col = c("gray25", "blue"), pch=c(10, 16), cex=0.8, bty = "n")

PlotDegree(graph1 = bit_graph, graph2 = graph_from_adjacency_matrix(test0, mode = "directed",
                  weighted = NULL, diag = TRUE ), direction = "in")

par(mfrow=c(1,2)) 
par(mar = c(2,2,2,2))
sort_plot(g = bit_graph)
sort_plot(g = sim_bit_g)

# Age vs. degree
out_t_real <- match(as.integer(V(bit_graph)$name), bit_edge[, 1])
in_t_real <-  match(as.integer(V(bit_graph)$name), bit_edge[, 2])
V(bit_graph)$out_age <- dim(bit_edge)[1] -  out_t_real
V(bit_graph)$in_age <- dim(bit_edge)[1] - in_t_real

out_t_sim <- match(as.integer(V(sim_bit_g)$name), sim_bit$edge_list[, 1])
in_t_sim <-  match(as.integer(V(sim_bit_g)$name), sim_bit$edge_list[, 2])
V(sim_bit_g)$out_age <- dim(sim_bit$edge_list)[1] -  out_t_sim
V(sim_bit_g)$in_age <- dim(sim_bit$edge_list)[1] - in_t_sim

par(mfrow=c(1,2)) 
par(mar = c(3,3,3,3))
plot(y = in_degree_real, x = V(bit_graph)$in_age, col = "red", log = "xy", cex = 0.1, cex.lab = 0.1, 
     xlab = "", ylab = "", main = "Node In-Degree Age vs. In-Degree (Log)", cex.main = 0.7, cex.axis = 0.5)
points(y = in_degree_sim, x = V(sim_bit_g)$in_age, log = "xy", col = "blue", cex = 0.1, cex.lab = 0.1)
title(ylab="In-Degree", xlab = "In-Degree Age", line=1.5, cex.lab=0.7)
legend("topleft", legend=c("Simulated", "Real"), col = c("blue", "red"), pch=16, cex=0.3, bty = "n")

plot(y = out_degree_real, x = V(bit_graph)$out_age, col = "red", log = "xy", cex = 0.1, cex.lab = 0.1, 
     xlab = "", ylab = "", main = "Node Out-Degree Age vs. Out-Degree (Log)", cex.main = 0.7, cex.axis = 0.5)
points(y = out_degree_sim, x = V(sim_bit_g)$out_age, log = "xy", col = "blue", cex = 0.1, cex.lab = 0.1)
title(ylab="Out-Degree", xlab = "In-Degree Age", line=1.5, cex.lab=0.7)
legend("topleft", legend=c("Simulated", "Real"), col = c("blue", "red"), pch=16, cex=0.3, bty = "n")

plot(y = in_degree_sim, x = V(sim_bit_g)$in_age, col = "blue",  cex = 0.1, cex.lab = 0.1, 
     xlab = "", ylab = "", main = "Node Out-Degree Age vs. Out-Degree (Log)", cex.main = 0.7, cex.axis = 0.5)
points(y = out_degree_sim, x = V(sim_bit_g)$out_age,col = "red", cex = 0.1, cex.lab = 0.1)

plot(y = real$in_degree, x = real$total_age, log="xy", col = "red",  cex = 0.1, cex.lab = 0.1, 
     xlab = "", ylab = "", main = "Node Out-Degree Age vs. Out-Degree (Log)", cex.main = 0.7, cex.axis = 0.5)
points(y = sim$in_degree, x = sim$total_age, log = "xy", col = "blue", cex = 0.1, cex.lab = 0.1)

# Estimate the age versus degree relationship in the simulated network, and use it to 
# predict the real netowkr and get the difference 
real <- data.frame(V(bit_graph)$out_age)
real[, "out_age"] <- data.frame(V(bit_graph)$out_age)
real[,1] <- NULL
real[,"out_degree"] <- out_degree_real
real[, "in_age"] <- V(bit_graph)$in_age
real[, "in_degree"] <- in_degree_real
real[, "total_age"] <- apply(real[, c(1,3)], 1, FUN=max)

sim <- data.frame(V(sim_bit_g)$out_age)
sim[, "out_age"] <- data.frame(V(sim_bit_g)$out_age)
sim[,1] <- NULL
sim[,"out_degree"] <- out_degree_sim
sim[, "in_age"] <- V(sim_bit_g)$in_age
sim[, "in_degree"] <- in_degree_sim
sim[, "total_age"] <- apply(sim[, c(1,3)], 1, FUN=max)
save(sim_bit_g, bit_graph, real, sim, file="/Users/MacUser/Desktop/Network/Final Project/updated_graphs.RData")

m_in_sim <- lm(I(log(in_degree))~ I(log(in_age)), data = sim)
m_out_sim <- lm(I(log(out_degree))~ I(log(out_age)), data = sim)

mse_in <- mean((real$in_degree- exp(predict(m_in_sim, newdata = real)))^2, na.rm = TRUE )
mse_out <- mean((real$out_degree- exp(predict(m_out_sim, newdata = real)))^2, na.rm = TRUE )

lm(I(log(out_degree))~ I(log(total_age)), data = real)
lm(I(log(out_degree))~ I(log(out_age)), data = real)

## Index/Parameter Dynamics 


source("/Users/MacUser/Desktop/Network/Final Project/Functions.R")
real_bit_part <- Partition(matrix = bit_edge, K = 26)
sim_bit_part <- Partition(matrix = sim_bit$edge_list, K = 20)

library(igraph); library(plyr)


real_par <- par_dynamics(partition = real_bit_part )
sim_par <- par_dynamics(partition = sim_bit_part)


real_bit_part1 <- Partition(matrix = bit_edge, K = 10)
sim_bit_part1 <- Partition(matrix = sim_bit$edge_list, K = 10)
real_par1 <- par_dynamics(partition = real_bit_part1 )
sim_par1 <- par_dynamics(partition = sim_bit_part1)

real_dynamics <- dynamics(edge = bit_edge, K = 30)
sim_dynamics <- dynamics(edge = as.matrix.data.frame(sim_bit$edge_list), K = 30)

part_real_dynamics <- part_dynamics(edge = bit_edge, K = 30)
part_sim_dynmaics <- part_dynamics(edge = as.matrix.data.frame(sim_bit$edge_list), K = 30)

real_gini <- real_dynamics[, 3]
sim_gini <- sim_dynamics[,3]

par(mfrow=c(1,2)) 
par(mar = c(3,3,3,3))
plot(x = degree(bit_graph, mode = "out"),y = transitivity(bit_graph, "local"), log = "xy", cex = 0.1)
plot(x = degree(sim_bit_g, mode = "out"), y = transitivity(sim_bit_g, "local"), log = "xy", cex = 0.1)

par(mfrow=c(1,2)) 
par(mar = c(3,3,3,3))
plot(x = real_dynamics[, 1], y = real_dynamics[, 4],
     ylim = c(min(as.numeric(real_dynamics[,4]), as.numeric(sim_dynamics[,4]) ), 
              max(as.numeric(real_dynamics[,4]), as.numeric(sim_dynamics[,4]) )), type = "l")
lines(x = sim_dynamics[, 1], y = sim_dynamics[, 4], type = "l", col = "blue")

plot(x = real_dynamics[, 1], y = real_dynamics[, 5],
     ylim = c(min(as.numeric(real_dynamics[,5]), as.numeric(sim_dynamics[,5]) ), 
              max(as.numeric(real_dynamics[,5]), as.numeric(sim_dynamics[,5]) )), type = "l")
lines(x = sim_dynamics[, 1], y = sim_dynamics[, 5], type = "l", col = "blue")

plot(x = real_dynamics[, 1], y = real_dynamics[, "delta_out"], type = "l")
lines(x = sim_dynamics[, 1], y = sim_dynamics[, "delta_out"], type = "l", col = "blue")

save(real_dynamics, sim_dynamics, part_real_dynamics, 
     part_sim_dynmaics, file="/Users/MacUser/Desktop/Network/Final Project/dynamics.RData")

plot(x = real_dynamics[, 1], y = real_dynamics[, 12],
     ylim = c(min(as.numeric(real_dynamics[,12]), as.numeric(sim_dynamics[,12]) ), 
              max(as.numeric(real_dynamics[,12]), as.numeric(sim_dynamics[,12]) )), type = "l")
lines(x = sim_dynamics[, 1], y = sim_dynamics[, 12], type = "l", col = "blue")

plot(x = real_dynamics[, 1], y = real_dynamics[, 11],
     ylim = c(min(as.numeric(real_dynamics[,11]), as.numeric(sim_dynamics[,11]) ), 
              max(as.numeric(real_dynamics[,11]), as.numeric(sim_dynamics[,11]) )), type = "l")
lines(x = sim_dynamics[, 1], y = sim_dynamics[, 11], type = "l", col = "blue")


## Gini Dynamics
real_gini_part <- gini_dynamics(partition = real_bit_part, direction = "out")
sim_gini_part <- gini_dynamics(partition = sim_bit_part, direction = "out")

plot(x = 1:18, y = sim_gini_part, col = "blue", type = "l")
points(x = 1:18, y = real_gini_part, col = "red", type = "l")

plot(x = 1:21, y = real_gini, col ="red", type = "l")
points(x = 1:21, y = sim_gini, col = "blue", type = "l")




## Cluster inspection
sim_cluster <- cluster_infomap(sim_bit_g)
real_cluster <- cluster_infomap(bit_graph)

hist(sizes(sim_cluster), breaks = 3000)
hist(sizes(real_cluster), breaks = 3000)

d <- (0:(max(degree(french, mode = "in"))-1))
ind <- (dd.sample_g!=0)
plot(d[ind], dd.sample_g[ind], log = "xy", col = "blue", 
     xlab = "Log-Degree", ylab = "Log-Frequency",  
     main = "Log-Log Degree Distribution")


plot(x = degree(bit_graph, mode = "in"), y = knn(bit_graph, vids = V(bit_graph))$knn, log = "xy",
     cex = 0.1)
plot(x = degree(simplify(sim_bit_g), mode = "in"), y = k2$knn, log = "xy", cex = 0.1)

## Use this if neccesary 
k <- knn(bit_graph, vids = V(bit_graph))
k2 <- knn(simplify(sim_bit_g))
par(mfrow=c(1,1)) 
par(mar = c(2.5,2.5,2.5,2.5))
plot(x = 1:max(degree(bit_graph, mode = "total")), y = k$knnk, cex = 0.2, col = "red",
     xlim = c(0,600), ylim = c(0, 300), cex.lab = 0.1, cex.axis = 0.5, 
     cex.main = 0.5, main = "Average Degree of Nearest Neighbors as a Function of Degree")
points(x = 1:max(degree(simplify(sim_bit_g), mode = "total")), y = k2$knnk, cex = 0.2, col = "blue")
abline(a = 0, b =1, lty = 3)
title(xlab="Total Degree", ylab = "Average Total Degree of Nearest Neighbors", line=1.5, cex.lab=0.5)
legend("topright", legend=c("Simulated", "Real"), col = c("blue", "red"), pch=16, cex=0.5, bty = "n")




