Linear Preferential Attachment models are frequently used to theoretically explain why scale-free networks exist, or networks which have power-law as its degree distribution, likelihood-based approaches are seldom done to test how well model explains data. I fitted such a model using a novel likelihood-based method developed in the literature on a bit-coin network in R and devised a set of informal simulation studies to examine how well model matches up with data. I developed the code from bottom-up for this project. 

`Final_Project.pdf` is the paper. `Kuo_Final_Presentation.pdf` are the slides for an earlier version of the paper. 

`Functions.R`, `Analysis.R`, and `LinearAP.R` are the main R codes which process the data, simulate the network, and estimate the model. 

`Bit.RData`, `SampleSim.RData`, `updated_graphs.RData` and `dynamics.RData` are the processed datasets. `soc-sign-bitcoinotc.csv` is the raw dataset. 

## Some Pictures

### Degree Distribution
A key aspect of the project is to see if degree distribution of the real network is well-captured by the model. The model hypothesized that both the in-degree and out-degree distributions of the real network will follow power law, meaning there will be a straight line in the log-log plot. Shown below are log-log plots of the real and simulated network. Though the degree distributions of the simulated network seem to visually match that of the real network, a more rigirous statistical test (Kolmogorov-Smirnov test) shows that the two distributions are different, so the model fails in this aspect

![](https://github.com/james-kuo/fitting-network-models/blob/master/degree_distribution.png)


### A Picture of the Network
This is what the model says how the bitcoin network should look like. 

![](https://github.com/james-kuo/fitting-network-models/blob/master/SimLAP.png)
