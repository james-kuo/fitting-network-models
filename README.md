Linear Preferential Attachment models are frequently used to theoretically explain why scale-free networks exist, or networks which have power-law as its degree distribution, likelihood-based approaches are seldom done to test how well model explains data. I fitted such a model using a novel likelihood-based method developed in the literature on a bit-coin network in R and devised a set of informal simulation studies to examine how well model matches up with data. I developed the code from bottom-up for this project. 

`Final_Project.pdf` is the paper. `Kuo_Final_Presentation.pdf` are the slides for an earlier version of the paper. 

`Functions.R`, `Analysis.R`, and `LinearAP.R` are the main R codes which process the data, simulate the network, and estimate the model. 

`Bit.RData`, `SampleSim.RData`, `updated_graphs.RData` and `dynamics.RData` are the processed datasets. `soc-sign-bitcoinotc.csv` is the raw dataset. 
