# Disruptive RIS for Enhancing Key Generation and Secret Transmission in Low-Entropy Environments
This repository contains the MATLAB codes for all the figures in the paper authored by myself Hibatallah Alwazani and by Dr. Anas Chaaban, linked here "[https://arxiv.org/pdf/2409.15303](https://arxiv.org/abs/2409.15303)"

## Abstract
Key generation, a pillar in physical-layer security( PLS), is the process of the exchanging signals from two legitimate users (Alice and Bob) to extract a common key from the random, common channels. The drawback of extracting
keys from wireless channels is the ample dependence on the dynamicity and fluctuations of the radio channel, rendering the key vulnerable to estimation by Eve (an illegitimate user) in low-entropy environments because of insufficient randomness. Added to that, the lack of channel fluctuations lower the secret key rate (SKR) defined as the number of bits of key generated per channel use. In this work, we aim to address this challenge by using a reconfigurable intelligent surface (RIS) to produce random phases at certain, carefully curated intervals such that it disrupts the channel in low-entropy environments. We propose an RIS assisted key generation protocol, study its performance, and compare with benchmarks to observe the benefit of using an RIS while considering various important metrics such as key mismatch rate and secret key throughput. Furthermore, we characterize a scaling law as a function of the rate of change of RIS phase switching for the average secret information rate under this protocol. Then, we use both the key throughput and information rate to optimize the overall secrecy rate. Simulations are made to validate our theoretical findings and effectiveness of the proposed scheme showing an improvement in performance when an RIS is deployed.

## Code Information
There are three function files used within the figure files, one is called Scenario.m where the hyperparameters of the wireless system is fixed such as pathloss parameters and locations. Another is the channels.m where all the channel realizations are set, and Match_probability.m which is used to calculation the key match probability.



*Go ahead and play with this code, but kindly if you use it in your research, reference our paper in your work.*
