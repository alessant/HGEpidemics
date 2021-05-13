# Epidemic Dynamics via Time-Varying Hypergraphs

Experiment related to the paper *A Design-Methodology for Epidemic Dynamics via Time-Varying Hypergraphs*, A. Antelmi, G. Cordasco, C. Spagnuolo and V. Scarano.

**Reproduce the experiments**
If you want to reproduce the results of the paper:
* Download the Tokyo dataset at the following [link](https://sites.google.com/site/yangdingqi/home/foursquare-dataset?authuser=0);
* Run the script `src/experiments/spreading/AAMAS20/AAMAS20_spreading_experiment.jl`.

**Abstract**
In epidemiology science, the importance to explore innovative modeling tools for acutely analyzing epidemic diffusion is turning into a big challenge considering the myriad of real-world aspects to capture. Typically, equation-based models, such as SIS and SIR, are used to study the propagation of diseases over a population. Improved approaches also include human-mobility patterns as network information to describe contacts among individuals. However, there still is the need to incorporate in these models information about different types of contagion, geographical information, humans habits, and environmental properties. In this paper, we propose a novel approach that takes into account: 1. direct and indirect epidemic contagion pathways to explore the dynamics of the epidemic, 2. the times of possible contagions, and 3. human-mobility patterns. We combine these three features exploiting time-varying hypergraphs, and we embed this model into a design-methodology for agent-based models (ABMs), able to improve the correctness in the epidemic estimations of classical contact-network approaches. We further describe a diffusion algorithm suitable for our design-methodology and adaptable to the peculiarities of any disease spreading policies and/or models. Finally, we tested our methodology by developing an ABM, realizing the SIS epidemic compartmental model, for simulating an epidemic propagation over a population of individuals. We experimented the model using real user-mobility data from the location-based social networkFoursquare, and we demonstrated the high-impact of temporal direct and indirect contagion pathways.


> @inproceedings{10.5555/3398761.3398774,
author = {Antelmi, Alessia and Cordasco, Gennaro and Spagnuolo, Carmine and Scarano, Vittorio},
title = {A Design-Methodology for Epidemic Dynamics via Time-Varying Hypergraphs},
year = {2020},
isbn = {9781450375184},
publisher = {International Foundation for Autonomous Agents and Multiagent Systems},
address = {Richland, SC},
booktitle = {Proceedings of the 19th International Conference on Autonomous Agents and MultiAgent Systems},
pages = {61–69},
numpages = {9},
location = {Auckland, New Zealand},
series = {AAMAS '20}
}


