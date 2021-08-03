# Epidemic Dynamics via Time-Varying Hypergraphs

If you want to reproduce the experiments related to the paper *A Design-Methodology for Epidemic Dynamics via Time-Varying Hypergraphs*, A. Antelmi, G. Cordasco, V. Scarano, and C. Spagnuolo, follow these instructions:

**Reproduce the experiments**
* Download the Tokyo dataset at the following [link](https://sites.google.com/site/yangdingqi/home/foursquare-dataset?authuser=0);
* Run the script `src/experiments/spreading/AAMAS20/AAMAS20_spreading_experiment.jl`.


If you want to reproduce the experiments related to the paper *Modeling and Evaluating Epidemic Control Strategies with High-Order Temporal Networks* currently submitted to the Jounal IEEE Access, A. Antelmi, G. Cordasco, V. Scarano, and C. Spagnuolo, follow these instructions:

**Reproduce the experiments**
* [Exp1 - NPIs application] Run the script `src/experiments/NPIs/spreading_experiment_with_NPIs.jl`, specifying as input parameter the configuration file (e.g., see `src/experiments/NPIs/BLE/configs/npis/nsga/data_params.json`).
* [Exp2 - NSGA] Run the script `src/experiments/NSGA/NPIs_optimization.jl`, specifying as input parameter the configuration file (e.g., see `src/experiments/NSGA/BLE/configs/ble_params.json`).

