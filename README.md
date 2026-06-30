# Colon Crypt Stem Cell Lineage Simulation

## Introduction
Human intestinal stem cells (hISCs) at the base of each crypt are responsible for regenerating the crypt lining as often as every three to five days in the 10 million crypts covering the inner surface of the human colon.  Residing in a stem cell niche, hISCs divide to renew the stem cells (SCs) in the niche and to populate the transit-amplifying (TA) cell pool to form the crypt lining.  TA cells proliferate, differentiate as they migrate upward, and are shed into the colon lumen after reaching the crypt orifice.  It has been suggested that each colon crypt, which consists of ~1,500 cells, becomes a single stem cell lineage monoclone in early life as a result of neutral drift. We simulated a simplified model of colon crypt development to explore various developmental scenarios and their variant allele frequency profiles over successive generations.

## Simulation details
We model the cell proliferation dynamics of two compartments, the stem cell niche (`SCniche`) where active SCs reside, and the remaining colon crypt (`Crypt`) containing TA cells. 

In the `SCniche`, the number of SCs are set at a fixed **_N_**. Over every SC cell cycle generation, each SC divides into two daughter SCs (total **_2N_** SCs), with each daughter SC inheriting all the mutations of the parent SC, plus accumulating new mutations at a rate of **_M_** mutations per genome per cell division. Half of the **_2N_** daughter SCs will then be randomly selected to remain in the `SCniche`, while the other half will be pushed into the `Crypt` as TA cells. 

The `Crypt` is modelled with a carrying capacity of 2048 TA cells where the positional location of each cell is recorded and trackable (i.e. stored in the form of a horizontal data array of 2048 slots). The **_N_** SCs randomly chosen from the `SCniche` to be pushed into the `Crypt` as TA cells will enter and occupy one end of the array (i.e. the left-most **_N_** slots of the array), pushing along all other extant TA cells in the `Crypt` further along (i.e. rightwards) down the array. With every TA cell cycle/doubling time, each TA cell divides into two daughter TA cells _in situ_, now occupying two adjacent slots and therefore also pushing along the TA cells on the right further down the array. Each daughter TA cell inherits all the mutations of the parent TA cell, and also accumulates new mutations at a similar rate of **_M_** mutations per genome per cell division. Therefore, with each SC and TA cell cycle generation, there will be TA cells that are pushed rightwards and out of the 2048 available array slots, representing the shedding of TA cells out of the `Crypt`.      
        
#### Additional model assumptions ####
- Cell cycle proliferation of SC and TA cells are synchronized, with the cell cycle/doubling time of TA cells (~12 hrs) being half the cell cycle/doubling time of SCs (~24hrs).
- Mutations accumulate at unique locations, thus eliminating back-mutations and reduced observations of mutation counts.

We ran the simulations with **_N_** from 4-16 and **_M_** from 1-3, each over 500 generations. Mutation allele frequencies are calculated every 10 generations. Results of these simulations are presented in the form of allele frequency distribution plots in an R Shiny App (see next section).   

## Explore simulation results
1. For an interactive exploration of the simulation results, please click [here](https://eddieloh-usc.github.io/eddieloh-usc-pages/ColonCryptSCLineageSymmetricSim/index.html).  
2. If the link above is not accessible, the interactive app can also be run from your local R/RStudio software with the code below.  
      ```
      # Install the Shiny package (if you haven't done so)
      install.packages("shiny")

      # Load the package and run the app
      library("shiny")
      runUrl("https://github.com/eddieloh-usc/ColonCryptSim/raw/refs/heads/main/ShinyApp.zip")
      ```

## Citation
[to be filled...]
