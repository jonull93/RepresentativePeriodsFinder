# Methods

The optimisation and clustering methods which are implemented here are adapted from in [Long term storage in generation expansion planning models with a reduced temporal scope by Gonzato, S., Bruninx, K., Delarue, E., 2021.](https://www.mech.kuleuven.be/en/tme/research/energy-systems-integration-modeling/pdf-publications/wp-esim2021-1).

## Optimisation

The optimisation problem which is solved can be written as follows:



## Clustering

Ward's hierarchical clustering method is used. This method is well described [here](https://www.mech.kuleuven.be/en/tme/research/energy-systems-integration-modeling/pdf-publications/wp-esim2021-1) and [here](https://www.researchgate.net/publication/325459549_Chronological_Time-Period_Clustering_for_Optimal_Capacity_Expansion_Planning_With_Storage).

One particularity with the algorithm is the ability to include mandatory periods. This is done by fixing the dissimilarity for the mandatory period with respect to the all other periods to `Inf`. This ensures that it is not merged. It is not possible to assign this period any other weight than one however. 

Another particularity is the ability to choose the dissimilarity function (ordering error) freely (see [Selecting representative periods](@ref)).
