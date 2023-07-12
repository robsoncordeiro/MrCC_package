# The Method Halite for Correlation Clustering

<p>
    The algorithm Halite is a fast and scalable density-based clustering algorithm for moderate-to-high-dimensionality data able to analyze large collections of complex data elements. It creates a multi-dimensional grid all over the data space and counts the number of points lying at each hyper-cubic cell provided by the grid. A hyper-quad-tree-like structure, called the Counting-tree, is used to store the counts. The tree is thereafter submitted to a filtering process able to identify regions that are, in a statistical sense, denser than its neighboring regions regarding at least one dimension, which leads to the final clustering result. The algorithm is fast and it has linear or quasi-linear time and space complexity regarding both the data size and the dimensionality.
</p>
<p>
    <b>Remark:</b> A first implementation of Halite was initially named as the method MrCC (after Multi-resolution Correlation Clustering) in an earlier <a href="http://doi.ieeecomputersociety.org/10.1109/ICDE.2010.5447924"><b>Conference Publication</b></a> of this work. Later, it was renamed to Halite for clarity, since several improvements on the initial implementation were included into a <a href="http://doi.ieeecomputersociety.org/10.1109/TKDE.2011.176"><b>Journal Publication</b></a>.
</p>
<p>
    <a href="http://doi.ieeecomputersociety.org/10.1109/ICDE.2010.5447924"><b>Conference Publication</b></a><br/>
    <a href="http://doi.ieeecomputersociety.org/10.1109/TKDE.2011.176"><b>Journal Publication</b></a>
</p>

<h3>MrCC - basic information</h3>

*******************************************************************

input data format:
 
axis_1 axis_2 axis_3 ... axis_d groundTruthCluster
axis_1 axis_2 axis_3 ... axis_d groundTruthCluster
				.
				.
				.
axis_1 axis_2 axis_3 ... axis_d groundTruthCluster				

Example: databases\12.dat

Obs.: - the groundTruthCluster data is not used by MrCC.
      - the file databases\12d.sub contains the ground truth WRT the clusters' subspaces for our example dataset, but it's also not used by MrCC.
      - you may use this ground truth information to evaluate the quality of MrCC's results on our example dataset.

*******************************************************************

output data format: 

ClusterResult1 relevanceOfAxis_1 relevanceOfAxis_2 ... relevanceOfAxis_d
ClusterResult2 relevanceOfAxis_1 relevanceOfAxis_2 ... relevanceOfAxis_d
				.
				.
				.
ClusterResultk relevanceOfAxis_1 relevanceOfAxis_2 ... relevanceOfAxis_d
LABELING
PointId clusterId
PointId clusterId
	.
	.
	.
PointId clusterId

Obs.: - relevanceOfAxis_j = 0 means that axis j is irrelevant to the corresponding cluster, while relevanceOfAxis_j = 1 means that axis j is relevant to this cluster.
      - PointId identifies one point by referring to the line in which this point is found in the input data file. For soft clustering, each PointId value can be related to more than one clusterId value.
      - clusterId = 0 means that the corresponding point is an outlier.

Example: results\result12d.dat

*******************************************************************
compiling MrCC:

First, you must install two third-part software:
   - Oracle Berkeley DB: "http://www.oracle.com/technetwork/database/berkeleydb/overview/index.html"
   - OpenCV: "http://opencv.willowgarage.com/"

Then, compile the code using any standard c++ compiler.

*******************************************************************

running MrCC:

MrCC \alpha H hardClustering initialLevel

Obs.: - the input/output specs are defined in "arboretum/ioSpecs.h".
      - the default value for \alpha is 1e-10.
      - the default value for H is 4.
      - hardClustering = 1 means that the result will be a dataset partition (one point belongs to at most one cluster).
      - hardClustering = 0 means that the algorithm will do soft-clustering (one point can belong to more than one cluster).
      - the default value for initialLevel is 1.
      - distinct forms of memory management are possible by tuning the "memory" parameter in the MrCC.cpp file. "memory==0" means that both the dataset and the tree will be put in main memory. "memory==1" means that only the tree will be put in main memory. "memory==2" means that neither the database nor the tree will be put in main memory.

Example: make demo (Linux or Mac OS), runExample.bat (Windows).
