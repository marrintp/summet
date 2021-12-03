# Description #
The Subspace Mean and Median Evaluation Toolkit (SuMMET) is a Matlab-based software package that is meant as a companion for the paper:  

"Finding the Subspace Mean or Median to Fit Your Need," that has been accepted for the 2014 IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR 2014).

The author of this software package is Tim Marrinan, written while he was a Ph.D. student working under the direction of Michael Kirby.  

The toolkit has three main components: (1) core Matlab functions that compute the Karcher mean, the L2-median, the extrinsic manifold mean, and the flag mean for collections of linear subspaces, (2) experimental data in Matlab format made up of short video clips depicting "micro-actions" taken from longer videos associated with DARPA's Mind's Eye program (and associated action labels), and (3) example Matlab scripts that use the core functions and experimental data to produce the results/figures from the associated paper.

The intention of this software package is that researchers will be able to compare the properties of the subspace averages in a straightforward way to develop an intuition about which method is most appropriate for a specific task. Additionally, it provides computational results that can be compared to those of subspace averages not included here when necessary.


# Contents #

The matlab package includes the following directories and files:

\
01. README.md  
02. LICENSE.txt

\src\
03.  ExManifoldMean.m  
04.  FlagMean.m  
05.  GrDist.m  
06.  GrExp.m  
07.  GrLog.m  
08.  KarcherMean.m  
09.  L2Median.m

\src\helper\
10.  aboxplot.m  
11.  ClusterPurity.m  
12.  ExemplarSelection.m  
13.  fig.m  
14.  GrKmeans.m  
15.  LabelCounter.m  
16.  Plot2DMeans.m  
17.  quartile.m

\examples\
18. CVPR_KmeansScript.m  
19. CVPR_PrototypeSelectionScript.m  
20. CVPR_SyntheticIllustrationScript.m  


# Abstract #

Many computer vision algorithms employ subspace models to represent data. Many of these approaches benefit from the ability to create an average or prototype for a
set of subspaces. The most popular method in these situations is the Karcher mean, also known as the Riemannian center of mass. The prevalence of the Karcher mean may
lead some to assume that it provides the best average in all scenarios. However, other subspace averages that appear less frequently in the literature may be more appropriate for certain tasks. The extrinsic manifold mean, the L2-median, and the flag mean are alternative averages that can be substituted directly for the Karcher mean in many applications.  

This paper evaluates the characteristics and performance of these four averages on synthetic and real-world data. While the Karcher mean generalizes the Euclidean mean to the Grassman manifold, we show that the extrinsic manifold mean, the L2-median, and the flag mean behave more like medians and are therefore more robust to the presence of outliers among the subspaces being averaged. We also show that while the Karcher mean and L2-median are computed using iterative algorithms, the extrinsic manifold mean and flag mean can be found analytically and are thus orders of magnitude faster in practice. Finally, we show that the flag mean is a generalization of the extrinsic manifold mean that permits subspaces with different numbers of dimensions to be averaged. The result is a ”cookbook” that maps algorithm constraints and data properties to the most appropriate subspace mean for a given application.


# File Usage #

TODO



# Contact #

In case of questions, suggestions, problems etc. please send an email.

Tim Marrinan:  
marrinat@oregonstate.edu

This matlab package is also hosted at:  
http://www.tmarrinan.com/research-interests-code/  
https://www.cs.colostate.edu/~vision/summet/

