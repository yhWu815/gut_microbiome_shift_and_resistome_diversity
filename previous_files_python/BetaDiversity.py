from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix
import unifrac
import os
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np

class CustomDistanceMatrix:
    def __init__(self, otu_table):
        self.otu_table = otu_table
        self.index = otu_table.index
    
    def EasyDistanceMatrix(self, metric):
        """
        Calculte Bray-Curtis, Euclidean, Canberra, Jaccard distance matrix.
        metric: "braycurtis", "canberra", "euclidean", "Jaccard"
        """
        #from scipy.spatial.distance import pdist, squareform
        #from skbio.stats.distance import DistanceMatrix
        if metric in ["braycurtis", "canberra", "euclidean", "jaccard"]:
            full_distance_matrix = squareform(pdist(self.otu_table, metric=metric))
            full_distance_matrix = DistanceMatrix(full_distance_matrix, ids=self.index)
            return full_distance_matrix
        else:
            print(f"Do not support {metric}。Please choose UnifracDistanceMatrix method.")
    
    def UnifracDistanceMatrix(self, metric, biom_format_data, phylogenic_tree, threads_num):
        """
        Calculate weighted_unifrac、unweighted_unifrac distance matrix。
        biom_format_data: biom_format otu table
        phylogenic_tree: otu_tre created by mafft and fasttree
        threads_num: the number of threads used to calculate the distance matrix
        """
        #import unifrac
        #import os
        #from skbio.stats.distance import DistanceMatrix
        if "unweighted" in metric:
            full_distance_matrix = unifrac.unweighted(table=biom_format_data, phylogeny=phylogenic_tree, threads=threads_num)
            full_distance_matrix = DistanceMatrix(full_distance_matrix, ids=self.index)
            return full_distance_matrix
        else:
            full_distance_matrix = unifrac.weighted_normalized(table=biom_format_data, phylogeny=phylogenic_tree, threads=threads_num)
            full_distance_matrix = DistanceMatrix(full_distance_matrix, ids=self.index)
            return full_distance_matrix

class CustomPCOA:
    def __init__(self, distance_matrix):
        self.distance_matrix = distance_matrix
    
    def PCOA(self):
        """
        reduce the dimensionality of the distance matrix using PCoA method.
        """
        #from skbio.stats.ordination import pcoa
        pcoa_results = pcoa(self.distance_matrix)
        return pcoa_results

    @staticmethod
    def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **Kwargs):
        """
        add confidence ellipse to the chart.
        parameters:
        x, y -- the coordinates of the data points
        ax -- the axis of the chart
        n_std -- the multiple of the standard deviation, which determines the size of the ellipse.
        n_std=3 which makes the ellipse enclose 98.9% of the points if the data is normally distributed
        facecolor -- the fill color of the ellipse
        **kwargs -- other parameters of the Ellipse function
        """
        #from matplotlib.patches import Ellipse
        #import matplotlib.transforms as transforms
        #import numpy as np
        if x.size != y.size:
            raise ValueError("x and y must be the same size")
        
        cov = np.cov(x, y)
        pearson = cov[0,1]/np.sqrt(cov[0,0] * cov[1,1])

        # Using a special case to obtain the eigenvalues of this
        # two-dimensionl dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse((0,0), width=ell_radius_x * 2, height=ell_radius_y*2,
                        facecolor=facecolor, **Kwargs)
        
        # Calculating the stdandard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0,0]) * n_std
        #scale_y = np.sqrt(cov[1,1] * n_std) there's a mitake
        #FIXME: the correct way to calculate the scale_y
        scale_y = np.sqrt(cov[1,1]) * n_std

        transf = transforms.Affine2D() \
            .rotate_deg(45) \
            .scale(scale_x, scale_y) \
            .translate(np.mean(x), np.mean(y))
        
        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse)