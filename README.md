# NBD-GMMC

Finding network biomarkers from gene co-expression networks (GCNs) has attracted a lot of research interest. A network biomarker is a topological module, i.e., a group of densely connected nodes in a GCN, in which the gene expression values correlate with sample labels. Compared with biomarkers based on single genes, network biomarkers are not only more robust in separating samples from different categories, but also able to better interpret the molecular mechanism of the disease. The previous network biomarker detection methods either employ distance based clustering methods or search for cliques in a GCN to detect topological modules. The first strategy assumes that the topological modules should be spherical in shape, and the second strategy requires all nodes to be fully connected. However, the relations between genes are complex, as a result, genes in the same biological process may not be directly, strongly connected. Therefore, the shapes of those modules could be oval or long strips. Hence, the shapes of gene functional modules and gene disease modules may not meet the aforementioned constraints in the previous methods. Thus, previous methods may break up the genes belonging to the same biological process into different topological modules due to those constraints. To address this issue, we propose a novel network biomarker detection method by using Gaussian mixture model clustering which allows more freedom in the shapes of the topological modules. We have evaluated the performance of our method on a set of eight TCGA cancer datasets. The results show that our method can detect network modules that possess better discriminate power, and provide biological insights.

This is our implementation of network biomarker detection using Gaussian mixture model.

The input are TCGA RNA sequence file that can be accessed from the official website of TCGA. 

To use the please remove the second line of the TCGA files, the code will automatically generate sample labels based on the barcode of the samples.

If you have any problem in running our code, please do not heastite to contact me at hanzhang89@qq.com

