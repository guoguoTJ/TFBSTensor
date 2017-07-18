README for TFbsTensor
Introduction
The identification of DNA-binding sites for transcription factors (TFs) is vital for deciphering transcriptional regulation. Chromatin immunoprecipitation followed by high-throughput sequencing (ChIP-seq) enables the genome-wide identification of in vivo TF binding sites; however, it is not possible to map every TF in every cell line because of limitations such as the availability of sample material and experimental cost. This restriction limits the discovery power in gene regulation analysis. Here, we propose a tensor decomposition-based approach (named TFbsTensor) for imputation of transcription factor binding sites that works by exploiting the multiple correlations of TF binding signals from ChIP-seq datasets.
To predict TF binding profiles, TFbsTensor models the existing TF binding datasets as a 3-mode tensor, where the three modes represent the TF, cell line, and genomic locus. Then tensor decomposition is used to effectively capture the structural properties of the TF binding tensor by modeling multi-linear interactions among three latent factors, and the imputation of TF binding profiles is formulated as a tensor completion problem. 
TFbsTensor is implemented in MATLAB 2014b. GPU-accelerated computing is used in the implementation, thus a hardware supporting GPU computing is needed to use our code. Our algorithm can efficiently handle the large dataset.

Environment and Dependencies to use TFbsTensor:
•	MATLAB.
•	A hardware supporting GPU computing is needed. 
Usage of TFbsTensor
•	The file ‘TFbsTensor.m’ is the implemented method to predict TF-binding data. 
•	Details on using the method and the parameters are given in the comments in the code.  
Examples:
•	Input datasets: DataTest1.mat. 
•	Running test_main.m to test the TFbsTensor method. 
If having any questions about the TFbsTensor method, please feel free to contact guoweili_henu@126.com. 


