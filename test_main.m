
% Load the data used for test. 
load('./DataTest1.mat')
%      'Matrix_o' -  It is the matrix of the observed datasets. Each row represents the genome-wide TF-binding profile  of a  TF_cell
%       sample. The rows are the samples, the comlums are the positions in the genome for each sample.
%      'index_M' - index_M is the sample index for the matrix Matrix_o. Each row of index_M represents the {TF,cell}  of a sample. Specifically, index_M=[TF, cell];
%      'test_Index' - For the samples you want to predict, test_Index represents the  {TF,cell}  of the samples.  Specificlly, test_Index=[TF, cell];
%   Notice that in the 'index_M' and 'test_Index', both TF and cell  are labeled by interger such as 1,2,3,4.....

%      'R' - R is the rank of the tensor. That is, the number of components in the  CP model
%      'maxiters' - Maximum number of iterations. The default value is  50.

% The main function used for prediction.
 Mat_p=TFbsTensor(Matrix_o,index_M, test_Index,R,maxiters);
