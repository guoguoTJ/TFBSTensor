 %% This is an  example to illustrate how to use our developed method to predict TF binidng data. 
 % Here we provide an example of predicting the TF binding data for an observed sample based on the developed method. And we also provide  evaluation metrics to valid the predicted data. 

clear;
load('./Data_example.mat');
% Parameters:
%   'Matrix_o' -  It is the matrix of the observed datasets. Each row represents a genome-wide TF-binding profile  of a  TF_cell
%    sample. The rows are the samples, the comlums are  genome positions.
%   'index_M' - index_M is the sample index for the matrix Matrix_o. Each row of index_M represents the {TF,cell}  of a sample. Specifically, index_M(i,:)=[TF, cell];
%   'test_Index' - For the samples you want to predict, test_Index represents the  {TF,cell}  of the samples.  Specificlly, test_Index=[TF, cell];
%    Notice that in  'index_M' and 'test_Index', both TF and cell  are labeled by interger such as 1,2,3,4.....
%    'value_test' - The real value of the  sample to be predicted.
%   'R' - R is the rank of the tensor. 
%   'maxiters' - Maximum number of iterations. The default value is  50.

% The main function  for prediction.
 value_p=TFbsTensor(Matrix_o,index_M, test_Index,R,maxiters);
 
 %% evaluation
 % Six metrics for evaluation. Metrics= [Catch1obs, Catch1imp, AUPRC1, AUPRC5, CorrinPeak, MAD ];
 % Parameters:
 %   'peak' - it represents all  purtative peak regions in the genome,  in which 1
 %    represents  peak region in the position and 0 indicates that this position is not a peak  region.It is used for the calculation of metric CorrinPeak.
  
 result_eval=evaluation(value_test,value_p,peak);
  
