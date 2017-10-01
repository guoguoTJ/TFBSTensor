function [result] = evaluation(M,M_p,peak_position)
%% Metrics used to evaluate the predicted data.

%   parameters :
%     'M'  is the observed data for a sample and  'M_p' is the  predicted data for the corresponding sample.
%     'peak_position' -  a vector composed of {0,1} represents the potential peak regions in the genome, in which 1
%     represents  peak region  and 0 indicates that the position is not a peak region.

%The mian function for evaluation.
%% Metric Catch1obs and Catch1imp
 Catch1obs=overlap_caculation(M,M_p,0.01,0.05);
 Catch1imp= overlap_caculation(M_p,M,0.01,0.05);
 
 %% Metric AUPRC1 and AUPRC5
    [~,order] = sort(-M);
    label_1 = zeros(size(M));
    label_1(order(1:ceil(length(M)/100))) = 1;
    label_2 = zeros(size(M));
    label_2(order(1:ceil(length(M)/20))) = 1;
    
    AUPRC1 = AUPRC_calculation( M_p, label_1 );
    AUPRC5 = AUPRC_calculation( M_p, label_2 );

  %% Metric CorrinPeak
 [CorrinPeak] = CorrinPeak_calculation(M,M_p,peak_position);
 
 %%  Metric MAD
MAD= MAD_calculation( M, M_p );

result= [Catch1obs  Catch1imp   AUPRC1 AUPRC5 CorrinPeak MAD ];


function[fraction]= overlap_caculation(value1,value2,percent1,percent2)
% The function  for calculating  metrics Catch1obs and Catch1imp.

%  parameters:
%     'value1' and 'value2' represent the imputed and observed data, respectively. 
%     'percent1'  represents  that the top  'percent1' bins in  the imputed  data are labeled as positive.  'percent2' is for the  observed data.

%% the main function
numtop1 = ceil(length(value1)*percent1);
numtop2=ceil(length(value2)*percent2);

[~, order1] = sort(value1,'descend');
index1=order1(1:numtop1);
[~, order2] = sort(value2,'descend');
index2=order2(1:numtop2);
fraction = length(intersect(index1,index2))/numtop1;

function [AUPRC] = AUPRC_calculation( value, label )
%The function is for calculating metric AUPRC1 and AUPRC5.
%%
[~,~,~,AUPRC] = perfcurve(label, value, 1, 'xCrit', 'reca', 'yCrit', 'prec');

function [CorrinPeak] = CorrinPeak_calculation(M,M_p,peak_position)
% Metric CorrinPeak.

%   parameters for the method:
%     'M'- is the observed data for a sample and the 'M_p' is the  predicted data for the corresponding sample.
%     'peak_position' - indicates the peak regions for the sample, in which 1
%     represents  peak region in the position and 0 indicates that this position is not a peak  region.
%% the main function for calculate matric CorrinPeak
index = find(peak_position>0);
    M=M(:,index);
    M_p=M_p(:,index);
    CorrinPeak= corr(M',M_p','type','Spearman');
    
function [MAD] = MAD_calculation( value, value_p )
% Metric MAD
%%
MAD= mean(abs(value-value_p));

