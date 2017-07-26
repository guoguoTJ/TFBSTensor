%function [Mat_return]=TFbsTensor(Matrix_o,index_M, test_Index,R,maxiters)
%   TFbsTensor imputes the missing TF-biding data based on the observed datasets.
%   It models the existing TF-binding datasets as a 3-mode tensor, where the three modes represent the TF, cell line, and genomic locus. Then it fits a weighted CP
%    model to the  tensor  with missing values via optimization. Finally, the
%   missing entries can be predicted.
%
%   Mat_return = TFbsTensor(Matrix_o,index_M, test_Index,R) fits an R-component weighted CANDECOMP/PARAFAC
%   (CP) model to the tensor that models the observed datasets Matrix_o,
%   and predicts the data for samples that indexed by test_Index. The result
%   Mat_return is the predicted data for the sample( or samples) indexed by test_Index.
%
%   Mat_return=TFbsTensor(Matrix_o,index_M, test_Index,R,maxiters)
%   The papameters and the details are as follows:
%      'Matrix_o' -  It is the matrix of the observed datasets. Each row represents the genome-wide TF-binding profile  of a  TF_cell
%       sample. The rows are the samples, the comlums are the positions in the genome for each sample.
%      'index_M' - index_M is the sample index for the matrix Matrix_o. Each row of index_M represents the {TF,cell} index of a sample. Specifically, index_M=[TF_index, cell_index];
%      'test_Index' - For the samples you want to predict, test_Index represents the index of the samples.  Specificlly, test_Index=[TF_index, cell_index];
%      'R' - R is the rank of the tensor. That is, the number of components in the  CP model
%      'maxiters' - Maximum number of iterations. The default value is  50.
%

%% Extract  dimensions of tensor and norm of obseved matrix.
normMat=norm(Matrix_o,'fro');% Frobenious norm of obseved matrix Matrix_o;

TF_index=index_M(:,1);
cell_index=index_M(:,2);
num_g=size(Matrix_o,2);
num_m=length(unique(TF_index));
num_c=length(unique(cell_index));
num_observe= length(TF_index);
TF_uni=unique(TF_index);% the TF index of the tensor.
cell_uni=unique(cell_index); % the cell index of the tensor .
% get the index of the observed data in the tensor.
O=zeros(size(index_M));
for i=1:size(index_M,1)
    O(i,1)=find(TF_uni== index_M(i,1));
    O(i,2)=find(cell_uni ==index_M(i,2));
end

% get the index of the data to predict in the tensor.
P_index=zeros(size(test_Index));
for i=1:size(test_Index,1)
    P_index(i,1)=find(TF_uni== test_Index(i,1));
    P_index(i,2)=find(cell_uni ==test_Index(i,2));
end

%% Initialization
X=zeros(num_observe,R);
digital = 1e-8;
fit=[];fitold=0;
w=1e-9;
%% CP_als, the main function
for iters=1:8
    %% CP_als, the main function
    %%Fit CP using cp_als by ignoring missing entries
    % Initialization for factor matrices.
    M= normrnd(0,digital,num_m,R);
    C = normrnd(0,digital,num_c,R);
    G= normrnd(0,digital,num_g,R);
    I=eye(R);

    % Main Loop of the ALS algorithm: Iterate until convergence
    for iter = 1:maxiters
        % Step1: Given factor matrix M and C, update G.
        for i=1:num_observe
            X(i,:)=(M(O(i,1),:).*C(O(i,2),:));
        end
        G=(((X')*X+w*I)\(X')*Matrix_o)';
        PP = G'*Matrix_o';
        
        %Step 2: Given fact matrix C and G, update M.
        A1 = ((G'*G));
        for i = 1:num_m
            o= O(O(:,1)==i,:);
            A=0; B=0;
            if sum(O(:,1)==i)==0
                continue;
            end
            
            P = PP(:,O(:,1)==i);
            len=size(o,1);
            for j = 1:len
                B = B + (C(o(j,2),:))'.*P(:,j);
                tmp = (C(o(j,2),:))';
                A = A+diag(tmp)*A1*diag(tmp);
            end
            M(i,:)=((A+w*I)\B)';
        end
        
        %Step 3:Given factor matrix M and G, update C.
        for i=1:num_c
            o=O(O(:,2)==i,:);%
            A=0;B=0;
            if sum(O(:,2)==i)==0
                continue;
            end
            P = PP(:,O(:,2)==i);
            len=size(o,1);
            for j = 1:len
                B = B + (M(o(j,1),:))'.*P(:,j);
                tmp = (M(o(j,1),:))';
                A = A + diag(tmp)*A1*diag(tmp);
            end
            C(i,:)=((A+w*I)\B)';
        end
 
    end

    
    % Get the predicted matrix of Matrix_o from the predicted tensor.
    TF_o_index=O(:,1);
    cell_o_index=O(:,2);
    len=length(TF_o_index);
    Mat_o_impute = zeros(num_g,len,'single');%
    MC=M(TF_o_index,:).*C(cell_o_index,:);
    GMC=gpuArray(single(MC));
    GG=gpuArray(single(G));
    for k = 1:len
        Mat_o_impute(:,k)=gather(sum((repmat(GMC(k,:),num_g,1)).*GG,2));
    end
    Mat_o_impute=Mat_o_impute';
    
    % Cumpute fitness.
    res(iter)=norm((Mat_o_impute-Matrix_o),'fro')/normMat; %fitniss of the predicted data on observed entries.
    
    % Get the predicted data of the test sample.
    TF_index_test=P_index(:,1);
    cell_index_test=P_index(:,2);
    len=length(TF_index_test);
    Mat_p = zeros(num_g,length(TF_index_test),'single');%
    MC=M(TF_index_test,:).*C(cell_index_test,:);
    GMC=gpuArray(single(MC));
    GG=gpuArray(single(G));
    for k = 1:len
        Mat_p(:,k)=gather(sum((repmat(GMC(k,:),num_g,1)).*GG,2));%
    end
    Mat_p=Mat_p';
    
    % choose the best Mat_p
    if iters==1
        Mat_return=Mat_p;
    else
        if  fit >fitold
            Mat_return=Mat_p;
        end
    end
    fitold=fit;
end



