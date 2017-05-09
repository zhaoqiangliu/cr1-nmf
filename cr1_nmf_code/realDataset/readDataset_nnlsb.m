% [ V,K_est,timelimit ] = readDataset_nnlsb( dataset_name )
%
% Input.
%   dataset_name   : a string; the name of the dataset
%
% Output.
%   V              : the data matrix obtained
%   K_est          : latent dimensionality
%   timelimit      : maximal running time

function [ V,K_est,timelimit ] = readDataset_nnlsb( dataset_name )

if(strcmp(dataset_name,'CK'))
    % can be downloaded from http://www.consortium.ri.cmu.edu/ckagree/
    % Database Description
    % Image Data:Cohn-kanade.tgz (1.7gb): zipped file of image data in png format 
    % -585 directories (subject + session)
    % -97 subject directories
    % -8795 image files
    X = load('../dataset/CK49times64.mat');V=X.V; 
    K_est = 97;
    timelimit = 1200;
elseif(strcmp(dataset_name,'faces94'))
    % can be downloaded from http://cswww.essex.ac.uk/mv/allfaces/faces94.html
    % Database Description
    % Number of  individuals: 152
    % Image resolution: 200 by 180 pixels (portrait format)
    % Contains images of male and female subjects in separate directories
    % 3040 images; 20 images per person
    X = load('../dataset/faces94.mat');V=X.V;
    K_est = 152;
    timelimit = 8000;
elseif(strcmp(dataset_name,'Georgia Tech'))
    % can be downloaded from http://www.anefian.com/research/face_reco.htm
    % 750 images; 15 images per person
    % original size: 480*640=307200; reshaped to size: 49*64=3136
    X = load('../dataset/Georgia_Tech.mat');V=X.V;
    K_est = 50;
    timelimit = 5000;
elseif(strcmp(dataset_name,'PaviaU'))
    % can be downloaded from http://www.ehu.eus/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes
    X = load('../dataset/PaviaU.mat'); V_temp = X.paviaU;
    V = reshape(V_temp,size(V_temp,1)*size(V_temp,2),size(V_temp,3));
    K_est = 9;
    timelimit = 300;
else
    fprintf('Error, no such dataset!\n');
    return;
end
