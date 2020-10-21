# Bayesian Causal Structural Learning with Zero-Inflated Poisson Bayesian Networks

This repository contains the implementation of a zero-inflated Poisson Bayesian network (ZIPBN) proposed by "Bayesian Causal Structural Learning with Zero-Inflated Poisson Bayesian Networks" by Junsouk Choi, Robert Chapkin, and Yang Ni. Specifically, our code in this repository will reproduce the results corresponding to Table 2 in the paper. We hope that the provided code is helpful to give a complete description of the procedure we did. 

## Requirements

Our implemenation requires some dependencies. Please run the following codes to the dependencies:

``` r
pkgs <- c("igraph", "pscl", "glmnet", "MXM", "foreach", "doParallel", "doRNG")
sapply(pkgs, install.packages, character.only = TRUE)
```

## Training

To train the model(s) in the paper, run this command:

```train
python train.py --input-data <path_to_data> --alpha 10 --beta 20
```

>📋  Describe how to train the models, with example commands on how to train the models in your paper, including the full training procedure and appropriate hyperparameters.

## Evaluation

To evaluate my model on ImageNet, run:

```eval
python eval.py --model-file mymodel.pth --benchmark imagenet
```

>📋  Describe how to evaluate the trained models on benchmarks reported in the paper, give commands that produce the results (section below).

## Pre-trained Models

You can download pretrained models here:

- [My awesome model](https://drive.google.com/mymodel.pth) trained on ImageNet using parameters x,y,z. 

>📋  Give a link to where/how the pretrained models can be downloaded and how they were trained (if applicable).  Alternatively you can have an additional column in your results table with a link to the models.

## Results

Our model achieves the following performance on :

### [Image Classification on ImageNet](https://paperswithcode.com/sota/image-classification-on-imagenet)

| Model name         | Top 1 Accuracy  | Top 5 Accuracy |
| ------------------ |---------------- | -------------- |
| My awesome model   |     85%         |      95%       |

>📋  Include a table of results from your paper, and link back to the leaderboard for clarity and context. If your main result is a figure, include that figure and link to the command or notebook to reproduce it. 
