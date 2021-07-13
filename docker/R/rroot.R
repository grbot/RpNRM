library("phangorn")
library("ape")
library("devtools")
library("magrittr")
library("ggplot2")
library(dplyr)
#library("ggtree")
library("shiny")
library("phytools")
library(shiny)
library(shinydashboard)
library(phytools)
library(caper)
library(geiger)
library("DT")

reroot.all<-function(tree,unroot=TRUE){
  # n<-length(tree$tip.label)
  #E<-2*n-3
  if(unroot) tree<-unroot(tree)
  if(is.null(tree$edge.length)){
    edge.lengths<-FALSE
    tree$edge.length<-rep(1,nrow(tree$edge))
  } else edge.lengths<-TRUE
  trees<-vector(mode="list",length=nrow(tree$edge))
  for(i in 1:nrow(tree$edge)){
    trees[[i]]<-reroot(tree,tree$edge[i,2],
                       0.5*tree$edge.length[i])
    if(!edge.lengths) trees[[i]]$edge.length<-NULL
    filename <- paste (i ,"0.nwk", sep ="")
    write.tree(trees[[i]],file = filename)
    
  }
  class(trees)<-"multiPhylo"
  return(trees)
}

