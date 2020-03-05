###########rerooting application UI############UI= user interface######
#setwd("/home/gerrit/gerrit.botha@uct.ac.za/machine/projects/rita/code")
setwd("/srv/shiny-server/app/")
library("phangorn")
library("ape")
library("devtools")
library("magrittr")
library("ggplot2")
library("ggtree")
library("shiny")
library("phytools")
library(shiny)
library(shinydashboard)

frow1<-fluidRow( valueBoxOutput("value1"), valueBoxOutput("value22"))
frow2<-fluidRow( box( title = "Tree mLikelihood distribution",status = "primary", solidHeader=TRUE, collapsible=TRUE, plotOutput("treebylikelihood",height ="300px")),box(    title = "Best Tree"    ,status = "primary"    ,solidHeader = TRUE     ,collapsible = TRUE     ,plotOutput("besttree", height = "300px") ))
ui<-dashboardPage(dashboardHeader(title = "Phylogenetic Analysis"),dashboardSidebar(sidebarMenu(menuItem("Dashboard",tabName = "dashboard",icon("dashboard")))),dashboardBody(frow1,frow2) )

server<-function(input, output){
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
    trees
  }


  #input<-readline(prompt = "Enter tree in nwk")


  tree<-read.tree("HIVM.fasta.nwk")
  reroot.all(tree)
  system("bash REROOT.sh",wait = TRUE)

  b<-read.csv("data.csv",header = FALSE)
  for (i in 1:length(b)) {
    b<-read.csv("data.csv",header = FALSE)
    b1<-read.tree("data1.nwk")
    which.max(b$V1)
    write.tree(b1[i],file = "bestTree.nwk")
    T<-read.tree("bestTree.nwk")
    output$besttree <- renderPlot({plot(T) })


  }
  output$value1<-renderValueBox( {
    valueBox(formatC((2*length(tree$tip.label))-3 , format = "d",big.mark = ','), paste('Number of Trees'),icon = icon("stats",lib='glyphicon')      ,color = "purple")
  })


  b<-read.csv("data.csv",header = FALSE)
  b<-data.frame(b)
  output$treebylikelihood <- renderPlot({ ggplot(data=b, aes(x=b$V1))+geom_density(fill="#0073C2FF", color = "#0073C2FF")+
      geom_vline(aes(xintercept = mean(b$V1)),linetype = "dashed", size = 0.6)+geom_rug()+xlab("Tree Log-Likelihoods")+
      theme_classic()})

 # output$besttree <- renderPlot({tree<-read.tree("bestTree.nwk"), plotTree(tree) })




}

#barplot(b$V1)
shinyApp(ui,server)
