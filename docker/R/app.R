## Only run examples in interactive R sessions
###########rerooting application UI############UI= user interface######
##setwd("/home/rita/Desktop/UCT-WORK/ROOTING/")
#setwd("/srv/shiny-server/app/")
#setwd("~/Downloads/rita/docker/data")
library("phangorn")
library("ape")
library("devtools")
library("magrittr")
library("ggplot2")
library("Biostrings")
library("ggtree")
library(dplyr)
#library("ggtree")
library("shiny")
library("phytools")
library(shinydashboard)
library(caper)
library(geiger)
library("DT")
library(dashboardthemes)
if (interactive()) {
  
  
  #req({mydata})
  #tree<-mydata
  #reroot.all(tree)
  
  frow1<-fluidRow( valueBoxOutput("value1"), valueBoxOutput("value22"))
  frow2<-fluidRow( box( title = "Tree Maximum Likelihoods Distribution Plot",status = "primary", solidHeader=TRUE, collapsible=FALSE,plotOutput("treebylikelihood",height ="300px")),box(    title = "Rooted Phylogenetic Tree Plot"    ,status = "primary"   ,solidHeader = TRUE    ,collapsible = TRUE     ,plotOutput("besttree", height = "300px") ))
  
  header <- dashboardHeader(titleWidth = 1250,title = "Phylogenetic Analysis: Rooting Phylogenetic Trees using Non-Reversible Nucleotide Substitution Models")
  
  sidebar <- dashboardSidebar(menuItem("OWNERSHIP", tabName = "model",
                                        box(title="Developer.",background = "green",width = 150,
                                            "The web application has been developed by Rita Sianga (PhD student) and professor Darren Martin, University of Cape Town.")),
                                    
                                              
                                       menuItem("ABOUT",tabName="model",box(title="Application Information",background = "green",width = 150,
                                                    "The application roots phylogenetic trees using non-reversible nucleotide substitution stocastic models. AT the core of this rooting method is a twelve (12) rate parameter model known as NREV12 which is defined using a software package HyPhy (hypothesis testing using phylogenies). For any uploaded phylogenetic tree, the application uses R phytools function reroot that re-roots a phylogenetic tree at an arbitrary position along every edge. ")),
                                                    
                                       menuItem("INSTRUCTIONS",tabName = "model",box(title="Instructions",
                                                    background="green", width = 150,
                                                    "Upload a phylogenetic tree in newick format and a sequence alignment titled seq.fasta. Once the upload is complete, Run the application by clicking on the run button.")),menuItem("DATA INPUT", tabName = "model",
                                                                   fileInput("file1", "Choose newick file", accept = c(".nwk"), multiple = FALSE),fileInput("file2", "Choose a fasta file saved  seq.fasta", accept = c(".fasta"), multiple = FALSE),actionButton("Run",label = "Run"),
                                                                   checkboxInput("header", "Header", TRUE)),width=350,tabItem(tabName="import",mainPanel(
                                                                                                                                                                              tableOutput("results"),textOutput("system")
                                                                                                                                                                            )))
  
  body <- dashboardBody(downloadButton("downloadData", "Download the rooted tree"),shinyDashboardThemes(
    theme = "blue_gradient"
  ),frow1,frow2,tags$head(
    tags$style(HTML("
          .content-wrapper {
            background-color: white !important;
          }
          .main-sidebar { font-size: 20px; }
          .main-sidebar { font-color: green !important; }
          .main-sidebar {
            background-color: green !important;
          }
          .main-header .logo {
        font-family: Georgia, Times, Times New Roman, serif;
        font-weight: bold;
        font-size: 24px;
          }
        #.btn-file {  
         #    background-color:red; 
          #   border-color: red; 
           #  }

             .progress-bar {
             background-color: red;
             }
        "))
  ))
  
  ui<-dashboardPage(header, sidebar, body,skin = "black")
  
  
  
  source("rroot.R")
  server <- function(input, output,session) {
    
    read_file<- reactive({
      req(input$file1)
      infile<-input$file1
      return(infile$datapath)
      
      
      
      
    })
    read_newick_tree_and_run <- eventReactive(input$Run,{
      req(input$file1$datapath)
      tree_load<-read_file()
      tree<-read.tree(tree_load)
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
      output$treebylikelihood <- renderPlot({ ggplot(data=b, aes(x=b$V1))+geom_density(fill="#868686FF", color = "#EFC000FF")+
          geom_vline(aes(xintercept = mean(b$V1)),linetype = "dashed", size = 0.6)+geom_rug()+xlab("Tree Log-Likelihoods")+
          theme_classic()})
      return(reroot.all(tree))
      
      
      
      
    })
    
    
    output$results <- renderTable({read_newick_tree_and_run ()$datapath})
    #output$system <- renderText({system("bash REROOT.sh",wait = TRUE)})
    output$downloadData<-downloadHandler(
      filename = function(){paste("bestTree","nwk",sep = ".")}, 
      content = function(RootedTree){
        file.copy("bestTree.nwk", RootedTree)
      },
      contentType = "text/nwk"
    )
    #system("bash REROOT.sh",wait = TRUE)
    
    
    
    
    
    
  }  
  shinyApp(ui, server)
}
