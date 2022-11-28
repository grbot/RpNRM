
   

##rpnrm APPLICATION
#change to right directory
setwd("/srv/shiny-server/rpnrm")
options(shiny.sanitize.errors = FALSE)
library("phangorn")
library("ape")
library("devtools")
library("magrittr")
library("ggplot2")
library("Biostrings")
library("ggtree")
library("dplyr")
#library("ggtree")
library("shiny")
library("phytools")
library("shinydashboard")
library("caper")
library("geiger")
library("DT")
library("dashboardthemes")




frow1<-fluidRow( valueBoxOutput("value1"), valueBoxOutput("value22"))
frow2<-fluidRow( box( title = "Tree Maximum Likelihoods Distribution Plot",status = "primary", solidHeader=TRUE, collapsible=FALSE,plotOutput("treebylikelihood",height ="300px")),box(    title = "Rooted Phylogenetic Tree Plot"    ,status = "primary"   ,solidHeader = TRUE    ,collapsible = TRUE     ,plotOutput("besttree", height = "300px") ),box(    title = "Model Test Results"    ,status = "primary"   ,solidHeader = TRUE    ,collapsible = TRUE     ,textOutput("test") ))
frow3<-fluidRow(valueBoxOutput("aic"))
header <- dashboardHeader(titleWidth = 1250,title = "Phylogenetic Analysis: Rooting Phylogenetic Trees using a Non-Reversible Nucleotide Substitution Model")

sidebar <- dashboardSidebar(menuItem("OWNERSHIP", tabName = "model",
                                     box(title="Developer.",background = "green",width = 150,
                                         "The web application has been developed by Rita Sianga-Mete (PhD student), University of Cape Town.")),
                            
                            
                            menuItem("ABOUT",tabName="model",box(title="Application Information",background = "green",width = 150,
                                                                 "RpNRM is a web application that roots phylogenetic trees using a 12-rate non-reversible nucleotide substitution stocastic model (NREV12). For any uploaded phylogenetic tree, the application will first run a model test using hyphy (hypothesis testing using phylogenies) to compare AIC scores for GTR (6-rate reversible model), NREV6 (a six-rate non-reversible model and NREV12 model. If NREV12 is the best suited model for the dataset, RpNRM will then root the phylogenetic tree on the most porbable branch. If NREV12 is not the best suited model for your dataset, we advise you root the tree using other softwares that make use of the reversible models.")),
                            
                            menuItem("INSTRUCTIONS",tabName = "model",box(title="Instructions",
                                                                          background="green", width = 150,
                                                                          "Upload a phylogenetic tree in newick format and a sequence alignment in fasta format. Once the upload is complete, Run the application by clicking on the run button.")),menuItem("DATA INPUT", tabName = "model",
                                                                                                                                                                                                                                                             fileInput("file1", "Choose newick file", accept = c(".nwk"), multiple = FALSE),fileInput("file2", "Choose a fasta file", accept = c(".fasta"), multiple = FALSE),actionButton("Run",label = "Run"),
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
#source("~/rpnrmcorrectfile/rpnrm/docker/R/rroot.R")
server <- function(input, output,session) {
  
  read_file<- reactive({
    req(input$file1)
    infile<-input$file1
    return(infile$datapath)
  })
  read_file2<- reactive({
    req(input$file2)
    infile2<-input$file2
    return(infile2$datapath)
    
    
    
  })
  read_newick_tree_and_run <- eventReactive(input$Run,{
    req(input$file1$datapath)
    file.copy(input$file1$datapath, "/srv/shiny-server/rpnrm")
    tree_load<-read_file()
    tree<-read.tree(tree_load)
    reroot.all(tree)
    req(input$file2$datapath)
    ##copy from r datapath to direc
    file.copy(input$file2$datapath, "/srv/shiny-server/rpnrm")
    #oldname<-input$file2$name      #input$file2$datapath
    #seq <- file.path(dirname(input$file2$datapath), basename(input$file2$name))
    #file.rename(from=oldname, to=seq.fasta)
    #input$file2$datapath<-seq.fasta
    #file.copy(input$file2$datapath, "~/rpnrmcorrectfile/rpnrm/docker/data")
    system("bash /srv/shiny-server/rpnrm/aic.sh",wait = TRUE)
    v<-read.csv("/srv/shiny-server/rpnrm/aic.csv", header = FALSE, sep=";")
    if(v[3,] > v[1,] || v[3,] > v[2,]){
      #output$aic<-renderText({"MODEL TEST RESULTS: NREV12 is not the best fitting model for your dataset. RpNRM will not root the input file."})
      output$test<-renderText({"MODEL TEST RESULTS: NREV12 is not the best fitting model for your dataset. RpNRM will not root the input file."})
      stop()} else{
        #output$aic<-renderText({"MODEL TEST RESULTS: NREV12 is the best fitting model for your dataset, please wait while your input tree gets rooted."})
        output$test<-renderText({"MODEL TEST RESULTS: NREV12 is the best fitting model for your dataset, please wait while your input tree gets rooted."})
        #renderText("NREV12 is the best fitting model for your dataset, wait while your phylogenetic tree gets rooted.")
        
        
      }   
    system("rm /srv/shiny-server/rpnrm/0.nwk", wait = TRUE)  
    system("bash /srv/shiny-server/rpnrm/REROOT.sh",wait = TRUE)
    file.copy(input$file1$datapath, "/srv/shiny-server/rpnrm")
    
    b<-read.csv("/srv/shiny-server/rpnrm/data.csv",header = FALSE)
    for (i in 1:length(b)) {
      b<-read.csv("/srv/shiny-server/rpnrm/data.csv",header = FALSE)
      #b1<-read.tree("/srv/shiny-server/rpnrm/data1.nwk")
      LRT<-which.max(b$V1)
      m<-read.tree("0.nwk")
      rooted<-reroot(m,m$edge[LRT,2],
                         0.5*m$edge.length[LRT])
      #if(!edge.lengths) trees[[i]]$edge.length<-NULL
      write.tree(rooted,file = "rooted.nwk")
      T<-read.tree("rooted.nwk")
      output$besttree <- renderPlot({plot(T,no.margin=TRUE,edge.width=2) })
      
      
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
    filename = function(){paste("rooted","nwk",sep = ".")}, 
    content = function(RootedTree){
      file.copy("/srv/shiny-server/rpnrm/rooted.nwk", RootedTree)
    },
    contentType = "text/nwk"
  )
  #system("bash REROOT.sh",wait = TRUE)
  
  
  
  
  
  
}  
shinyApp(ui, server)

