
library(shiny)
library(shinythemes)
library(markdown)
library(shinyWidgets)
library(DT)
library(data.table)
library(ggplot2)
library(gridExtra)
setwd("C:/Users/maja/Documents/proba/")
Sys.setlocale(category = "LC_ALL", locale = "Croatian")
allcells <- fread("Sampleinfo.txt")
allcells[,scaledvalue:=value/mean(value),.(celltype,variable)]

alltssfiles <- readRDS("allTss_smaller.RDS")
resultstfbs <- fread("results_TFBS.csv")

plotitTFBS1 <- function(x){
    b <- ggplot(x, aes(V1,V4, group=sampleName, color=sampleGroup))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=sampleGroup), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    b    
}
plotitTFBS <- function(x){
    a <- ggplot(x, aes(x,V4, group=sampleName, color=sampleGroup))+geom_line()+
        ggtitle(unique(x$transcriptionFactor))+theme_light()+ xlab("position ordered by value")+
        ylab("signal value")+geom_smooth( aes(group=sampleGroup), color="black", lty=2)+
        theme(legend.position = "none")
    
    b <- ggplot(x, aes(V1,V4, group=sampleName, color=sampleGroup))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=sampleGroup), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    print(grid.arrange(a,b))
    
}
sebytfbs <- fread("cfDNA_coverage_scaledper10kbandmean_allTFs_PER_SEsample.txt")
sebytfbs[order(value),orderedPosition:=1:.N,.(sesample,sample)]
kw_se_tf_bysample <- fread("kw_se_tf_bysample.txt")
sebysample <- readRDS("cfDNA_coverage_scaledper10kbandmean_PER_SEsample.RDS")
sebysample[order(scaledmeanascaledval),orderedPosition:=1:.N,.(sesample,sample)]
kw_se_bysample <- fread("kw_se_bysample.txt")

plotitSEbyregions <- function(x){
    a <- ggplot(x, aes(orderedPosition,scaledmeanascaledval, group=sample, color=condition))+geom_line()+
        ggtitle(paste(unlist(unique(x[,.(sesample,tissue_type,tissue,sample_origin)])),collapse = ";"))+
        theme_light()+ xlab("position ordered by value")+
        ylab("signal value")+geom_smooth( aes(group=condition), color="black", lty=2)+
        theme(legend.position = "none")
    
    b <- ggplot(x, aes(position,scaledmeanascaledval, group=sample, color=condition))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=condition), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    print(grid.arrange(a,b))
    
}

plotitSEbyTFBSregions <- function(x){
    a <- ggplot(x, aes(orderedPosition,value, group=sample, color=condition))+geom_line()+
        ggtitle(paste(unlist(unique(x[,.(sesample,tissue_type,tissue,sample_origin)])),collapse = ";"))+
        theme_light()+ xlab("position ordered by value")+
        ylab("signal value")+geom_smooth( aes(group=condition), color="black", lty=2)+
        theme(legend.position = "none")
    
    b <- ggplot(x, aes(position,value, group=sample, color=condition))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=condition), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    print(grid.arrange(a,b))
    
}

plotit <- function(x){
    a <- ggplot(x, aes(orderedPosition,value, group=variable, color=sampleGroup))+geom_line()+
        ggtitle(unique(paste(x$cell,x$identifier,sep="-")))+theme_light()+ xlab("position ordered by value")+
        ylab("signal value")+geom_smooth( aes(group=sampleGroup), color="black", lty=2)+
        theme(legend.position = "none")
    
    b <- ggplot(x, aes(position,value, group=variable, color=sampleGroup))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=sampleGroup), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    print(grid.arrange(a,b))
    
}

plotitscaled <- function(x){
    a <- ggplot(x, aes(orderedPosition,scaledvalue, group=variable, color=sampleGroup))+geom_line()+
        ggtitle(unique(paste(x$cell,x$identifier,sep="-")))+theme_light()+ xlab("position ordered by value")+
        ylab("signal value")+geom_smooth( aes(group=sampleGroup), color="black", lty=2)+
        theme(legend.position = "none")
    
    b <- ggplot(x, aes(position,scaledvalue, group=variable, color=sampleGroup))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=sampleGroup), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    print(grid.arrange(a,b))
    
}
plotitscaledat0 <- function(x){
    a <- ggplot(x, aes(orderedPositionat0,scaledvalueat0, group=variable, color=sampleGroup))+geom_line()+
        ggtitle(unique(paste(x$cell,x$identifier,sep="-")))+theme_light()+ xlab("position ordered by value")+
        ylab("signal value")+geom_smooth( aes(group=sampleGroup), color="black", lty=2)+
        theme(legend.position = "none")
    
    b <- ggplot(x, aes(position,scaledvalueat0, group=variable, color=sampleGroup))+geom_line()+theme_light()+
        xlab("position relative to TFBS")+
        ylab("signal value")+
        geom_smooth( aes(group=sampleGroup), color="black", lty=2, span=0.02, method="loess")+
        theme(legend.position = "bottom")
    
    
    print(grid.arrange(a,b))
    
}
# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = shinytheme("united"),
                   
    # Application title
    titlePanel("cfDNA analysis"),

    navbarPage("approaches tried",
               navbarMenu("TFBS by TF",
                          tabPanel("Introduction",
                                   tabsetPanel(
                                       tabPanel("cfDNA",
                                                mainPanel(
                                                    img(src = "cfDNA_UCSC.png", width=900, height=1000)),
                                       ),
                                       tabPanel("Profile",
                                            sidebarPanel(
                                                selectInput("transcriptionFactor", "Select transcription factor:", 
                                                   choices=sort(unique(alltssfiles$transcriptionFactor)), selected = "CTCF"),
                                                helpText("Ordered: CTCF, CTCFL, KLF
                                                Not ordered: CLOCK, ZGPAT, MYST")
                                       
                                            ),
                                   
                                   
                                            mainPanel(
                                                h3("cfDNA coverage around transcription factor binding sites"),
                                                plotOutput("TFBSplot")  
                                                ) 
                                            )
                                       )
                                   ),
                          
                          tabPanel("Comparison",
                                   sidebarPanel(
                                       selectInput("transcriptionFactorc", "Select transcription factor:", 
                                                   choices=sort(unique(alltssfiles$transcriptionFactor)), selected = "CTCF"),
                                       helpText("Ordered: CTCF, CTCFL, KLF
                                                Not ordered: CLOCK, ZGPAT, MYST")
                                       
                                   ),
                                   
                                   # Create a spot for the barplot
                          mainPanel(
                                       h3("cfDNA coverage around transcription factor binding sites"),
                                       
                                       plotOutput("TFBSplot2")  
                                   )    
                          ),
                          tabPanel("Results",
                                   # Create a spot for the barplot
                                   mainPanel(
                                       h3("Results"),
                                       DT::dataTableOutput("mytable")  
                                   )    
                          )
               ),
               navbarMenu("Superenhancers",
                            tabPanel("Introduction",
                              mainPanel(
                                       tabsetPanel(
                                            tabPanel("problem",
                                               h3("cfDNA coverage on superenhancer regions"),
                                               img(src = "cfDNA_SE_1.png", width=900, height=500)),
                                            tabPanel("solution1",
                                                     h3("cfDNA coverage on TFBS overlapping superenhancers"),                    
                                                     img(src = "cfDNA_SE_3.png", width=900, height=500)
                                                     ),
                                            tabPanel("solution2",  
                                                 h3("cfDNA coverage on superenhancer regions"),                    
                                                 img(src = "cfDNA_SE_2.png", width=900, height=500))
                                       )
                            )
                          ),
                          tabPanel("TF SE overlaps",
                           tabsetPanel(
                                      
                            tabPanel("profiles",
                                sidebarPanel(
                                       selectInput("sesampletype1", "Select tissue type :", 
                                                   choices=unique(sebytfbs$tissue_type)),
                                       selectInput("sesampletissue1", "Select tissue:", 
                                                   choices="", selected = ""),
                                       selectInput("sesample1", "Sample:", 
                                                   choices="", selected = ""),
                                       hr(),
                                       helpText("Cell types available")
                                   ),
                                   
                                   # Create a spot for the barplot
                                   mainPanel(
                                       h3("Superenhancers that overlap SE per sample"),
                                       plotOutput("sebytfbsregions")  
                                   )    
                                ),
                            tabPanel("results",
                                     # Create a spot for the barplot
                                     mainPanel(
                                         h3("Superenhancers that overlap SE per sample"),
                                         DT::dataTableOutput("kw_se_tf_bysample")  
                                         
                                     )    
                            )
                          )
                         ),
                         tabPanel("SE regions",
                                  tabsetPanel(
                                      
                                      tabPanel("profiles",
                                               sidebarPanel(
                                                   selectInput("sesampletype2", "Select tissue type :", 
                                                               choices=unique(sebysample$tissue_type)),
                                                   selectInput("sesampletissue2", "Select tissue:", 
                                                               choices="", selected = ""),
                                                   selectInput("sesample2", "Sample:", 
                                                               choices="", selected = ""),
                                                   hr(),
                                                   helpText("Cell types available")
                                               ),
                                               
                                               # Create a spot for the barplot
                                               mainPanel(
                                                   h3("cfDNA profile in SE regions per sample"),
                                                   plotOutput("sebysample")  
                                               )    
                                      ),
                                      tabPanel("results",
                                               # Create a spot for the barplot
                                               mainPanel(
                                                   h3("Superenhancers that overlap SE per sample"),
                                                   DT::dataTableOutput("kw_se_bysample")  
                                                   
                                               )    
                                      )
                                  )
                         )

                          ),
                
                tabPanel("Work in progress",
                         titlePanel("DNAse hypersensitivity sites"),
                         sidebarPanel(
                             selectInput("description2", "Select type:", 
                                         choices=unique(allcells$Description)),
                             selectInput("celltype2", "Sample:", 
                                         choices="", selected = ""),
                             hr(),
                             helpText("Cell types available")
                         ),
                         mainPanel(
                             tabPanel("scaled",
                                h2("cfDNA around TFBS that overlap DNAse hypersensitivity sites"),
                                span("Recently published "),a("paper from ENCODE ",href="https://www.nature.com/articles/s41586-020-2528-x?proof=t"), 
                                span(" describes DNAse hypersensitivity sites in 243 human cell and tissue types and states."),
                                hr(),
                                p(" The plot shows cfDNA average coverage over positions relative to TFBS.
                                The cfDNA coverage is scaled for total library size and divided by mean for each sample. 
                                "),
                                plotOutput("scaledplot")
                         )
                )
                
                )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
    
    
    output$TFBSplot <- renderPlot({
        x <- alltssfiles[transcriptionFactor==input$transcriptionFactor]
        plotitTFBS1(x)
        
    }, height = 400, width = 600)
    
    output$TFBSplot2 <- renderPlot({
        x <- alltssfiles[transcriptionFactor==input$transcriptionFactorc]
        plotitTFBS(x)
        
    }, height = 600, width = 600)
    output$mytable = DT::renderDataTable({
        resultstfbs
    })
    output$kw_se_tf_bysample = DT::renderDataTable({
        kw_se_tf_bysample
    })
    output$kw_se_bysample = DT::renderDataTable({
        kw_se_bysample
    })
    
    output$phonePlot <- renderPlot({
        x <- allcells[celltype==input$celltype]
        plotit(x)
            }, height = 600, width = 600)
    ###########################################################################
    # SE
    output$sebytfbsregions <- renderPlot({
        x <- sebytfbs[sesample==input$sesample1]
        plotitSEbyTFBSregions(x)
        
    }, height = 600, width = 600)
    output$sebysample <- renderPlot({
        x <- sebysample[sesample==input$sesample2]
        plotitSEbyregions(x)
        
    }, height = 600, width = 600)
    ###########################################################################
    output$scaledplot <- renderPlot({
        x <- allcells[celltype==input$celltype2]
        plotitscaled(x)
        
    }, height = 600, width = 600)
    output$scaledat0plot <- renderPlot({
        x <- allcells[celltype==input$celltype3]
        plotitscaledat0(x)
        
    }, height = 600, width = 600)
    ############################################################################
    
    observeEvent(
        input$sesampletype1,
        updateSelectInput(session, "sesampletissue1", "Tissue:", 
                          choices = sort(unique(sebytfbs$tissue[sebytfbs$tissue_type==input$sesampletype1]))))
    observeEvent(
        input$sesampletissue1,
        updateSelectInput(session, "sesample1", "Sample:", 
                          choices = sort(unique(sebytfbs$sesample[sebytfbs$tissue==input$sesampletissue1]))))
    ##################
    observeEvent(
        input$sesampletype2,
        updateSelectInput(session, "sesampletissue2", "Tissue:", 
                          choices = sort(unique(sebysample$tissue[sebysample$tissue_type==input$sesampletype2]))))
    observeEvent(
        input$sesampletissue2,
        updateSelectInput(session, "sesample2", "Sample:", 
                          choices = sort(unique(sebysample$sesample[sebysample$tissue==input$sesampletissue2]))))
    observeEvent(
        input$description2,
        updateSelectInput(session, "celltype2", "Celltype", 
                          choices = sort(unique(allcells$celltype[allcells$Description==input$description2]))))

    ############################################################################
}

# Run the application 
shinyApp( ui = ui, server = server)
