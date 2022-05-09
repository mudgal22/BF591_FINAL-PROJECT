
library(shiny)
library(ggplot2)
library(tidyverse)
library(DT)
library(colourpicker)
library(dplyr)
#options(shiny.maxRequestSize=30*1024^2)


ui <- 
  fluidPage(
    title = 'FinalProject591-Urvy Mudgal',
    
    #Samples
    tabPanel('Samples',
             sidebarLayout(
               sidebarPanel(
                 p('Load Sample series Matrix'),
                 fileInput('file_sample_info', 'select files to upload',
                           accept = c(
                             'text/csv',
                             'text/comma-separated-values',
                             '.csv'),
                           placeholder = "sample_info.csv")),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", dataTableOutput(outputId = 'sample_summary')),
                   tabPanel("Data", dataTableOutput(outputId = 'sample_data')),
                   tabPanel("Plot", plotOutput(outputId = 'sample_bar'),
                            plotOutput(outputId = 'sample_histogram')))
               )
               
             )
    ),
    
    #Counts Matrix Tab 
    tabPanel('Counts',
             sidebarLayout(
               sidebarPanel(
                 p('Load normalized_Counts Matrix file here'),
                 fileInput('file_counts', 'select file to upload',
                           accept = c(
                             'text/csv',
                             'text/comma-separated-values',
                             '.csv'),
                           placeholder = "normalized.csv"),
                 
                 radioButtons('pca_x', "select first PC",
                              choices = c('1','2','3','4'),
                              selected = '1'),
                 radioButtons('pca_y', "select second PC",
                              choices = c('1','2','3','4'),
                              # choices = c(1,2,3,4),
                              selected = '2'),
                 
                 
                 sliderInput('slider_variance', "filter genes based on % variance",
                             value = 40,
                             min = 0,
                             max = 100),
                 sliderInput('slider_zero', "filter # of genes with zero variance",
                             value = 0,
                             min = 0,
                             max = 68),
                 submitButton(text = "Apply Changes", width = "100%")),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", dataTableOutput(outputId = 'counts_summary')),
                   tabPanel("PCA", plotOutput("counts_pca")))
               )
               
             )
    ),
    
    #Required - Differential Expression Section 
    tabPanel('Differential Expression',
             sidebarLayout(
               sidebarPanel(
                 p('Load differential expression results'),
                 fileInput('de_input', 'Choose file to upload',
                           accept = c(
                             'text/csv',
                             'text/comma-separated-values',
                             '.csv'),
                           placeholder = "deseq_results.csv"),
                 radioButtons('xaxis', "choose column for the x axis",
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue", "padj"),
                              selected = "log2FoldChange"), 
                 
                 radioButtons('yaxis', "choose column for the y axis",
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue", "padj"),
                              selected = "pvalue"), 
                 
                 colourInput('accent_col',"Select threshold color","#69b3a7"), 
                 
                 colourInput('base_col',"Select base color","#404080"),
                 
                 sliderInput('update_slider', "select the magnitude of p-adjusted coloring:",
                             value = -5,
                             min = -23, 
                             max = 0), 
                 submitButton(text = "Apply Changes", width = "100%")
               ),
               
               # Show the volcano plot
               mainPanel(tabsetPanel(
                 tabPanel("Plot", {
                   plotOutput("volcano")
                 }),
                 tabPanel("Table",
                          dataTableOutput("results_table"))
               ))
             )), 
    
    
  ) 
#define server logic
server <- function(input, output) {
  
  #Sample tab
  
  #reactive function to load data
  load_sample<- reactive({
    read_sample <- read.delim(file = input$file_sample_info$datapath, sep = ",", header = TRUE) %>%
      data.frame() %>%
      return()
  })
  
  
  bar_plot <-function(df1) {
    bar <- df1 %>%
      ggplot(df1, mapping = aes(x=PMI, y=RIN, fill=Diagnosis)) +
      geom_bar(stat="identity", alpha = 0.74) +
      scale_fill_manual(values=c("#69b3a2", "#404080")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90)) 
    
    return(bar)
  }
  
  
  histogram_plot <- function(df1){
    hist<- df1 %>%
      ggplot(aes(x=RIN, fill = RIN)) +
      geom_histogram(color ="#e9ecef", fill = "#69b3a2",  alpha=0.5, position = 'identity') +
      theme_minimal()+
      scale_fill_manual(values=c("#69b3a2", "#404080")) +
      labs(fill="")
    
    return(hist)
  }
  
  #summary tab - summarize type and value for each column
  
  Column_Name <- c('Diagnosis','PMI', 'Age of Death', 'RIN')
  Type <- c( 'factor', 'double', 'double', 'double' )
  Mean_Value <- c('Control, HD', '15.2', '65.8','7.71')
  
  summary_data <- data.frame(Column_Name, Type, Mean_Value)
  
  output$sample_summary<-renderDataTable({
    req(input$file_sample_info$datapath)
    
    summary_data
    
  })
  
  
  summary_table <- function(data){ 
    sample_data <- read_csv('cbind_sample_info.csv')
    #   sample_data$Diagnosis <- as.factor(sample_data$Diagnosis)
    col_names <- colnames(sample_data)
    new_data <- drop_na(sample_data)
    RIN_mean <- mean(new_data$RIN)
    death_mean <- round(mean(new_data$Age_of_Death),3)
    PMI_mean <- round(mean(new_data$PMI),3)
    condition <- c('Neurologically_Normal', 'Huntington Disease')
    
    summary <- data.frame(col_names)
    all_means <- c(RIN_mean,death_mean,PMI_mean, Diagnosis)
    type <- c('Numeric', 'Intger', 'Numeric','Character') 
    combine <- cbind(summary, all_means)
    sample_summary <- cbind(summary,type)
    names(sample_summary) <- c("Column Name", "Mean/ Distinct Type", "Type") 
    return(sample_summary)
    
  }
  
  
  #output sample info table
  output$sample_data<-renderDataTable({
    req(input$file_sample_info$datapath)
    
    load_sample()
  })
  
  
  output$sample_bar <- renderPlot({
    req(input$file_sample_info$datapath)
    
    bar_plot(df1 = load_sample())
  })
  
  
  output$sample_histogram <- renderPlot({
    req(input$file_sample_info$datapath)
    
    histogram_plot(df1 = load_sample()) 
  })
  
  
  
  #Counts tab
  
  #load counts data
  norm_counts<- reactive({
    read_counts <- read.delim(file = input$file_counts$datapath, sep = ",", header = TRUE) %>%
      # if (is.null(read_counts))
      #   return(NULL)
      return(read_counts)
  })
  
  
  
  filter_zero <- function(data, slider1){
    nonzeros <- data[rowSums(data > 0) >= slider1, ]
    return(nonzeros)
  }
  
  #filter counts table 
  draw_counts_table <- function(data, slider1, slider2) {
    d <- data[c(-1)]
    nonzeros <- d[rowSums(d > 0) >= slider1, ]
    variances <- apply(nonzeros, 1, var)
    table <- nonzeros[variances < slider2, ]  
    
    return(table)
  }
  
  
  plot_pca <- function(data, meta_info, title="", a, b) {
    pca <- prcomp(t(data))
    plot_data <- meta_info
    #change string to integer, bc radio buttons need the input to be a string 
    a <- strtoi(a)
    b<- strtoi(b)
    plot_data$PC1 <- pca$x[ , a]
    plot_data$PC2 <- pca$x[ , b]
    percent_variance <- pca$sdev^2 / sum( pca$sdev^2 )
    print(percent_variance)
    pca_plot <- ggplot(plot_data,aes(x=PC1, y=PC2, col=condition)) +
      geom_point() +
      xlab(paste0("PC", a ,":",round(percent_variance[a] * 100),"% variance")) +
      ylab(paste0("PC", b , ": ",round(percent_variance[b] * 100),"% variance")) +
      ggtitle(title)
    return(pca_plot)
  }
  
  
  
  #output results -- counts matrix summary table
  output$counts_summary<-renderDataTable({
    req(input$file_counts$datapath)
    
    
    draw_counts_table(norm_counts(),slider1 = input$slider_zero, slider2 = input$slider_variance) 
  }) 
  
  #pca plot 
  output$counts_pca <- renderPlot({
    req(input$file_counts$datapath)
    
    
    plot_pca(norm_counts(), 
             meta_data, 
             title='Normalized Counts PCA',
             input$pca_x,
             input$pca_y)
    
  })
  
  # Differential Expression 
  
  #reactive function -- load data
  load_data<- reactive({
    inFile <- read.delim(file = input$de_input$datapath, sep = ",", header = TRUE) %>%
      rename(ENSEMBL_ID = X) %>%
      return()
  })
  
  #function for generating volcano plot
  volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    volc_plot <- dataf %>%
      select(c(x_name,y_name)) %>%
      mutate(volcano_status = ifelse(.[,y_name] <10^slider, 'TRUE', 'FALSE')) %>%
      ggplot()+
      geom_point(aes(x= !!sym(x_name),
                     y= -log10(!!sym(y_name)),
                     color = factor(volcano_status, levels = c("TRUE", "FALSE")))) +
      scale_color_manual(values = c(color1,color2)) +
      theme_minimal()+
      theme(legend.position = "bottom") +
      labs(x = x_name,
           y = paste0('-log10(',y_name,')'),
           color = paste0(y_name, ' < 1 x 10^',slider))
    
    return(volc_plot)
  }
  
  
  #Function to filter data frame table to rows that are above the slider selection
  draw_table <- function(dataf, slider) {
    new_table <- dataf %>%
      arrange(pvalue) %>%
      filter(padj < 10^slider) %>%
      mutate(pvalue = formatC(.$pvalue, digits = 2, format = 'e'),
             padj = formatC(.$padj, digits =2, format = 'e'))
    
    return(new_table)
  }
  
  
  #output results
  output$results_table<-renderDataTable({
    (req(input$de_input$datapath))
    
    draw_table(dataf = load_data(),
               slider = input$update_slider)
  })
  
}  
shinyApp(ui,server)
