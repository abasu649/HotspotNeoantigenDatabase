library(shinythemes)
library(ggplot2)
library(dplyr)

# merged_data <- read.delim("./merged_file.tsv")
# df <- as.data.frame(merged_data)
# df <- df[complete.cases(df), ]

merged_data <- read.delim("./smaller_merged_file.tsv")
df <- as.data.frame(merged_data)
# df <- df[complete.cases(df), ]


tcga <- read.delim("./average_gene_expression_matrix_chr1-5.tsv")
rownames(tcga) <- tcga$X
tcga = subset(tcga, select = -X)

shinyApp(
  ui = tagList(
    navbarPage(
      theme = "paper", 
      "Hotspot Neoantigen Database",
      tabPanel("Table",
               sidebarPanel(
                 selectInput("gene", "Gene:", unique(merged_data$Gene.Name)),
                 sliderInput("ba", "Maximum Binding Affinity IC50 (nM):", 0, 500, 50),
                 sliderInput("freq", "Minimum HLA Frequency (%):", 0, 24, 12),
                 checkboxGroupInput("alleles", label = "HLA Alleles",
                                    choices = c("Class I" = 1, "Class II" = 2), selected = c(1, 2)),
                 actionButton("submit", "Submit", class = "btn-primary")
               ),
               mainPanel(
                 fluidRow(
                   column(
                     dataTableOutput(outputId = "table"), width = 12)
                 )
               )
      ),
      tabPanel("Plots",
               sidebarPanel(
                 selectInput("plotgene", "Gene:", unique(merged_data$Gene.Name)),
                 selectInput(inputId="plotvar", label="Variant:", choices=unique(merged_data$Variant.Type), multiple=TRUE, selected = unique(merged_data$Variant.Type)),
                 sliderInput("plotba", "Maximum Binding Affinity IC50 (nM):", 0, 500, 50),
                 sliderInput("plotfreq", "Minimum HLA Frequency (%):", 0, 24, 12),
                 actionButton("plotsubmit", "Submit", class = "btn-primary")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("User Defined Plot",
                            fluidRow(
                              column(width = 6,
                                     plotOutput("scatterplot1", height = 300)
                              ),
                              column(width = 6,
                                     plotOutput("scatterplot2", height = 300)
                              )
                            )
                   ),
                   tabPanel("Gene Expression", 
                            plotOutput("barplot"))
                 )
               )
      )
    )
  ),
  
  server = function(input, output) {
    
    reactiveBA <- eventReactive(input$submit, {
      input$ba
    })
    
    reactiveAlleles <- eventReactive(input$submit, {
      input$alleles
    })
    
    reactiveGene <- eventReactive(input$submit, {
      input$gene
    })
    
    reactivePlotGene <- eventReactive(input$plotsubmit, {
      input$plotgene
    })
    
    reactivePlotBA<- eventReactive(input$plotsubmit, {
      input$plotba
    })
    
    reactivePlotFreq <- eventReactive(input$plotsubmit, {
      input$plotfreq
    })
    
    reactiveFrequencies <- eventReactive(input$submit, {
      input$freq
    })
    
    reactivePlotVar <- eventReactive(input$plotsubmit, {
      input$plotvar
    })
    
    output$table <- renderDataTable({
      reactiveBA()
      reactiveAlleles()
      reactiveGene()
      reactiveFrequencies()
      plot_df <- df
      plot_df <- plot_df[plot_df$Gene.Name == isolate(input$gene),]
      plot_df <- plot_df[plot_df$Median.MT.Score <= isolate(input$ba),]
      plot_df <- plot_df[plot_df$HLA.Allele.Frequency >= isolate(input$freq),]
      
      if (length(isolate(input$alleles)) == 1) {
        if (isolate(input$alleles)[1] == 1) {
          plot_df <- plot_df[grepl("HLA-A", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-B", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-C", plot_df$HLA.Allele, fixed=TRUE),]
        } else {
          plot_df <- plot_df[!(grepl("HLA-A", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-B", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-C", plot_df$HLA.Allele, fixed=TRUE)),]
        }
      }
      colnames(plot_df) = c('Gene', 'Hotspot Mutation', 'Protein Change',
                       'Protein Position', 'Variant Type', 'MT Sequence',
                       'WT Epitope Sequence', 'HLA Allele', 'Global HLA Allele Frequency (%)',
                       'MT Median Binding Affinity IC50 (nM)', 'WT Median Binding Affinity IC50 (nM)',
                       'Median Fold Change', 'Peptide Length', 'Sub peptide Position', 'Mutation Position')
      plot_df
    }, options = list(scrollX = TRUE, searching = FALSE, pageLength = 10)
    )
    
    output$barplot <- renderPlot({
      reactivePlotGene()
      if (isolate(input$plotgene) %in% row.names(tcga)) {
        exp <- tcga[isolate(input$plotgene),]
        names(exp) <- NULL
        exp <- unlist(c(exp))
      } else {
        exp <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
      }
      
      
      cancer_name <- c("Adrenocortical carcinoma", "Bladder Urothelial Carcinoma", "Diffuse Large B-cell Lymphoma",
                       "Uterine Corpus Endometrial Carcinoma", "Skin Cutaneous Melanoma", "Head and Neck sqq cell carcinoma",
                       "Prostate adenocarcinoma", "Kidney renal papillary cell carcinoma", "Pancreatic adenocarcinoma",
                       "Sarcoma", "Cervical sqq and endocervical adenocarcinoma", "Colon adenocarcinoma",
                       "Lung sqq cell carcinoma", "Rectum adenocarcinoma", "Kidney renal clear cell carcinoma",
                       "Liver hepatocellular carcinoma", "Breast invasive carcinoma", "Ovarian serous cystadenocarcinoma",
                       "Uterine Carcinosarcoma", "Glioblastoma multiforme", "Kidney Chromophobe", "Thyroid carcinoma", 
                       "Brain Lower Grade Glioma", "Lung adenocarcinoma", "Mesothelioma", "Pheochromocytoma and Paraganglioma",
                       "Testicular Germ Cell Tumors", "Uveal Melanoma", "Thymoma", "Cholangiocarcinoma", "Esophageal carcinoma", 
                       "Stomach adenocarcinoma", "Acute Myeloid Leukemia")
      # Create data
      data <- data.frame(
        expression=exp,
        cancer_name = cancer_name
      )
      
      # Barplot
      ggplot(data, aes(x=reorder(cancer_name, expression), y=expression)) + 
        geom_bar(stat = "identity", fill=rgb(66, 139, 202, max = 255)) + 
        ggtitle(paste("Plot of TCGA PanCancer Atlas Expression Data for ", isolate(input$plotgene), sep = "")) +
        xlab("Cancer Type") + ylab("log2(expression + 1)") + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        coord_flip()
      
    }, height = 450)
    
    output$scatterplot1 <- renderPlot({
      reactivePlotGene()
      reactivePlotBA()
      reactivePlotFreq()
      reactivePlotVar()
      plot_df <- df
      plot_df <- plot_df[plot_df$Gene.Name == isolate(input$plotgene),]
      plot_df <- subset(plot_df, Variant.Type %in% as.list(isolate(input$plotvar)))
      plot_df <- plot_df[plot_df$Median.MT.Score <= isolate(input$plotba),]
      plot_df <- plot_df[plot_df$HLA.Allele.Frequency >= isolate(input$plotfreq),]
      plot_df <- plot_df[grepl("HLA-A", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-B", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-C", plot_df$HLA.Allele, fixed=TRUE),]
      ggplot(plot_df, aes(x=Median.MT.Score, y=HLA.Allele.Frequency)) +
        geom_point() + ggtitle(paste("Class I Neoantigens for ", isolate(input$plotgene), sep = "")) +
        xlab("Median MT Binding Affinity IC50 (nM)") + ylab("Global HLA Allele Frequency (%)") 
    })
    
    output$scatterplot2 <- renderPlot({
      reactivePlotGene()
      reactivePlotBA()
      reactivePlotVar()
      plot_df <- df
      plot_df <- plot_df[plot_df$Gene.Name == isolate(input$plotgene),]
      plot_df <- subset(plot_df, Variant.Type %in% as.list(isolate(input$plotvar)))
      plot_df <- plot_df[plot_df$Median.MT.Score <= isolate(input$plotba),]
      plot_df <- plot_df[plot_df$HLA.Allele.Frequency >= isolate(input$plotfreq),]
      plot_df <- plot_df[!(grepl("HLA-A", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-B", plot_df$HLA.Allele, fixed=TRUE) | grepl("HLA-C", plot_df$HLA.Allele, fixed=TRUE)),]
      ggplot(plot_df, aes(x=Median.MT.Score, y=HLA.Allele.Frequency)) +
        geom_point() + ggtitle(paste("Class II Neoantigens for ", isolate(input$plotgene), sep = "")) +
        xlab("Median MT Binding Affinity IC50 (nM)") + ylab("Global HLA Allele Frequency (%)") 
    })
    
    
  }

)