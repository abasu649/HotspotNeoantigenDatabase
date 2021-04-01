
practice <- read.delim("/Users/anamikabasu/Code/Database/Shiny/practice.tsv")

function(input, output) {

  # Filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    data <- practice
    if (input$allele != "All") {
      data <- data[data$HLA.Allele == input$allele,]
    }
    if (input$peptide != "All") {
      data <- data[data$peptide == input$peptide,]
    }
    data
  }))

}

