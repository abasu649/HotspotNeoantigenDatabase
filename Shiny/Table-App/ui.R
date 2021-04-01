
practice <- read.delim("/Users/anamikabasu/Code/Database/Shiny/practice.tsv")

fluidPage(
  titlePanel("Hotspot Neoantigen Database"),

  # Create a new Row in the UI for selectInputs
  fluidRow(
    column(4,
        selectInput("allele",
                    "HLA Allele:",
                    c("All",
                      unique(as.character(practice$HLA.Allele))))
    ),
    column(4,
        selectInput("peptide",
                    "Peptide Length:",
                    c("All",
                      unique(as.character(practice$Peptide.Length))))
    )
  ),
  # Create a new row for the table.
  DT::dataTableOutput("table")
)
