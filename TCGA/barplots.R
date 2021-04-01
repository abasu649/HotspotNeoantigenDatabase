# Load ggplot2
library(ggplot2)

tcga <- read.delim("/Users/anamikabasu/Code/Database/TCGA/average_gene_expression_matrix_chr1-5.tsv")

casp <-tcga[1,-1]
names(casp)<- NULL
casp <- unlist(c(casp))

cancer_name <- c("Adrenocortical carcinoma", "Bladder Urothelial Carcinoma", "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                 "Uterine Corpus Endometrial Carcinoma", "Skin Cutaneous Melanoma", "Head and Neck squamous cell carcinoma",
                 "Prostate adenocarcinoma", "Kidney renal papillary cell carcinoma", "Pancreatic adenocarcinoma",
                 "Sarcoma", "Cervical sqq cell carcinoma and endocervical adenocarcinoma", "Colon adenocarcinoma",
                 "Lung squamous cell carcinoma", "Rectum adenocarcinoma", "Kidney renal clear cell carcinoma",
                 "Liver hepatocellular carcinoma", "Breast invasive carcinoma", "Ovarian serous cystadenocarcinoma",
                 "Uterine Carcinosarcoma", "Glioblastoma multiforme", "Kidney Chromophobe", "Thyroid carcinoma", 
                 "Brain Lower Grade Glioma", "Lung adenocarcinoma", "Mesothelioma", "Pheochromocytoma and Paraganglioma",
                 "Testicular Germ Cell Tumors", "Uveal Melanoma", "Thymoma", "Cholangiocarcinoma", "Esophageal carcinoma", 
                 "Stomach adenocarcinoma", "Acute Myeloid Leukemia")
# Create data
data <- data.frame(
  cancer_id=colnames(tcga)[-1],  
  expression=casp,
  cancer_name = cancer_name
)

# Barplot
ggplot(data, aes(x=reorder(cancer_name, -expression), y=expression)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Polt of TCGA PanCancer Atlas Expression Data") +
  xlab("Cancer Type") + ylab("log2(expression + 1)")

