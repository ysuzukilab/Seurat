library(rmarkdown)

inputFile  <- commandArgs(trailingOnly = TRUE)[1] #input .Rmd file
outputFile <- commandArgs(trailingOnly = TRUE)[2] #output .html file

render(inputFile, output_format = "html_document", output_file = outputFile)








