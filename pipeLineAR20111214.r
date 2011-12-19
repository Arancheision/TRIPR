 ######################################################################
#This is the pipeline to get the mapping information of the transposase 
#cassette integrations
#The input is a txt file containg the first 5 columns of the bam file,
#The output is a list containing:
                               #integration points
                               #barcode
                               #color
######################################################################

source('TRIP/function_getTheBarcodes')

list.promoterX= get.the.barcodes('filename')

#to add the color we do the following:

color.domains= read.delim('GSE22069_Drosophila_chromatin_domains.txt', as.is=TRUE)
#To have the same names as in our list:
color.domains$seqname= sub('chr','', color.domains$seqname)

source('TRIP/function_getTheColor')
list.promoterX$color= get.the.color(list.promoterX)

#I'll have 7 lists, one per promoter


#In order to compare the mapping data with the RNAseq I'll order the
#iPCR data by barcode and color in only one list:
      
iPCR= list(barcode=list(HSP=unlist(list.HSP$barcode), p14= unlist(list.p14$barcode), 
  p20= unlist(list.p20$barcode), p18= unlist(list.p18$barcode), p10= unlist(list.p10$barcode), 
  p15= unlist(list.p15$barcode), p7= unlist(list.p7$barcode)), color= list( HSP=unlist(list.HSP$color),
  p14= unlist(list.p14$color), p20= unlist(list.p20$color), p18= unlist(list.p18$color), 
  p10= unlist(list.p10$color), p15= unlist(list.p15$color), p7= unlist(list.p7$color)))
