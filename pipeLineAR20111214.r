 ######################################################################

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




