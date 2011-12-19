      
      
#load the data:
      
      
gDNA= read.delim('1610-CCG-count-min5.txt', sep='\t', as.is= TRUE, fill=TRUE, flush=TRUE)
names (gDNA) = c('barcode', 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')
       
RNA= read.delim('1610-TTA-count-min5.txt', sep='\t', as.is= TRUE, fill=TRUE, flush=TRUE)
names (RNA) = c('barcode', 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')
      
#order the iPCR data by barcode:
      
iPCR= list(barcode= list(HSP=unlist(list.HSP$barcode), p14= unlist(list.p14$barcode), 
  p20= unlist(list.p20$barcode), p18= unlist(list.p18$barcode), 
  p10= unlist(list.p10$barcode), p15= unlist(list.p15$barcode), 
  p7= unlist(list.p7$barcode)), color= list( HSP=unlist(list.HSP$color),
  p14= unlist(list.p14$color), p20= unlist(list.p20$color), p18= unlist(list.p18$color), 
  p10= unlist(list.p10$color), p15= unlist(list.p15$color), p7= unlist(list.p7$color)))
      
##########################################################################################
#in order to compare the barcodes in the 3 data sets (iPCR, gDNA, RNA) I make a list 
#with all of them, by promoter.
 
barcode.list= list()
promoters= c( 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')

for (i in promoters) {
    
    barcode.list[[i]]= list(iPCR= unlist(iPCR$barcode[[i]]), 
    gDNA= gDNA$barcode[gDNA[i] > 0], RNA= RNA$barcode[RNA[i] >0])
}
# with a venn diagram I can see how many barcodes are common within the 3 sets
pdf(file='vennDiagrams.pdf')
venn(barcode.list$HSP)
venn(barcode.list$p7)
venn(barcode.list$p10)
venn(barcode.list$p14)
venn(barcode.list$p15)
venn(barcode.list$p18)
venn(barcode.list$p20)
dev.off()
###########################################################################################

   


#To get the barcodes that match either in 2 or the 3 data sets:

  

  promoters= c( 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')
  gDNA.mapped= list()

  for (i in promoters) {
  gDNA.mapped$barcode[[i]]= gDNA$barcode[gDNA$barcode %in% iPCR$barcode[[i]]]
  gDNA.mapped$color[[i]]= iPCR$color[[i]][iPCR$barcode[[i]]%in%gDNA.mapped$barcode[[i]]]
  gDNA.mapped$counts[[i]]= gDNA[(gDNA$barcode %in% iPCR$barcode[[i]]),i]

  }
  
   mayor.que.cero= function (x) {
    x > 0
   }
   
  idx=lapply(gDNA.mapped$counts, FUN=mayor.que.cero)

  for (i in promoters) {

  gDNA.mapped$barcode[[i]]= gDNA.mapped$barcode[[i]][idx[[i]]]
  gDNA.mapped$color[[i]]= gDNA.mapped$color[[i]][idx[[i]]]
  gDNA.mapped$counts[[i]]= gDNA.mapped$counts[[i]][idx[[i]]]

  
}
 #I do the same for the RNA
 
RNA= read.delim('1610-TTA-count-min5.txt', sep='\t', as.is= TRUE, fill=TRUE, flush=TRUE)

names (RNA) = c('barcode', 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')
 
 
  promoters= c( 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')
  
  RNA.mapped= list()

  for (i in promoters) {
  RNA.mapped$barcode[[i]]= RNA$barcode[RNA$barcode %in% iPCR$barcode[[i]]]
  RNA.mapped$color[[i]]= iPCR$color[[i]][iPCR$barcode[[i]]%in%RNA.mapped$barcode[[i]]]
  RNA.mapped$counts[[i]]= RNA[(RNA$barcode %in% iPCR$barcode[[i]]),i]

  }
  
   mayor.que.cero= function (x) {
    x > 0
   }
   
  idx=lapply(RNA.mapped$counts, FUN= mayor.que.cero)

  for (i in promoters) {

  RNA.mapped$barcode[[i]]= RNA.mapped$barcode[[i]][idx[[i]]]
  RNA.mapped$color[[i]]= RNA.mapped$color[[i]][idx[[i]]]
  RNA.mapped$counts[[i]]= RNA.mapped$counts[[i]][idx[[i]]]

  }
  
    pdf(file='RNAplotsCounts.pdf')
    
    plot(x= as.factor(RNA.3match$color$p18), y= RNA.3match$counts$p18, 
    col= levels(as.factor(RNA.3match$color$p18)), ylab= 'counts', main= 'P18')
    
    plot(x= as.factor(RNA.mapped$color$p15), y= RNA.mapped$counts$p15, 
    col= levels(as.factor(RNA.mapped$color$p15)), ylab= 'counts', main= 'P15')
    
    plot(x= as.factor(RNA.mapped$color$HSP), y= RNA.mapped$counts$HSP, 
    col= levels(as.factor(RNA.mapped$color$HSP)), ylab= 'counts', main= 'HSP')
    
    plot(x= as.factor(RNA.mapped$color$p7), y= RNA.mapped$counts$p7, 
    col= levels(as.factor(RNA.mapped$color$p7)), ylab= 'counts', main= 'P7')
    
    plot(x= as.factor(RNA.3match$color$p10), y= RNA.3match$counts$p10, 
    col= levels(as.factor(RNA.3match$color$p10)), ylab= 'counts', main= 'P10')
    
    plot(x= as.factor(RNA.mapped$color$p14), y= RNA.mapped$counts$p14, 
    col= levels(as.factor(RNA.mapped$color$p14)), ylab= 'counts', main= 'P14')
    
    plot(x= as.factor(RNA.mapped$color$p20), y= RNA.mapped$counts$p20, 
    col= levels(as.factor(RNA.mapped$color$p20)), ylab= 'counts', main= 'P20')
    
    dev.off()
    
  #I make a list with all the barcodes that are common in the 3 sets
    
  promoters= c( 'HSP', 'p14', 'p20', 'p18', 'p10', 'p15', 'p7')
  RNA.3match= list()
 
 for ( i in promoters) {
   
    RNA.3match$barcode[[i]]= RNA.mapped$barcode[[i]][RNA.mapped$barcode[[i]]%in%gDNA.mapped$barcode[[i]]]
    RNA.3match$color[[i]]= RNA.mapped$color[[i]][RNA.mapped$barcode[[i]]%in%gDNA.mapped$barcode[[i]]]
    RNA.3match$counts[[i]]= RNA.mapped$counts[[i]][RNA.mapped$barcode[[i]]%in%gDNA.mapped$barcode[[i]]]
    
 } 

 
 
 
 
 
 