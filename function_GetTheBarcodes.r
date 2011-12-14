  ##########################################################################
 #get.the.barcodes is a function that: aggregate the integrations by barcode,
                                #removes the colisions
                                #cluster the integrations
                                #and gets the real barcodes
 ###########################################################################

                         #input, is the bam file converted to text
                         #The output is a list containg: barcodes
                                                         #integration points
                                                         #counts

get.the.barcodes = function(filename) {

  raw.data = read.table(filename, sep = "\t", as.is= TRUE, flush= TRUE, fill= TRUE)

    # to remove all the characters and leave only the barcode we use sub:
    raw.data$V1 = sub(".*:", "", raw.data$V1)
    
    # we remove the unmapped and low quality data,
    #making a cutoff 10 for MAPPING QUALITY (V5 column) and V2= 16 or 0
    
    sel=(raw.data$V2%in%c(0,16))&(raw.data$V5 >10)
    raw.data= subset(raw.data, sel)


    # Replacement 1: quality -> unique ids.
    raw.data$V5 = 1:nrow(raw.data)
    
    #in per.barcode we aggregate the integrations by barcode, so for a given
    # barcode there might be more than one integration point
    
    raw.data = aggregate(raw.data, by=list(raw.data$V1), FUN=unique)
    
    #Replacement 2: unique IDs -> counts per barcode
    raw.data$V5= sapply(raw.data$V5, length)
    
    #to check wether the same barcode is present in two different positions,
    # further than 6bp
    
    check.spread= function (x) {
       sum(diff(x))>6
    }
    
    colisions= sapply(raw.data$V4, FUN= check.spread)
    # All cases in which one barcode is integrated in two different positions

    #I remove those:

   raw.data= raw.data[!colisions,]
   
   #there are also cases in which there is the same barcode in more than one chromosome:
   
   chr.colisions=sapply(raw.data$V3, FUN= length)
   
   
   #I remove those:
   
   raw.data= raw.data[!chr.colisions >1 ,]
   
   #I calculate the mean of the positions of each island:  
   
    raw.data$V4= sapply(raw.data$V4, FUN=mean)
    
    
    #To get the correct input for get.the.clusters,
    #I make a list containg: barcodes, positions and counts by chromosome:
    raw.data$V2= unlist(raw.data$V2)
    raw.data$V3= unlist(raw.data$V3)
    raw.data$V4= unlist(raw.data$V4) 

    barcodes= tapply(raw.data$V1, INDEX=raw.data$V3, FUN=I)
    positions= tapply(raw.data$V4, INDEX=raw.data$V3, FUN=I)
    counts= tapply(raw.data$V5, INDEX=raw.data$V3, FUN=I)
    
    listname= list(barcode= barcodes, position= positions, counts= counts)

    #to get the clusters, and the corresponding barcodes:

    distance = c()
    clust = list()
    ct = list()
    integrations= list()
    totalcounts= list()
    barcodes= list()
    chrom.names = c("2L", "2R", "3L", "3R", "4", "X")
    
    
    for (chrom in chrom.names) {
    
      if(length(listname$position[[chrom]])<= 1) {
        clust[[chrom]] = NA
      } else {
        distance = dist(listname$position[[chrom]])
        clust[[chrom]]= hclust(distance)
        ct[[chrom]]= cutree(clust[[chrom]], h=5)
        integrations[[chrom]]=  tapply(ct[[chrom]], X= listname$position[[chrom]], FUN= mean)
        totalcounts[[chrom]] = tapply (ct[[chrom]], X= listname$counts[[chrom]], FUN= I)
        barcodes [[chrom]]= tapply (ct[[chrom]], X=listname$barcode[[chrom]], FUN= unique)
      }
    }
    
    listname = list('integrations'= integrations, 'counts'= totalcounts, 'barcodes'= barcodes)


  #Due to mutations during PCR or seq errors, there is more than one barcode in the same cluster,
  #the 'real' one would be that with the highest count
  
  nintegrations= c()
  real.barcodes= list()
  real.counts= list()
  chrom.names = c("2L", "2R", "3L", "3R", "4", "X")

  choice = sapply(listname$counts, sapply, which.max) #To get the maxcount
  
  for (chrom in chrom.names) {
  
     nintegrations= length(listname$barcodes[[chrom]])
     
     
     
     for (i in 1:nintegrations) {

          real.barcodes[[chrom]][[i]]= listname$barcodes[[chrom]][[i]][choice[[chrom]][[i]]]
          real.counts[[chrom]][[i]]= listname$counts[[chrom]][[i]][choice[[chrom]][[i]]]
         
      }
 }
    return(list(integrations=listname$integrations, barcode= real.barcodes, counts= real.counts))

}

