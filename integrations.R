# get.integrations, is a function that clusters the integrations within 5 nucleotides proximity 
   #and gives a color to each of them 
    
      # the input is a text file containing only the first 5 columns of the .bam file:
      
          #V1 contains a given number of characters with the barcode at the end
             #this column will be modify to contain only the barcode sequence

          #V2 (FLAG) column,
          #Bit 0 = The read was part of a pair during sequencing
          #Bit 1 = The read is mapped in a pair
          #Bit 2 = The query sequence is unmapped
          #Bit 3 = The mate is unmapped
          #Bit 4 = Strand of query (0=forward 1=reverse)
          #for 16 =2^4 = 10000 = FALSE, FALSE, FALSE,FALSE, reverse
          #for 4 = 2^2 = 100 = FALSE, FALSE, TRUE= seq unmapped
          #for 20= 2*5= 10100= FALSE, FALSE, TRUE,FALSE, reverse  = seq unmapped
            # It will be modify to contain only 16 or 0 values (Rev or Fw)
            
          #V3 Chromosome name
          
          #V4 Starting position

          #V5 Mapping quality
      
      
      # the output is a list containg: the integrations (clustered) per chromosome, 
                                       #color of each integration per chromosome
                                       #counts per integration per chromosome
                                       
get.integrations= function (filename){

   raw.data = read.table(filename, sep = "\t", as.is= TRUE)


    # to remove all the characters and leave only the barcode we use sub:
    raw.data$V1 = sub(".*:", "", raw.data$V1)  

    # we remove the unmapped and low quality data, 
    #making a cutoff 10 for MAPPING QUALITY (V5 column) and V2= 16 or 0

    sel=(raw.data$V2%in%c(0,16))&(raw.data$V5 >10)
    raw.data= subset(raw.data, sel)

    # In order to know the number of barcodes we have after cleaning up the data:
    barcodes=unique(raw.data$V1)

    #Counts per barcode:
    counts= table(raw.data$V1)


    # Replacement 1: quality -> unique ids.
    raw.data$V5 = 1:nrow(raw.data)

    #in per.position we have aggregate the integrations by position, for a given integration point there is more than one
    #barcode, in most cases this is due to mutations in the barcode

    per.position = aggregate(raw.data, by=list(raw.data$V4), FUN=unique) 

    inscount = sapply(per.position$V5, length)
    # Replacement 2: unique ids -> number of ids.
    per.position$V5 = inscount


    #in per.barcode we have aggregate the integrations by barcode, so for a given barcode there might be more than 
    #one integration point
    per.barcode = aggregate(raw.data, by=list(raw.data$V1), FUN=unique)
    bbcounts = sapply(per.barcode$V5, length)
    # Same as replacement 2 above.
    per.barcode$V5= bbcounts





 


    #Because per.position has a weird format (it's a data frame which first column is a list) I'll make a
    #real list with the data

    Ch2L= per.position[(per.position$V3== '2L'),5:6]
    Ch3L= per.position[(per.position$V3== '3L'),5:6]
    Ch2R= per.position[(per.position$V3== '2R'),5:6]
    Ch3R= per.position[(per.position$V3== '3R'),5:6]
    Ch4= per.position[(per.position$V3== '4'),5:6]
    ChX= per.position[(per.position$V3== 'X'),5:6]

    positions= list('2L'= Ch2L$V4, '3L'= Ch3L$V4, '2R'= Ch2R$V4, '3R'=Ch3R$V4, '4'= Ch4$V4, 'X'= ChX$V4)
    counts=  list('2L'= Ch2L$V5, '3L'= Ch3L$V5, '2R'= Ch2R$V5, '3R'=Ch3R$V5, '4'= Ch4$V5, 'X'= ChX$V5)

    # To group the integrations that differ by 5 nucleotides, we clust them, and use a for loop
    # for all chromosomes
 
    distance = c()
    clust = list()
    ct = list()
    integrations= list('2L'= c(), '3L'= c(), '2R'= c(), '3R'= c(), '4'= c(), 'X'= c())
    totalcounts= list('2L'= c(), '3L'= c(), '2R'= c(), '3R'= c(), '4'= c(), 'X'= c())

    for (i in 1:length(positions)) {
        distance = dist(positions[[i]])
        clust[[i]]= hclust(distance)
        ct[[i]]= cutree(clust[[i]], h=5)
        integrations[[i]]=  tapply(ct[[i]], X= positions[[i]], FUN= mean)
        totalcounts[[i]] = tapply (ct[[i]], X= counts[[i]], FUN= sum)
    }

    # To plot the integrations per chromosome:

    pdf(file='integrationsPlot.pdf')

    for (i in 1:6) {
    
        plot (x= integrations[[i]], y= totalcounts[[i]], type= 'h', main= paste('Chromosome',names(integrations[i])), xlab= 'Position (bp)', ylab='Counts')
    }

    dev.off()


    #Overlap with Chromatin colors
      #load the data, consisting in a data frame with 4 columns
      # Seqname(chromosome name), Start, End, Chromatin(color)

    color.domains= read.delim('GSE22069_Drosophila_chromatin_domains.txt', as.is=TRUE)

    #Remove the Chr from the seqname, in order to get the same name as in our 
    # integration list
    color.domains$seqname= sub('chr','', color.domains$seqname)


    #To asign a color for each integration:


    fly.chroms = c("2L", "2R", "3L", "3R", "4", "X") #Creates a vector with 
    #elements equals to our seqnames in the integration and color.domains data
    integration.color = list()

    for (chrom in fly.chroms) {        
            
      color= rep(NA, length(integrations[[chrom]])) # Modify on every iteration.
      this.domains = subset(color.domains, color.domains$seqname == chrom)
      #to use as index later on

      for (i in 1:nrow(this.domains))  {

          idx <- (integrations[[chrom]] >= this.domains$start[i]) &
          (integrations[[chrom]] <= this.domains$end[i])
          color[idx] <- this.domains$chromatin[i]
      }
      integration.color[[chrom]] = color #it stores the chromatin color in each iteration
    }

   
    pdf(file='integrationsVsColorPlot.pdf')
    
    for (chrom in fly.chroms) {
        plot (x= integrations[[chrom]], y= totalcounts[[chrom]], type= 'h',
        xlab= 'Position (bp)', ylab='Counts', col= integration.color[[chrom]], 
        main= paste('Chromosome',names(integrations[chrom])))
    }
    dev.off()



    return(list('integrations'=integrations, 'integration.color'=integration.color,
     'totalcounts'=totalcounts))
} 
