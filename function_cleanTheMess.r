
 ####################################################################################
 #Clean the mess is a function that returns a list with unique integrations points
 #aggregated by position
 ####################################################################################
 
 
  
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



 clean.the.mess= function(filename) {

        raw.data = read.table(filename, sep = "\t", as.is= TRUE)


        # to remove all the characters and leave only the barcode we use sub:
        raw.data$V1 = sub(".*:", "", raw.data$V1)

        # we remove the unmapped and low quality data,
        #making a cutoff 10 for MAPPING QUALITY (V5 column) and V2= 16 or 0

        sel=(raw.data$V2%in%c(0,16))&(raw.data$V5 >10)
        raw.data= subset(raw.data, sel)


        positions= list()
        counts= list()
        chrom.names= sort(unique(raw.data$V3))
        
        for (chrom in chrom.names) {
              per.chromosome =subset( raw.data, raw.data$V3 == chrom ) 
              counts.table =table ( per.chromosome$V4) 
              positions[[chrom]]= as.integer(names(counts.table))
              counts[[chrom]] = counts.table
        }

        return(list('positions'=positions, 'counts'=counts))
}




