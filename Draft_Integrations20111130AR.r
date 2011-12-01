  ####################################################################################
 #Clean the mess is a function that returns a list with unique integrations points
 #aggregated by position
 ####################################################################################
 
 #The .bam file has been transform into a text file containing only the first 5 columns that are
 #the ones we need for the mapping, the comands to do that are:
            #samtools view <filename> | awk -v OFS="\t" '{print $1, $2, $3, $4, $5}' > 'outputfilename'


  # the text file contains:

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

        return(listname=list('positions'=positions, 'counts'=counts))
   }

 #############################################################################
 #get.the.clusters is a function that clusters the integrations within
 #5 nucleotides proximity
 ##############################################################################
      
 # The input is a list containing:
    #positions, it's a list containg as many elements as chromosomes
    #Counts, it's a list containg as many elements as chromosome


    get.the.clusters = function(listname) {


        distance = c()
        clust = list()
        ct = list()
        integrations= list()
        totalcounts= list()
        chrom.names= sort(names(listname$positions))


        for (chrom in chrom.names) {

          if(length(listname$positions[[chrom]])<= 1) {
          clust[[chrom]] <- NA
          } else {
          distance = dist(listname$positions[[chrom]])
          clust[[chrom]]= hclust(distance)
          ct[[chrom]]= cutree(clust[[chrom]], h=5)
          integrations[[chrom]]=  tapply(ct[[chrom]], X= listname$positions[[chrom]], FUN= mean)
          totalcounts[[chrom]] = tapply (ct[[chrom]], X= listname$counts[[chrom]], FUN= sum)
        }
      }

      return (cluster.list=list('integrations'= integrations, 'counts'= totalcounts))
    }
    

        #Output is  a list containg
            #Integrations, list of as many elements as chromosomes containing the integrations
            #Counts, list of as many elements as chromosomes containing the counts of each
            #integration

