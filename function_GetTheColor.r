    
    
    #####################################################################################
     #get.the.color is a function that overlaps the integrations with the color domains
    #####################################################################################
    
    
    
      
     #Input:
          #List of 2 elements (integrations and counts) each of them containg another 6
          #elements (chromosome)

     
          #the color domains should be previously load in the sesion
          # color.domains= read.delim('GSE22069_Drosophila_chromatin_domains.txt', as.is=TRUE)
          # it should be a data frame with 4 columns containig:
                                            # Seqname(chromosome name),
                                            #Start
                                            #End
                                            #Chromatin(color)
            #The file is modify in order to have the same names for the chromosomes
            #color.domains$seqname= sub('chr','', color.domains$seqname)

     #Output:
            #is a list with as many elements as chromosomes containg the color of each integration
            
            

get.the.color= function(listname,color.domains) {

    fly.chroms = c("2L", "2R", "3L", "3R", "4", "X")
    integration.color = list()

    for (chrom in fly.chroms) {

      color= rep(NA, length(listname$integrations[[chrom]])) # Modify on every iteration.
      this.domains = subset(color.domains, color.domains$seqname == chrom)
      #to use as index later on

      for (i in 1:nrow(this.domains))  {

          idx <- (listname$integrations[[chrom]] >= this.domains$start[i]) &
          (listname$integrations[[chrom]] <= this.domains$end[i])
          color[idx] <- this.domains$chromatin[i]
      }
      integration.color[[chrom]] = color #it stores the chromatin color in each iteration
    }

    return(integration.color)
}
