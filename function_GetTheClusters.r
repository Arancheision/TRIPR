
   #############################################################################
  #get.the.clusters is a function that clusters the integrations within
  #5 nucleotides proximity
  ##############################################################################
      
      # The input is a list containing:
                         #positions, it's a list containg as many elements as chromosomes
                         #Counts, it's a list containg as many elements as chromosome
      #Output is also a list containg
                         #Positions per chromosome
                         #Total countsCounts per chromosome
   

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
    
    return (list('integrations'= integrations, 'counts'= totalcounts))
}


