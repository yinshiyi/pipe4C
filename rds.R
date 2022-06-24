# load the package
require(rtracklayer) 
require("caTools")

# this genome file depending of which genome is used to generated 4C
require("BSgenome.Hsapiens.UCSC.hg38")
require("BSgenome.Hsapiens.UCSC.hg19")

#Get the chromosome info required for the bigwig
assign( 'genome', base::get(as.vector('BSgenome.Hsapiens.UCSC.hg38')))


#function to normalize the reads
norm4C <- function( readsGR, nReads=1e6, nTop=2, wSize=21 ) {
  readsGR$normReads <- 0
  sumTop <- sum( -sort( -readsGR$reads )[ 1:nTop ] )
  wNorm <- nReads/( sum( readsGR$reads )-sumTop )
  readsGR$normReads <- wNorm*readsGR$reads
  readsGR$norm4C <- runmean( x=readsGR$normReads, k=wSize, endrule="mean" )
  return( readsGR )
}

#function to make the bw file
makeBW<-function(file, OutFile, BSgenome, blind=TRUE, nReads=1e6, nTop=2, wSize=21){
  message("Converting RDS to BigWig file: ",file)
  GR<-readRDS(file)$reads
  if(blind==FALSE){
    GR <- subset( GR, type == 'non_blind')
    # non_blind type is 4c type
  }
  GR<- norm4C(GR, nReads=nReads, nTop=nTop, wSize = wSize)
  
  start(GR)<-GR$pos
  end(GR)<-GR$pos
  GR$score<-GR$norm4C
  
  
  
  seqinfo(GR)<-keepStandardChromosomes(seqinfo(genome))  
  
  
  export(GR, OutFile)
}





#Commando to convert the RDS to bw

makeBW(file="test.rds", OutFile="test.bw", BSgenome=genome)



