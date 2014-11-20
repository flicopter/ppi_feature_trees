args = commandArgs( trailingOnly=T )

directory <- args[1] 
result_filename <- args[2]

get_score <- function( file ) {
	d<-read.table( file, skip=4 )
	scores <- d$V7[1:2000]
	nscores<-(scores - mean(scores)) / sd(scores)
	result <- nscores[1]
        result
}

outfiles <- list.files( directory, pattern="*.out", full.names=TRUE )  
filenames <- basename( outfiles )
filenames <- gsub( ".out$", "", gsub( "^.*-","",filenames ) )  
scores <- sapply( outfiles, get_score )
results <- data.frame( Ligand=filenames, Score=scores )
row.names( results ) <- NULL
write.table( results, result_filename, sep="\t", quote=FALSE, row.names=FALSE )

