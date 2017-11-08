# MHH - Cluster analysis to work with index file- for DNA PAINT Clusters. 16/02/17

#read the data file which runs various programs
rows <- read.table("/Users/Mathew/Desktop/rplay.txt")


# Make a function to process each file

processFiles <- function(files) {
	
	filestoload <- paste(files,"fitresults.txt",sep="")			# File-name of the files to load. 
	filestosave <- paste(files,"cluster.txt",sep="")			# File-name of the files to save. 
	
	load<-filestoload											# Load the files. 
	
	tetf<- read.table(load,header=T,dec='.',sep='\t')			# Assign rows to tetf.
	
	#install or load package
  	if(!require(fpc)){install.packages('fpc')}; require('fpc')	# Install DBSCAN package.
  		
  	#cluster
  	db <- dbscan(tetf[,12:13], eps=2, MinPts = 4)				# Run cluster analysis. 
  		
  	plot(db,tetf[,12:13])										# Show plots while running. 
		
	write.table(db$cluster, file=filestosave, sep = ",",col.names = NA, qmethod = "double")		# Save the data. 
		
	}


# Find all .csv files
i=0
for(i in rows)
{
files <- rows[i,1]

print(files)

}

# Apply the function to all files.
result <- sapply(files,processFiles)







