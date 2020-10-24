
## Date - 20201008
## Find min, max, mean and mode for each coloum and for the row types

rm(list=ls())

## set directories

input <- "C:/Users/Suresh/Desktop"
output <- "C:/Users/Suresh/Desktop/bala_sample_output"

# read input file
setwd(input)
list.files()

example.file <- read.table("test_fmt_for_suresh.txt", header = T, sep = " ", stringsAsFactors = F)
names(example.file) <- gsub("^X", "", names(example.file))
head(example.file);nrow(example.file);ncol(example.file)

## calculate min, max, mean and mode for each column
number.functions <- c("min", "max", "mean", "median")
out.file <- NULL

j = 1
for(j in 1:length(unique(example.file[["Sampletype"]]))){
		uni.type <- grep(unique(example.file[["Sampletype"]])[j], example.file$Sampletype, value = T)
		tmp.type <- example.file[example.file[["Sampletype"]] %in% uni.type,]
				
		##result.type <- as.data.frame(result.type)
		##namesresult.type[["Sampletype"]] <- unique(uni.type)

		result.type <- NULL	
		result.type <- apply(tmp.type[2:length(tmp.type)], 2, min)
		result.type <- rbind(result.type, apply(tmp.type[2:length(tmp.type)], 2, max))
		result.type <- rbind(result.type, apply(tmp.type[2:length(tmp.type)], 2, mean))
		result.type <- rbind(result.type, apply(tmp.type[2:length(tmp.type)], 2, median))
		result.type <- as.data.frame(result.type)
		rownames(result.type) <- paste(uni.type, number.functions, sep = "_")
		#print(head(result.type))

	#out.file <- result.type
	out.file <- rbind(out.file, result.type)
	
}

out.file <- t(out.file)
print(out.file)
nrow(out.file)
rownames(out.file)
colnames(out.file)

##save output file and save to the desired directory

setwd(output)
write.table(out.file, file = "sample.file_bala.txt", sep = "\t", row.names = T, col.names = T)
				
				
			