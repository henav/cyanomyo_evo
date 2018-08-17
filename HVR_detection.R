#script to take blast results of GOS project, create recruitment plot to S-TIM4 genome, 
# add mutations on the genome to the plot. 
library(IRanges)

#1. make a table of reads/mates. 
#example for line in the original table (blast out)):
#JCVI_READ_1915587	938	91	242	104763	104612	TCCTTCATTTTCTGTTGCAAGTCAATGAGTTTGTCTGTTGTATCTGCGACTGCTTTGATTGTAGTTGCAGCGACTTCATATGCTCTTGCAGAATCTGATTCTTGTGCTAATTCTAATATACCATTCACTGCTTCTTGTCCCTTCTCAACTAA	TCTTTCATCTTTTTCTGAAGATCAATTAGTTTATCAGCAGTATCTGCTACGTTCTTAATCGTTGTTGCAGCAACTTCATACGCACGAGGATGATCTGATGCTTGTGCTACTTCCAATATACCATCCACTGCTTCTTGTCCTTTCATCACTAA	1e-20	107	152	75.66	115	37	115	0	0	16	16

blast_titles<-c("reads", "qlen", "qstart", "qend", "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "length", "pident", "nident", "mismatch", "positive", "gapopen", "gaps", "qcovs", "qcovhsp")
dataset <- read.table(file="", sep="\t",header=F,col.names = blast_titles)
pairs<-read.table(file="gos_reads_pairs.txt", header=TRUE, sep="\t") #file with the GOS pairs, from the GOS dataset

#2. A. filter the reads:
#     A1. merge the mates to the reads
dataset<-merge.data.frame(dataset,pairs, by.x="reads",by.y="reads", all.x = TRUE)   

#     A2. check the first conditon on the reads, add to column, Then check second condition, add to new column
dataset$read_first_cond <- ((dataset[,2]>=100 & dataset[,2]<300 & dataset[,3]<=20 & dataset[,4]>=(dataset[,2]-20)) | (dataset[,2]>=300 & dataset[,3]<=25 & dataset[,4]>=(dataset[,2]-25)))
dataset$read_second_cond <- (dataset[,12]>79.9 & (abs(dataset[,4]-dataset[,3])/dataset[,2])>0.8)

#     A3. check the two conditions for the MATES, then assign them to a new table that contains the mates and the condintions.
mates_data<-data.frame(dataset$mates,((dataset[,2]<300 & dataset[,3]<=20 & dataset[,4]>=(dataset[,2]-20)) | (dataset[,2]>=300 & dataset[,3]<=25 & dataset[,4]>=(dataset[,2]-25))), (dataset[,12]>79.9))
mates_data<-na.omit(mates_data)
colnames(mates_data) = c("mates","mate_first_cond","mate_second_cond")

#     A4. Merge the MATES data talbe with the dataset table, by read. so now, the info for the MATE of each read is added
dataset<-merge.data.frame(dataset,mates_data, by.x="reads",by.y="mates", all.x = TRUE)

#     A5. subset the reads acording to the conditions. 

newset<-(dataset[((dataset$read_first_cond) & (dataset$mate_first_cond) & (!is.na(dataset$mate_first_cond))) | 
        ((dataset$read_first_cond) & (dataset$mate_second_cond) & (!is.na(dataset$mate_first_cond))) | 
        ((dataset$read_second_cond) & (dataset$mate_first_cond) & (!is.na(dataset$mate_first_cond))) |
        ((dataset$read_first_cond) & (is.na(dataset$mates))),])

               
#2. B. replace the order of the start and end of subject, if it's on the reverese strnad (the length should be positive, for IRanges)

for (i in 1:nrow(newset)) {
	if (newset[i,6]<newset[i,5]) {
		temp_holder<-newset[i,6]
		newset[i,6]<-newset[i,5]
		newset[i,5]<-temp_holder
		rm(temp_holder)
	}
}

#3. do a recruitment plot.
#first, draw the plot area
x<-(1:175290)
y<-rep(1,175290)
pos1<-15500
pos2<-21500 #limits of the x axis in the genome
plot (x,y,ylim=c(0,100),xlim=c(pos1,pos2),type="n", ylab="Blast Identity",xlab="Position")
#now the actual recruited reads
segments(newset[,5],newset[,12],newset[,6],newset[,12], lwd=2,col="purple") #think about different colors for different sets, or color code acourding to %identity


#4. Add a coverage plot (using IRanges) 


ir<-IRanges(start<-newset[,5],end<-newset[,6])
cov <- coverage(ir)
cov <- as.vector(cov)
mat <- cbind(seq_along(cov), cov)
d <- diff(cov) != 0
mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
mat <- mat[order(mat[,1]),]
#lines(cov,col="blue", lwd=0.5)
plot(cov,col="blue", lwd=1.5, type="l", ylab="X-fold coverage", xaxt="n",xlab="")
lab<-c(0,10,20,30,40,50,60,70,80,90,100) # im not sure this is working
axis(4,at=lab, labels=lab)

#5. A. plot the mutations over the recruitment plot / plot the mutations on a separate panel 

MED4<-c()#MED4 populations files (files with the mutation frequencies in each population)...
MIT9515<-c()#MIT9515 populations files...
WH8102<-c()#WH8102 populations files...
colors<-c("light green","dark green","red")
bacs<-list(MED4,MIT9515,WH8102)

plot (x,y,ylim=c(0,100),xlim=c(pos1,pos2),type="n", ylab="Mutation Frequency",xlab="Position")

for (i in 1:length(bacs)) {
	for (file in 1:length(bacs[[i]])) {
		dataset1 <-read.table(bacs[[i]][file], sep="\t", header=TRUE)
		dataset1.del<-dataset1[grep("DEL", dataset1$type), ]
		dataset1.ins<-dataset1[grep("INS", dataset1$type), ]
		dataset1.snp<-dataset1[grep("SNP", dataset1$type), ]
		points (dataset1.snp[,6],dataset1.snp[,7]*100, bg=colors[i],col="black",pch=21,cex=1)
		points(dataset1.ins[,6],dataset1.ins[,7]*100, bg=colors[i],col="black",pch=23,cex=1)
		points(dataset1.del[,6],dataset1.del[,7]*100, bg=colors[i],col="black",pch=23,cex=1)
		
		rm(dataset1)
		rm(dataset1.snp)
		rm(dataset1.ins)
		rm(dataset1.del)
	}
}

#6. define areas which are hvr

treshold<-median(cov)/5
tot_count<-0
counter<-0
newset<-as.table(first,nrow=1)


for (pos in c(pos1:pos2)) {
	if (cov[pos]<=treshold) {
		counter<-counter+1
	}
	if (cov[pos]>treshold)  {
		if (counter>500) { #can change the counter - to look for HVR in different minimal length
			print (paste("counter:", counter, "position:", pos), quote=FALSE)
		  tot_count<-tot_count+counter
			x.cord<-c(pos-counter-1,pos-counter-1,pos-1,pos-1)
			y.cord<-c(0,110,110,0)
			polygon(x.cord,y.cord, col="#0101ff22") 
		}
		counter=0
	}	
}
tot_count
abline(h=treshold, lwd=0.5)
