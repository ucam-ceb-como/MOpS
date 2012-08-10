inputfile="pahtest2-bintree-serializer-postprocess-PAH(0.006s).csv"

psl<-read.csv(inputfile, header=TRUE)
# 7th colunm is mass of individual PAH (u)
c<-psl[7]
# remove the NA form c before calculation
c<-c[!is.na(c)]
# calculate the average size of the ten largest particle.
mea<-mean(c)

#write.table(mea,"pahtest2/stats.csv",na = "",append=FALSE,row.names = FALSE,quote = FALSE,sep = ",",col.names="median")
write.table(mea,"stats.csv",na = "",append=FALSE,row.names = FALSE,quote = FALSE,sep = ",",col.names="median")
