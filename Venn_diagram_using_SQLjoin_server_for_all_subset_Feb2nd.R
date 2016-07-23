#!/usr/bin/ Rscript

args <- commandArgs(TRUE)
# el<-read.table(args[1], sep="\t", header=T,stringsAsFactors = FALSE)
#data1<-read.csv(args[1], header=F,stringsAsFactors = FALSE)
#data2<-read.csv(args[2], header=F,stringsAsFactors = FALSE)
#data3<-read.csv(args[2], header=F,stringsAsFactors = FALSE)
#data4<-read.csv(args[2], header=F,stringsAsFactors = FALSE)
#data5<-read.csv(args[2], header=F,stringsAsFactors = FALSE)

library(stringr)
# str_extract('A12_FC_IGG.fa_unique_CDR1_CDR3.csv','CDR.*\\.')

name<-paste(str_extract(args[1],'^\\w*?_'),str_extract(args[1],'CDR.*\\.'),'_all_subset',sep='')

format_data<-function(nu) {
data<-read.csv(args[i], header=F,stringsAsFactors = FALSE)

data$V3<-str_replace(str_replace(str_extract(args[i],'\\w*?_IGG'),'_IGG',''),'A22_','')

data<-data[data$V1!='',]
colnames(data)<-c('dna','count','subset')

return (data)

}

for (i in 1:5){
print (i)
assign(paste('data',i,sep=''),format_data(i))
}

print (head(data1))
print (head(data2))
print (head(data3))
print (head(data4))
print (head(data5))


data<-Reduce(function(...) merge(..., by='dna', all=T), list(data1,data2,data3,data4,data5))

head(data)
library(gplots)

data2<-data[,c(3,5,7,9,11)]
  library(plyr)
  
  head(data2)
  str(data2)
  data2[is.na(data2)]<-0
  table(data2)
  result<-as.vector(apply(data2,2,function(x) max(as.character(x))))
  data2<-data[,c(3,5,7,9,11)]
  colnames(data2)<-result
  head(data2)
  data2[!is.na(data2)]<-1
  data2[is.na(data2)]<-0
  str(data2)
  for (i in 1:dim(data2)[2]) {
    data2[,i]<-as.numeric(data2[,i])
  }
  #pdf(file="A20_SHM_IgG_nonredudant_knokcinOnly.pdf", width=10, height=8, pointsize=12)
  pdf(file=paste('venn_diagram_all_subset',name,'.pdf',sep=''),width=10, height=8, pointsize=12)
  venn(data2)
  dev.off()
  venn_table<-venn(data2)
  write.csv(venn_table,paste('data_for_venn_diagram_all_subset_',name,'.csv',sep=''),quote=F)

  write.csv(data,paste('data_for_venn_diagram_all_subset_',name,'_raw.csv',sep=''),quote=F)

