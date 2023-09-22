library(devtools)
dev_mode()
install_github('mike-lawrence/ez')
library(ez)
Nperm = 10
SpanSpace_Word = read.csv('~/Documents/writen_paper/NER2020/Visuallization_matlab/Code/datainfo.csv',header=F)
SpanSpace_Word[,c(2:5)] <- lapply(SpanSpace_Word[,c(2:5)] , factor)
colnames(SpanSpace_Word)<-list('Contrast','region','subject','metric','task')
ezP_Word <- ezPerm(SpanSpace_Word,dv= Contrast,wid=.(region),within = .(subject, metric, task),perms = Nperm)

SpanSpace_Word = read.csv('~/Documents/writen_paper/NER2020/Visuallization_matlab/Code/datainfo.csv',header=F)
SpanSpace_Word[,c(2:4)] <- lapply(SpanSpace_Word[,c(2:4)] , factor)
colnames(SpanSpace_Word)<-list('Contrast','region','subject','metric')
# ezP_Word <- ezPerm(SpanSpace_Word,dv= Contrast,wid=.(region),within = .(subject, metric),perms = Nperm)
ezA_Word <- ezANOVA(SpanSpace_Word,dv= Contrast,wid=subject,within=.(metric), between=.(region),return_aov=T,type=3,detailed=T)



library(ez)
Nperm = 10
SpanSpace_contrast = read.csv('~/Documents/writen_paper/NER2020/Visuallization_matlab/Code/contrast_info.csv',header=F)
SpanSpace_contrast[,c(2:4)] <- lapply(SpanSpace_contrast[,c(2:4)] , factor)
colnames(SpanSpace_contrast)<-list('region','subject','metric','task')
ezP_contrast <- ezPerm(SpanSpace_contrast,dv= region,wid=.(subject),within = .(metric, task),perms = Nperm)

colnames(SpanSpace_contrast)<-list('region','subject','metric','task')
ezP_contrast <- ezPerm(SpanSpace_contrast,dv= region,wid=.(subject),within = .(metric, task),perms = Nperm)


library(ez)
Nperm = 10
SpanSpace_Word = read.csv('~/Downloads/NumPositive_Word.csv',header=F)
SpanSpace_Word[,c(2:4)] <- lapply(SpanSpace_Word[,c(2:4)] , factor)
colnames(SpanSpace_Word)<-list('NumPositive','Block','Filter','Subject')
ezP_Word <- ezPerm(SpanSpace_Word,dv=NumPositive,wid=.(Subject),within = .(Block,Filter),perms = Nperm)
