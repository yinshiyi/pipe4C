setwd(dirname(rstudioapi::getSourceEditorContext()$path))
a<-read.csv('4c_promoter.csv')

control<-a[which(a$samplename==a$control),]
samplelist<-a[which(a$samplename!=a$control),]

for (i in c(1:nrow(samplelist))) {
  print((samplelist$filelocation[i]))
  print((control$filelocation))
  print((samplelist$samplename[i]))
  print((samplelist$viewpoint[i]))
  }

  