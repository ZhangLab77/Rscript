#ggplot
{
  ggplot(df, aes(gp, y)) +
    geom_point() +
    geom_point(data = ds, aes(y = mean), colour = 'red', size = 3)
  
  # Same plot as above, declaring only the data frame in ggplot().
  # Note how the x and y aesthetics must now be declared in
  # each geom_point() layer.
  ggplot(df) +
    geom_point(aes(gp, y)) +
    geom_point(data = ds, aes(gp, mean), colour = 'red', size = 3)
  
  # Alternatively we can fully specify the plot in each layer. This
  # is not useful here, but can be more clear when working with complex
  # mult-dataset graphics
  ggplot() +
    geom_point(data = df, aes(gp, y)) +
    geom_point(data = ds, aes(gp, mean), colour = 'red', size = 3) +
    geom_errorbar(
      data = ds,
      aes(gp, mean, ymin = mean - sd, ymax = mean + sd),
      colour = 'red',
      width = 0.4
    )
  
  ggplot(datas,aes(x=Group,y=trans*100,fill=TC))+ #”fill=“设置填充颜色
    stat_boxplot(datas,aes(x=Group,y=trans*100,fill=Group),geom = "errorbar",width=0.15,aes(color="black"))+
    stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+#由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
    #geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
    geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
    #scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色
    #scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色
    ggtitle("")+ #设置总的标题
    theme_bw()+ #背景变为白色
    theme(legend.position="none", #不需要图例
          axis.text.x=element_text(angle = 45,colour="black",family="Times",size=14,hjust = 1), #设置x轴刻度标签的字体属性
          axis.text.y.left =element_text(colour="black",family="Times",size=12), #设置x轴刻度标签的字体属性
          axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
          axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
          axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
          plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
          panel.grid.major = element_blank(), #不显示网格线
          plot.margin = unit(c(1,2,0.2,0.5),'lines'))+ #设置margin，类似于par(mar)
    ylab("#readsInAlu/#readsInTpromoter")+xlab("")
  
  
  
  ggplot(datas,aes(x=Group,y=trans*100,fill=Group))+ #”fill=“设置填充颜色
    stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
    #geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
    geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
    #scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色
    #scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色
    ggtitle("")+ #设置总的标题
    theme_bw()+ #背景变为白色
    theme(legend.position="none", #不需要图例
          axis.text.x=element_text(angle = 45,colour="black",family="Times",size=14,hjust = 1), #设置x轴刻度标签的字体属性
          axis.text.y.left =element_text(colour="black",family="Times",size=12), #设置x轴刻度标签的字体属性
          axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
          axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
          axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
          plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
          panel.grid.major = element_blank(), #不显示网格线
          plot.margin = unit(c(1,2,0.2,0.5),'lines'))+ #设置margin，类似于par(mar)
    ylab("#readsInAlu/#readsInTpromoter")+xlab("")
  
  
}


#HEK293 enhancer, m5c
{
  setwd("/data/baix/bxDownload/hg19/HEK293/m5c/")
  
  xfile=as.matrix(read.table(file="H3k4me1_h3k27ac.bed",header = F))
  #chr19   29352727        29357209        Peak_1  1000    .       
  #chr19   29354140        29356877        Peak_3242       110   .29.29526        89.70232        86.26564        1392
  
  sortM=apply(xfile[,c(2:3,8:9)],1,function(x){sort(as.numeric(x))})
  
  len=sortM[,3]-sortM[,2]
  flag=which(len>=200)
  out=cbind(xfile[flag,1],sortM[flag,2:3])
  write.table(out,file ="me1_k27ac.intersect.txt",quote =F,sep = "\t",row.names = F,col.names = F )
  
  summary(len)
}

#DESeq2
{
  library(DESeq2)
  data=as.matrix(read.table("Tumor.control.Rep1andRep2.count",header = F)) ;  dataL=nrow(data)
  datas=matrix(as.numeric(data[1:(dataL-5),-1]),nrow=nrow(data[1:(dataL-5),-1]))
  datas=`row.names<-`(datas,data[1:(dataL-5),1])
  condition <- factor(c("Tumor","Control","Tumor","Control"))
  dds <- DESeqDataSetFromMatrix(datas, DataFrame(condition), ~ condition)
  dds <- dds[ rowSums(counts(dds)) > 1, ]   #过滤low count数据
  dds <- DESeq(dds)     #差异分析
  res <- results(dds)   #用result()函数获取结果
  summary(res)  #summary()函数统计结果
  res=res[order(res$padj),]  ;   resdata=merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
  
}

#2020.12.11
#FANTOM5
{
  setwd("/data/baix/bxDownload/hg19/FANTOM5/")
  a=as.matrix(read.table(file="human_permissive_enhancers_phase_1_and_2_expression_count_matrix_sampleName.txt",header = F))
  #Id CNhs14406 CNhs14407 CNhs14408 CNhs14410  ....
  
  b=as.matrix(read.table(file="Human.sample_name2library_id.txt",header = F))
  #acute_myeloid_leukemia_(FAB_M5)_cell_line:THP-1_(fresh) CNhs10722
  
  
  flag=c()
  for(i in 1:nrow(b)){
    flagi=which(a==b[i,2])
    if(length(flagi)!=0){
      flag=rbind(flag,flagi)
    }else{
      flag=rbind(flag,0)
    }
  }
  write.table(cbind(b,flag),file ="Human.sample_name2library_id_sampleColNum.txt",quote =F,sep = "\t",row.names = F,col.names = F )
  
}


#2020.12.13
#kidney circRNA, FANTOM5 enahncer
{
  setwd("~/bxDownload/hg19/kidney/circRNA-seq/GSE108735/")
  xfile=as.matrix(read.table(file="kidney.Enhancer_circRNA.txt",header=F,sep = "\t"))
  #chr1    1611804 1612115 0       2       0       0       
  #chr1    1601102 1666274 .       .       -       circBase        NM_001787       CDK11B    sense overlapping       65172
  
  len1=as.numeric(xfile[,2])-as.numeric(xfile[,9])
  len2=as.numeric(xfile[,3])-as.numeric(xfile[,10])
  flag1=which(abs(len1)<500 & abs(len2)<500 & xfile[,12]!="exon")
  
  flag2=which(len1<0 & len2>0 & xfile[,12]!="exon")
  out=xfile[c(flag1,flag2),]
  
  write.table(out,file ="kidney.Enhancer_candidate_circRNA.txt",quote =F,sep = "\t",row.names = F,col.names = F )
  
  eRNA.len=as.numeric(out[,3])-as.numeric(out[,2])
  hist(eRNA.len)
  
}

#plot fig2.eAlu_tAlu.age.pdf fig2.eAlu_tAlu.age.png
{
  #GetAluDgyperPvalue.R
  setwd("~/GM12878/2newGetEPpairsInTAD/CAGEhead/")
  afile <- as.matrix(read.table(file="/data/baix/bxDownload/hg19/Alu/AluInhg19_Repbase.bed",header=F,sep="\t"))
  #chr1 26791 27053  AluSp SINE/Alu +
  alu=table(afile[,4])
  aluN=rownames(alu)
  afileL=dim(afile)[1]
  
  #efile <- as.matrix(read.table(file="./eRNA.eAlu.txt",header=F,sep="\t"))
  efile <- as.matrix(read.table(file="~/GM12878/2newGetEPpairsInTAD/CAGEhead+-/GM12878.eRNA_eAlu.txt",header=F,sep="\t"))
  
  #chr10   103663844       103664554       GMIneRNA+389    GMIneRNA+       +       chr10   103664429       103664737       AluSc   SINE/Alu        +
  efile=as.matrix(efile)
  efileL=dim(efile)[1]
  balu=table(efile[,11])
  
  baluN=rownames(balu)
  a=round(alu/dim(afile)[1],6)
  b=round(balu/efileL,6)
  
  eoutp=c()
  for(i in 1:length(aluN)){
    bflag=which(baluN==aluN[i])
    if(length(bflag)!=0){
      p1=phyper(balu[bflag]-1,alu[i],afileL-alu[i],efileL,lower.tail = F)
      eoutp=rbind(eoutp,c(aluN[i],alu[i],balu[bflag],a[i],b[bflag],b[bflag]/a[i],round(p1,4)))
    }
  }
  
  subAlu=c("AluJo","AluJb","AluSx","AluSq","AluSg","AluSc","AluSp","AluY","AluYb8","AluYa5")
  Afile=afile[which(afile[,4]%in%subAlu),]
  alu=table(Afile[,4])
  aluN=rownames(alu)
  AfileL=dim(Afile)[1]
  Efile=efile[which(efile[,11]%in%subAlu),]
  EfileL=dim(Efile)[1]
  balu=table(Efile[,11])
  
  baluN=rownames(balu)
  a=round(alu/dim(Afile)[1],6)
  b=round(balu/EfileL,6)
  
  eoutp=c()
  for(i in 1:length(aluN)){
    bflag=which(baluN==aluN[i])
    if(length(bflag)!=0){
      p1=phyper(balu[bflag]-1,alu[i],AfileL-alu[i],EfileL,lower.tail = F)
      eoutp=rbind(eoutp,c(aluN[i],alu[i],balu[bflag],a[i],b[bflag],b[bflag]/a[i],round(p1,4)))
    }
  }
  
  
  bfile=read.table(file="./targetProteinCoding.tAlu.txt",header=F,sep="\t")
  #chr10   101684767       101687767       ENSG00000227695.1       DNMBP-AS1       +       antisense       chr10   101685168       101685471     AluSx1   SINE/Alu        -
  
  bfile=as.matrix(bfile)
  bfileL=dim(bfile)[1]
  balu=table(bfile[,11])
  baluN=rownames(balu)
  a=round(alu/dim(afile)[1],6)
  b=round(balu/bfileL,6)
  
  toutp=c()
  for(i in 1:length(aluN)){
    bflag=which(baluN==aluN[i])
    if(length(bflag)!=0){
      p1=phyper(balu[bflag]-1,alu[i],afileL-alu[i],bfileL,lower.tail = F)
      toutp=rbind(toutp,c(aluN[i],alu[i],balu[bflag],a[i],b[bflag],b[bflag]/a[i],round(p1,4)))
    }
  }
  
  subAlu=c("AluJo","AluJb","AluSx","AluSq","AluSg","AluSc","AluSp","AluY","AluYb8","AluYa5")
  
  rfile=eoutp #AluJb   117273  386     0.097995        0.159702        0
  pvalue=c()
  for(i in 1:length(subAlu)){
    flag=which(rfile[,1]==subAlu[i])
    if(length(flag)!=0){
      pvalue=c(pvalue,rfile[flag,7])
    }else{
      pvalue=c(pvalue,1.1)
    }
  }
  epvalue=as.numeric(pvalue)
  
  rfile=toutp#AluJb   117273  386     0.097995        0.159702        0
  pvalue=c()
  for(i in 1:length(subAlu)){
    flag=which(rfile[,1]==subAlu[i])
    if(length(flag)!=0){
      pvalue=c(pvalue,rfile[flag,7])
    }else{
      pvalue=c(pvalue,1.1)
    }
  }
  tpvalue=as.numeric(pvalue)
  
  
  #version5
  { 
    library(plotrix)
    layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2),ncol = 2,byrow = T))
    par(mar=c(5,5,2,1),cex=1)
    tp=-log10(tpvalue)
    tp[c(1,3)]=3.5
    plot(1:length(subAlu),-log10(epvalue),lwd=2,type="o",col="tan1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,3.5),main="",ylab = "-log10(pvalue)")
    par(new=T)
    plot(1:length(subAlu),tp,type="o",lwd=2,col="skyblue1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,3.5),ylab="")
    axis.break(2, 3.4)
    legend("topright",inset=0.01,legend = c("Alu on eRNA","Alu on Target"),col = c("tan1","skyblue1"),lty = c(1, 1), bg = "gray90",cex=0.9,pch=21)
    axis(1,at=1:length(subAlu),labels=subAlu, las=2,lty = 1, lwd = 1)
    axis(2,at=c(seq(0,3,0.5),3.5),las=2,labels=c(seq(0,3,0.5),">3"), lwd =1,lty = 1)
    
    par(mar=c(0,5,0,1))
    plot(seq(0,1,0.5),seq(0,1,0.5),type = "n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    arrows(0, 0.85, 1, 0.85, col = 1,lwd=5,code = 2,length=0.2)
    text(0.05,0.55,"80myr",cex=1)
    text(0.95,0.55,"5myr",cex=1)
    
  }  
  #version6
  { 
    
    subAlu=c("AluJo","AluJb","AluSx","AluSq","AluSg","AluSc","AluSp","AluY","AluYb8","AluYa5")
    
    rfile=eoutp #AluJb   117273  386     0.097995        0.159702        0
    e.odds=c()
    for(i in 1:length(subAlu)){
      flag=which(rfile[,1]==subAlu[i])
      if(length(flag)!=0){
        e.odds=c(e.odds,rfile[flag,6])
      }else{
        e.odds=c(e.odds,0)
      }
    }
    e.odds=as.numeric(e.odds)
    
    rfile=toutp#AluJb   117273  386     0.097995        0.159702        0
    p.odds=c()
    for(i in 1:length(subAlu)){
      flag=which(rfile[,1]==subAlu[i])
      if(length(flag)!=0){
        p.odds=c(p.odds,rfile[flag,6])
      }else{
        p.odds=c(p.odds,0)
      }
    }
    p.odds=as.numeric(p.odds)
    
    layout(matrix(c(1,1,1,1,2,2,2,2,2,2,2,2,3,3),ncol = 2,byrow = T))
    par(mar=c(0,5,2,1),cex=1)
    
    plot(1:length(subAlu),e.odds,lwd=2,type="o",col="tan1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,1.5),main="",ylab = "Odds ratio")
    par(new=T)
    plot(1:length(subAlu),p.odds,type="o",lwd=2,col="skyblue1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,1.5),ylab="")
    
    #legend("topright",inset=0.01,legend = c("Alu on eRNA","Alu on Target"),col = c("tan1","skyblue1"),lty = c(1, 1), bg = "gray90",cex=0.9,pch=21)
    #axis(1,at=1:length(subAlu),labels=subAlu, las=2,lty = 1, lwd = 1)
    axis(2,at=seq(0,1.5,0.5),las=2,labels=seq(0,1.5,0.5), lwd =1,lty = 1)
    abline(h=1,lty=2)
    
    
    library(plotrix)
    #layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2),ncol = 2,byrow = T))
    par(mar=c(5,5,0,1),cex=1)
    tp=-log10(tpvalue)
    tp[c(1,3)]=3.5
    plot(1:length(subAlu),-log10(epvalue),lwd=2,type="o",col="tan1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,3.6),main="",ylab = "-log10(pvalue)")
    par(new=T)
    plot(1:length(subAlu),tp,type="o",lwd=2,col="skyblue1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,3.6),ylab="")
    axis.break(2, 3.4)
    legend("topright",inset=0.01,legend = c("Alu on eRNA","Alu on Target"),col = c("tan1","skyblue1"),lty = c(1, 1), bg = "gray90",cex=0.9,pch=21)
    axis(1,at=1:length(subAlu),labels=subAlu, las=2,lty = 1, lwd = 1)
    axis(2,at=c(seq(0,3,0.5),3.5),las=2,labels=c(seq(0,3,0.5),">3"), lwd =1,lty = 1)
    
    par(mar=c(0,5,0,1))
    plot(seq(0,1,0.5),seq(0,1,0.5),type = "n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    arrows(0, 0.85, 1, 0.85, col = 1,lwd=5,code = 2,length=0.2)
    text(0.05,0.55,"80myr",cex=1)
    text(0.95,0.55,"5myr",cex=1)
    
    
    
    
  }  
  
  #version7 fig2.eAlu_tAlu.age_oddRatio.png fig2.eAlu_tAlu.age_pvalue.png
  { 
    
    subAlu=c("AluJo","AluJb","AluSx","AluSq","AluSg","AluSc","AluSp","AluY","AluYb8","AluYa5")
    
    rfile=eoutp #AluJb   117273  386     0.097995        0.159702        0
    e.odds=c()
    for(i in 1:length(subAlu)){
      flag=which(rfile[,1]==subAlu[i])
      if(length(flag)!=0){
        e.odds=c(e.odds,rfile[flag,6])
      }else{
        e.odds=c(e.odds,0)
      }
    }
    e.odds=as.numeric(e.odds)
    
    rfile=toutp#AluJb   117273  386     0.097995        0.159702        0
    p.odds=c()
    for(i in 1:length(subAlu)){
      flag=which(rfile[,1]==subAlu[i])
      if(length(flag)!=0){
        p.odds=c(p.odds,rfile[flag,6])
      }else{
        p.odds=c(p.odds,0)
      }
    }
    p.odds=as.numeric(p.odds)
    
    tp=-log10(tpvalue)
    tp[c(1,3)]=3.5
    ep=-log10(epvalue)
    t.pch=rep(1,length(tp))
    e.pch=rep(1,length(ep))
    t.pch[which(tp>(-log10(0.01)))]=19
    e.pch[which(ep>(-log10(0.01)))]=19
    
    #fig2.eAlu_tAlu.age_oddRatio.png
    layout(matrix(c(rep(1,12),2,2),ncol = 2,byrow = T))
    par(mar=c(4,3,1,1),cex=1)
    plot(1:length(subAlu),e.odds,lwd=2,type="o",pch=e.pch,col="tan1",xaxt="n",yaxt="n",xlab = "",ylim = c(0.3,1.5),main="",ylab = "Odds ratio")
    par(new=T)
    plot(1:length(subAlu),p.odds,type="o",lwd=2,pch=t.pch,col="skyblue1",xaxt="n",yaxt="n",xlab = "",ylim = c(0.3,1.5),ylab="")
    axis(2,at=seq(0,1.5,0.5),las=2,labels=seq(0,1.5,0.5), lwd =1,lty = 1)
    legend("topright",inset=0.01,legend = c("Alu on Enhancer","Alu on Tpromoter"),col = c("tan1","skyblue1"),lty = c(1, 1),cex=0.9,pch=21,bty="n")
    axis(1,at=1:length(subAlu),labels=subAlu, las=2,lty = 1, lwd = 1)
    
    par(mar=c(0,3,0,1))
    plot(seq(0,1,0.5),seq(0,1,0.5),type = "n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    arrows(0, 0.85, 1, 0.85, col = 1,lwd=5,code = 2,length=0.2)
    text(0.05,0.55,"80myr",cex=1)
    text(0.95,0.55,"5myr",cex=1)
    
    #fig2.eAlu_tAlu.age_pvalue.png
    library(plotrix)
    layout(matrix(c(rep(1,12),2,2),ncol = 2,byrow = T))
    par(mar=c(4,3,1,1),cex=1)
    tp=-log10(tpvalue)
    tp[c(1,3)]=3.5
    plot(1:length(subAlu),-log10(epvalue),lwd=2,pch=e.pch,type="o",col="tan1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,3.6),main="",ylab = "-log10(pvalue)")
    par(new=T)
    plot(1:length(subAlu),tp,type="o",lwd=2,pch=t.pch,col="skyblue1",xaxt="n",yaxt="n",xlab = "",ylim = c(0,3.6),ylab="")
    axis.break(2, 3.4)
    legend("topright",inset=0.01,legend = c("Alu on Enhancer","Alu on Tpromoter"),col = c("tan1","skyblue1"),lty = c(1, 1), bg = "gray90",cex=0.9,pch=21,bty="n")
    axis(1,at=1:length(subAlu),labels=subAlu, las=2,lty = 1, lwd = 1)
    axis(2,at=c(seq(0,3,0.5),3.5),las=2,labels=c(seq(0,3,0.5),">3"), lwd =1,lty = 1)
    
    par(mar=c(0,3,0,1))
    plot(seq(0,1,0.5),seq(0,1,0.5),type = "n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    arrows(0, 0.85, 1, 0.85, col = 1,lwd=5,code = 2,length=0.2)
    text(0.05,0.55,"80myr",cex=1)
    text(0.95,0.55,"5myr",cex=1)
    
    
    
    
  }  
  
  
  
}



data<-data.frame(Sample<-c(rep('control1',3),rep('control2',3),rep('control3',3),rep('treat1',3),rep('treat2',3),rep('treat3',3),rep('treat4',3)), 
                 contion<-rep(c('Cell','Tissue','Organ'),7), 
                 value<-c(503,264,148,299,268,98,363,289,208,108,424,353,1,495,168,152,367,146,48,596,143))

colnames(data)=c('sample',"contion","value")

ggplot(data,mapping = aes(Sample,value,fill=contion))+
  geom_bar(stat='identity',position='fill') +
  labs(x = 'Sample',y = 'frequnency') +
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





xfile=as.matrix(read.table(file="~/bxDownload/panTro5/Hi-C/GSM3685694_chimpanzee_panTro5_mapping_domin.list_50k.txt",header = T))
yfile=as.matrix(read.table(file="/data/baix/bxDownload/panTro5/TAD/chimp.TAD_liftOverFromGM12878TAD.txt",header = F))

len1=as.numeric(xfile[,3])-as.numeric(xfile[,2])
len2=as.numeric(yfile[,3])-as.numeric(yfile[,2])

par(mfrow=c(2,1),cex=1.2)
hist(log10(len1),xlab="log10(TAD length)",main="GSM3685694 TAD")
hist(log10(len2),xlab="log10(TAD length)",main="TAD from GM12878 Liftover ")


#2021.01.22
#caclulate panTro5 chromatin length  ->chrom_panTro5.sizes
{
  rfile=file("/data/baix/bxDownload/panTro5/refSeq/panTro5.fa","r")
  #rfile=file("~/bxDownload/panTro5/Hi-C/test/a.fa","r")
  oneline=readLines(rfile,n=100000)
  
  chr=c()
  len=c()
  while(length(oneline)!=0){
    firstLetter=substr(oneline,1,1)
    flag=which(firstLetter==">")
    if(length(flag)!=0){
      if(1%in%flag){
        seqs=cbind(flag+1,c(flag[-1]+1),length(oneline))
        for(i in 1:nrow(seqs)){
          chr=c(chr,strsplit(oneline[flag[i]],">")[[1]][2])
          leni=sum(nchar(oneline[seqs[i,1]:seqs[i,2]]))
          len=c(len,leni)
        }
      }else{
        seqs=cbind(c(1,flag+1),c(flag-1,length(oneline)))
        
        for(i in 1:nrow(seqs)){
          if(i==1){
            leni=leni+sum(nchar(oneline[seqs[i,1]:seqs[i,2]]))
            len[length(len)]=leni
          }else{
            chr=c(chr,strsplit(oneline[flag[i-1]],">")[[1]][2])
            leni=sum(nchar(oneline[seqs[i,1]:seqs[i,2]]))
            len=c(len,leni)
          }
          
        }
        
      }
    }else{
      seqs=c(1,length(oneline))
      leni=leni+sum(nchar(oneline[seqs[1]:seqs[2]]))
      len[length(len)]=leni
    }
    
    oneline=readLines(rfile,n=100000)
    
  }
  
  close(rfile)
  out=cbind(chr,len)
  
  chr_n=paste("chr",c(1,"2A","2B",3:22,"X","Y"),sep="")
  
  flag=which(chr%in%chr_n)
  
  
  write.table(out[flag,],file ="~/bxDownload/panTro5/Hi-C/chrom_panTro5.sizes",quote =F,sep = "\t",row.names = F,col.names = F )
  
  
}




#different resolution non-0 percentage
{
  xfile=as.matrix(read.table(file="/data/baix/bxDownload/hg19/refSeq/ChromSize.csv",header = F))
  #chr1    249250621
  yfile=as.matrix(read.table(file="/data/baix/bxDownload/hg19/GM12878/Hi-C/juicer_dump.out/bin.number.txt",header=F))
  #48283373        chr10
  resolutions=5000 
  percent=c()
  bins_num=c()
  interact=c()
  for(i in 1:nrow(xfile)){
    
    yflag=which(yfile[,2]==xfile[i,1])
    if(length(yflag)!=0){
      bins=as.numeric(xfile[i,2])/resolutions
      interactI= bins*(bins-1)/2+bins
      percentI=as.numeric(yfile[yflag,1])/interactI
      percent=rbind(percent, c(xfile[i,1],round(percentI,2)))
      bins_num=c(bins_num,as.numeric(yfile[yflag,1]))
      interact=c(interact,interactI)
    }
    
  }
  
  sum(bins_num)/sum(interact)
  

  write.table(percent,file ="non-0.percent_5kb.txt",quote =F,sep = "\t",row.names = F,col.names = F )
  
}

#identify eRNA  
{
  setwd("~/newCell4Type/1geteRNAwithCAGE/")
  cells=c(  "HepG2", "MCF-7", "K562","HeLa-S3")
  for(k in 1:2){
    k=4
    afile=as.matrix(read.table(file = paste("/data/baix/bxDownload/hg19/",cells[k],"/enhancerPeak/Dnase_me1_k27acPeak.txt",sep = ""),header = F))
    #chr10   100019755       100019905       chr10   100018601       100021066       chr10   100019690   100020862
    sortSE=t(apply(afile[,c(5:6,8:9)], 1, function(x){sort(as.numeric(x))}))
    
    out=cbind(afile[,1],sortSE[,2:3],sortSE[,c(1,4)])
    write.table(out,file =paste("~/newCell4Type/1geteRNAwithCAGE/",cells[k],"/1me1_k27ac.union_intersect.txt",sep = ""),quote =F,sep = "\t",row.names = F,col.names = F )
    
  }
  
  
  bfile=as.matrix(read.table(file = "./MCF-7/3me1_k27ac.union_intersect_withDHS_CAGE.txt",header = F))
  #chr10 100150001 100152061 100150474 100151512 chr10 100150896 100151026 chr10:100150896:100151026:-:0.65:0.812895 0 - 18.1383 3.93e-08 0
  #chr10   1033327 1033812 1029046 1035763 chr10   1033666 1034032 TSS_chr10_minus_1017096_1034408_pk1 1000     -       21.0
  
  flag=which(bfile[,11]=="+")
  out=rbind(bfile[flag,c(1,7,5,11,1:8,12)],bfile[-flag,c(1,4,8,11,1:8,12)])
  elen=as.numeric(out[,3])-as.numeric(out[,2])
  hist(elen)
  a=which(elen<0)
  out=unique(out)
  
  write.table(out,file ="./MCF-7/3sub.me1_k27ac.intersect_union.With_DHS.CAGE.txt",quote =F,sep = "\t",row.names = F,col.names = F )
  
  getCAGEeRNA=function(input,output){
    afile=as.matrix(read.table(file = input,header = F))
    #chr10 100164643 100165459 + chr10 100162717 100165706 100162892 100165459 chr10 100164643 100164759 0.2844
    cageIn1=as.numeric(afile[,2])-as.numeric(afile[,8])
    cageIn2=as.numeric(afile[,9])-as.numeric(afile[,3])
    flag=which(cageIn1>=0 &cageIn2>=0)
    cageIn.File=afile[flag,]
    paste.file=apply(cageIn.File[,c(1,4,8,9)],1,function(x){paste(x,collapse = ",")})
    Upaste.file=unique(paste.file)
    Upaste.fileL=length(Upaste.file)
    outflag=c()
    out=c()
    for(i in 1:Upaste.fileL){
      flagi=which(paste.file==Upaste.file[i])
      flagiL=length(flagi)
      if(flagiL==1){
        outflag=c(outflag,flagi)
      }else{
        ilen=as.numeric(cageIn.File[flagi,3])-as.numeric(cageIn.File[flagi,2])
        ilen.l=which(ilen<200)
        ilen.T=which(ilen>=200 & ilen<=2000)
        ilen.F=which(ilen>2000)
        
        if(length(c(ilen.l,ilen.T))==length(ilen)){
          flagii=which.max(as.numeric(cageIn.File[flagi[ilen.T],13]))
          outflag=c(outflag,flagi[ilen.T[flagii]])
        }else if(length(ilen.F)!=0){
          if(cageIn.File[flagi[1],4]=="+"){
            sorts=sort(as.numeric(cageIn.File[flagi[ilen.F],2]))
            if(length(ilen.T)!=0){
              flagii=which.max(as.numeric(cageIn.File[flagi[ilen.T],13]))
              outflag=c(outflag,flagi[ilen.T[flagii]])
              sorts=c(sorts,as.numeric(cageIn.File[flagi[ilen.T[flagii]],2]))
            }else{
              sorts=c(sorts,as.numeric(cageIn.File[flagi[1],3]))
            }
            num=1
            k=1
            while (num+k<=length(sorts)) {
              if((sorts[num+k]-1-sorts[num])>=200){
                if((sorts[num+k]-1-sorts[num])<=2000){
                  out=rbind(out,c(cageIn.File[flagi[1],1],sorts[num],sorts[num+k]-1,cageIn.File[flagi[1],c(4,8:9)],cageIn.File[flagi[ilen.F[num]],13]))
                  
                }
                num=num+k
                k=1
                
              }else{
                k=k+1
              }
            }
            
          }else{
            sorts=sort(as.numeric(cageIn.File[flagi[ilen.F],3]))
            if(length(ilen.T)!=0){
              flagii=which.max(as.numeric(cageIn.File[flagi[ilen.T],13]))
              outflag=c(outflag,flagi[ilen.T[flagii]])
              sorts=c(sorts,as.numeric(cageIn.File[flagi[ilen.T[flagii]],3]))
            }else{
              sorts=c(as.numeric(cageIn.File[flagi[1],2]),sorts)
            }
            num=length(sorts)
            k=1
            while ((num-k)>=1) {
              if((sorts[num]-sorts[num-k]-1)>=200){
                if((sorts[num]-sorts[num-k]-1)<=2000){
                  out=rbind(out,c(cageIn.File[flagi[1],1],sorts[num-k]+1,sorts[num],cageIn.File[flagi[1],c(4,8:9)],cageIn.File[flagi[ilen.F[num-1]],13]))
                  
                }
                num=num-k
                k=1
              }else{
                k=k+1
              }
            }
            
          }
          
        }
      }
      
      
    }
    
    Out=rbind(cageIn.File[outflag,c(1:4,8:9,13)],out)
    write.table(Out,file =output ,sep="\t",quote = F,row.names = F,col.names = F)
    
    hist(log(as.numeric(Out[,7])),breaks = 100)
    elen=as.numeric(Out[,3])-as.numeric(Out[,2])
    hist(elen)
    length(which(elen>2000))
  }
  
  input="./HeLa-S3/4me1_k27ac.union_intersect_withDHS_CAGE_movedTSS2kb.txt"
  output="./HeLa-S3/5enhancer_intersect.me1_k27ac.bed"
  getCAGEeRNA(input,output)
  
  
  setwd("~/newCell4Type/1geteRNAwithCAGE/HeLa-S3/")
  eRNAname=function(input,output,locat){
    afile=as.matrix(read.table(file = input,header = F))
    #chr1 101756744 101760265 + 101755845 101760265 0.2904
    
    aa=apply(afile[,c(1,5:6)], 1,function(x){paste(x,collapse = ",")})
    Uaa=unique(aa)
    name=rep(0,length(aa))
    for(i in 1:length(Uaa)){
      flagp=which(aa==Uaa[i] & afile[,4]=="+")
      flagm=which(aa==Uaa[i] & afile[,4]=="-")
      if(length(flagp)==1){
        name[flagp]=paste(locat,"+",i,sep="")
      }else{
        name[flagp]=paste(locat,"+",i,".",1:length(flagp),sep="")
      }
      if(length(flagm)==1){
        name[flagm]=paste(locat,"-",i,sep="")
      }else{
        name[flagm]=paste(locat,"-",i,".",1:length(flagm),sep="")
      }
      
      
    }
    
    Out=cbind(afile[,1:3],name,substr(name,1,9),afile[,c(4,7)])
    elen=as.numeric(Out[,3])-as.numeric(Out[,2])
    flag=which(elen>=200 & elen<=2000)
    write.table(Out[flag,],file =output ,sep="\t",quote = F,row.names = F,col.names = F)
    write.table(Out[-flag,],file =paste("super.",output,sep="") ,sep="\t",quote = F,row.names = F,col.names = F)
    
    
    
  }
  input="6Intergenic.enhancer_intersect.me1_k27ac.bed"
  output="7Extragenic.eRNA.bed"
  eRNAname(input,output,"GMExeRNA")
  
  input="6Intragenic.enhancer_intersect.me1_k27ac.bed"
  output="7Intragenic.eRNA.bed"
  eRNAname(input,output,"GMIneRNA")
  
  
  
  
  
  
}    


#supplymentary Figure-> expression
{
  xfileName=c("E+P+Alu","E+P-Alu","E-P+Alu","E-P-Alu","E-uniq.P+Alu")
  cells=c("HepG2", "MCF-7", "K562","HeLa-S3")
  #get E-P wether have Alu
  {
    setwd("~/newCell4Type/2getEP/")
   tfile=as.matrix(read.table(file="/data/baix/bxDownload/hg19/refSeq/gencode/Protein.promoter_withAlu.bed",header = F))
   targetwithAlu=unique(tfile[,4])
   for(k in 1:4){
 
     k=4
     etfile=as.matrix(read.table(paste(cells[k],"/eRNA.Trotein.txt",sep=""),header = F))
     #chr10   100162892       100163920       GMIneRNA-386    GMIneRNA-       -       chr10   100009780       100012780       ENSG00000230928.1       RP11-34A14.3    +       antisense    
     efile=as.matrix(read.table(paste(cells[k],"/eRNA.Alu.txt",sep=""),header = F))
     #chr10 105227199 105227486 GMExeRNA9- GMExeRNA9 - chr10 105227020 105227306 AluSx1 SINE/Alu +
     
     eRNAwithAlu=unique(efile[,4])
     {
      eRNA.et=apply(etfile[,1:4],1,function(x){paste(x,collapse = "-")})
       eRNA=unique(apply(efile[,c(1:3,6)],1,function(x){paste(x,collapse = "-")}))
       eflag=eRNA.et%in%eRNA
       tflag=etfile[,8]%in%targetwithAlu
       both=which(eflag & tflag)
       ehave=which(eflag & !tflag)
       phave=which(!eflag & tflag)
       noneHave=which(!eflag & !tflag)
       unq.phave=which(etfile[phave,8]%in%etfile[both,8])
       
       write.table(etfile[both,],file=paste(cells[k],"/",xfileName[1],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
       write.table(etfile[ehave,],file=paste(cells[k],"/",xfileName[2],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
       write.table(etfile[phave,],file=paste(cells[k],"/",xfileName[3],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
       write.table(etfile[noneHave,],file=paste(cells[k],"/",xfileName[4],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
       write.table(etfile[phave[-unq.phave],],file=paste(cells[k],"/",xfileName[5],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
       
     }
     eflag=etfile[,4]%in%eRNAwithAlu
     tflag=etfile[,10]%in%targetwithAlu
     both=which(eflag & tflag)
     ehave=which(eflag & !tflag)
     phave=which(!eflag & tflag)
     noneHave=which(!eflag & !tflag)
     unq.phave=which(etfile[phave,10]%in%etfile[both,10])
     
     write.table(etfile[both,],file=paste(cells[k],"/",xfileName[1],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
     write.table(etfile[ehave,],file=paste(cells[k],"/",xfileName[2],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
     write.table(etfile[phave,],file=paste(cells[k],"/",xfileName[3],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
     write.table(etfile[noneHave,],file=paste(cells[k],"/",xfileName[4],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
     write.table(etfile[phave[-unq.phave],],file=paste(cells[k],"/",xfileName[5],".txt",sep = ""),quote=F,sep="\t",row.name=FALSE,col.name=FALSE)
     
   }
       
  }
  
  
  {
    setwd("~/newCell4Type/2getEP/")
    xfileName=c("E+P+Alu","E+P-Alu","E-P+Alu","E-P-Alu","E-uniq.P+Alu")
    cells=c("HepG2", "MCF-7", "K562","HeLa-S3")
    c(0.5,0.75,0.4)
    ck=4
    
    yfile=as.matrix(read.table(file=paste("/data/baix/bxDownload/hg19/",cells[ck],"/RNA-seq/rep1/1gencode.v19.gene.TPM-FPKM.bed",sep=""),header=T))
   #chr1        11869           14412       +       ENSG00000223972.4       DDX11L1      1.07            1.67
    FPKMvalue=list()
    for(g in 1:5){
      xfile=as.matrix(read.table(file=(paste(cells[ck],"/maybe/",xfileName[g],".txt",sep="")),header=FALSE))
      #chr10 101742138 101744233 GMIneRNA10- GMIneRNA10 -  chr10 101684767 101687767 ENSG00000227695.1 DNMBP-AS1 + antisense
      gene=unique(xfile[,8])
      geneL=length(gene)
      flagy=which(yfile[,5]%in%gene)
      FPKMvalue[[g]]=as.numeric(yfile[flagy,8])
    }
    both=FPKMvalue[[1]]
    ehave=FPKMvalue[[2]]
    phave=FPKMvalue[[3]]
    nonhave=FPKMvalue[[4]]
    phave_moved=FPKMvalue[[5]]
    
    seqs=c(seq(0,50,5),30000)
    #seqs=c(seq(0,1,0.1),5,seq(10,100,10),20000)
    #par(mfrow=c(2,2))
    a1=hist(both,breaks = seqs,plot=F)
    a2=hist(ehave,breaks = seqs,plot=F)
    a3=hist(phave,breaks = seqs,plot=F)
    a4=hist(nonhave,breaks =seqs,plot=F)
    a5=hist(phave_moved,breaks =seqs,plot=F)
    a=cbind(a1$counts/length(both),a2$counts/length(ehave),a3$counts/length(phave),a4$counts/length(nonhave),a5$counts/length(phave_moved))
    b=c()
    for(k in 1:5){
      b=cbind(b,cumsum(a[,k]))
    }
    ymax=1
    ymin=0.4
    
    pdf(file=paste("~/eRNAprojectPlots/supplementaryPlot/",cells[ck],"geneExpression.pdf",sep=""),family = "ArialMT")
    
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    par(mar=c(4,4,4,1),cex=1.1,mgp=c(1.8,0.5,0))
    plot(1:dim(a)[1],b[,1],xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col="orangered",ylim=c(ymin,ymax),lwd=2,xlab = "Expression level of target genes(FPKM)",ylab = "Cumulative frequency",main = cells[ck])
    par(new=T)
    plot(1:dim(a)[1],b[,3],xaxt="n",yaxt="n",type="l",col="steelblue2",ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    par(new=T)
    plot(1:dim(a)[1],b[,2],xaxt="n",yaxt="n",type="l",col="goldenrod1",ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    par(new=T)
    plot(1:dim(a)[1],b[,4],xaxt="n",yaxt="n",type="l",col="darkseagreen3",ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    par(new=T)
    plot(1:dim(a)[1],b[,5],xaxt="n",yaxt="n",type="l",col="steelblue2",ylim=c(ymin,ymax),lty=2,lwd=2,xlab="",ylab ="",main="")
    
    legend("bottomright",inset=0.05,legend = c(paste("E+P+, n=",length(both),sep=""),paste("E-P+, n=",length(phave),sep=""),paste("E+P-, n=",length(ehave),sep=""),paste("E-P-, n=",length(nonhave),sep=""),paste("E-P+, n=",length(phave_moved),sep="")),col =cols[c(1:4,2)],lty = c(1, 1, 1,1,2),cex=1,lwd = 2,bty = "n")
    #axis(1,at=1:dim(a)[1],labels=seqs[-1], col = "black", lty = 1, lwd = 0.5)
    axis(1,at=seq(0,dim(a)[1],2),labels=seq(0,50,10), col = "black", lty = 1, lwd = 1)
    
    dev.off()
    
  }
  
  
  
  
}



a=rbind(c(97,819),c(35,3688))
fisher.test(a)

b=rbind(c(1208,9115),c(383,6086))
fisher.test(b)

en=rbind(c(75,757),c(30,532))
en=rbind(c(4,29),c(30,532))
fisher.test(en)

fisher.test(156/841)/(1903/11477)
(124/604)/(1903/11477)

a=rbind(c(134,675),c(35,3688))
chisq.test(a)
fisher.test(a)

b=rbind(c(1494,9043),c(383,6086))
chisq.test(b)
fisher.test(b)

en=rbind(c(108,503),c(30,532))
chisq.test(en)
fisher.test(en)

fisher.test(156/841)/(1903/11477)
(124/604)/(1903/11477)


# get rondom EPI in same TAD and same distance-> then used to check whether supported by ChIA-PET loop
{
  setwd("~/GM12878/5newWithHic/")
  
  #yfile=as.matrix(read.table(file = "EP_listBy_up.down_TAD.txt",header = F))
  yfile=as.matrix(read.table(file = "eRNA.targetProtein_TAD.txt",header = F))
  #chr1    916497  919497       ENSG00000187642.5       C1orf170        -      
  #chr1    974026  974300  GMIneRNA+1      GMIneRNA+   + protein_coding  chr1    915000  1005000 TAD1
  
  EP.names=apply(yfile[,c(4,10)],1,function(x){paste(x,collapse ="-")})
  EP.names.u=unique(EP.names)
  outflag=c()
  for(i in 1:length(EP.names.u)){
    yflag=which(EP.names==EP.names.u[i])
    outflag=c(outflag,yflag[1])

  }
  write.table(yfile[outflag,],file ="eRNA.targetProtein_TAD.txt",quote =F,sep = "\t",row.names = F,col.names = F )

  
  yfile=as.matrix(read.table(file = "EP_listBy_up.down_TAD_unique.txt",header = F))
  #chr1    916497  919497       ENSG00000187642.5       C1orf170        -      
  #chr1    974026  974300  GMIneRNA+1      GMIneRNA+   + protein_coding  chr1    915000  1005000 TAD1
  yfileL=nrow(yfile)
  down.e_up.s=as.numeric(yfile[,9])-as.numeric(yfile[,2])
  TAD.s=as.numeric(yfile[,15])
  TAD.e=as.numeric(yfile[,16])-down.e_up.s
  up.e_s=as.numeric(xfile[,3])-as.numeric(xfile[,2])
  down.s_up.s=as.numeric(xfile[,8])-as.numeric(xfile[,2])
  reps=1000
  up.start=c()
  for(i in 17165:yfileL){
    up.starti=sample(TAD.s[i]:TAD.e[i],reps,replace = T)
    up.start=rbind(up.start,cbind(1:reps,up.starti))
  }
  for(k in 1:reps){
    repi=yfile[,1:13]
    up.starti=which(up.start[,1]==k)
    repi[,2]=up.start[up.starti,2]
    repi[,3]=up.start[up.starti,2]+up.e_s
    repi[,8]=up.start[up.starti,2]+down.s_up.s
    repi[,9]=up.start[up.starti,2]+down.e_up.s
    write.table(repi,file =paste("./RandomEPInTAD/Random_EP_listBy_up.downInTAD",k,".txt",sep=""),quote =F,sep = "\t",row.names = F,col.names = F )
    
  }
 
  
  #random get enhancer region
  {
    xfile=as.matrix(read.table(file = "eRNA.targetProtein_TAD.txt",header = F))
    #chr1       974026          974300       GMIneRNA+1      GMIneRNA+       +       chr1       916497          919497       ENSG00000187642.5 C1orf170        -       protein_coding  
    #chr1       915000         1005000       TAD1
    xfileL=nrow(xfile)
    tad=unique(xfile[,17])
    up.e_s=as.numeric(xfile[,3])-as.numeric(xfile[,2])
    TAD.s=as.numeric(xfile[,15])
    TAD.e=as.numeric(xfile[,16])-up.e_s
    reps=1000
    up.start=c()
    for(i in 1:xfileL){
      up.starti=sample(TAD.s[i]:TAD.e[i],reps,replace = T)
      up.start=rbind(up.start,cbind(1:reps,up.starti))
    }
    for(k in 1:reps){
      repi=xfile[,1:13]
      up.starti=which(up.start[,1]==k)
      repi[,2]=up.start[up.starti,2]
      repi[,3]=up.start[up.starti,2]+up.e_s

      write.table(repi,file =paste("./RandomEnhancerInTAD/Random_EP_listBy_up.downInTAD",k,".txt",sep=""),quote =F,sep = "\t",row.names = F,col.names = F )
      
    }
    
    
    
    
    
    
    
    
    
  }
  
  
  setwd("~/GM12878/5newWithHic/")
  xfile=as.matrix(read.table(file ="~/GM12878/5newWithHic/EP_listBy_up.down_CHIApet.txt",header = F))
  
  yfile=as.matrix(read.table(file ="~/GM12878/5newWithHic/RandomEnhancerInTAD/bnum.txt",header = F))
  
  num=as.numeric(yfile[,1])
  
  zscore=round((nrow(xfile)-mean(num))/sd(num),2)#83.37404
  p=1-pnorm(5)
  
  format(1-pnorm(83.37404),scientific = T,digits = 3)
  
  pdf(file = "~/eRNAprojectPlots/supplementaryPlot/SuppotByChIAPetLoop",family = "ArialMT")
  par(mfrow=c(1,1),mar = c(4, 4, 3,1),mgp=c(1.8,0.5,0),cex=1)
  hist(num,breaks = 20,col = "skyblue1",xlim=c(550,850),xaxt="n",cex.lab=1.2,cex.axis=1,xlab="#EPI",main="Support by ChIA-PET loop")
  axis(1,at=seq(550,800,50),las=1,labels=c(seq(550,700,50),2500,2550), lty = 1, lwd = 1)
  axis.break(1,730)
  lines(c(nrow(xfile)-1750,nrow(xfile)-1750),c(-10,60),col=2,lwd=2)
  text(780,70,paste("Z=",zscore,sep=""))
  dev.off()  
  
  
  
}#192.168.77.39

name.out=cbind(paste("Random_EP_listBy_up.downInTAD",1:1000,".txt",sep=""),paste("2Random",1:1000,".txt",sep=""),paste("Random_EP_listBy_up.downInTAD",1:1000,"_ChIApet.txt",sep=""))
write.table(name.out,file ="namefile.txt",quote =F,sep = "\t",row.names = F,col.names = F )



setwd("/data/baix/bxDownload/hg19/GWAS_catalog/v1.0.3/")
xfile=as.matrix(read.table(file="gwas-catalog_chr_pos_rsname_maybeCorrelateDisease.txt",header = F))
yfile=as.matrix(read.table(file="trait_selected.txt",header = F))

flag=which(xfile[,6]%in%yfile)
out=xfile[-flag,]
write.table(out,file ="gwas-catalog_chr_pos_rsname_DiseaseAssociated.txt",quote =F,sep = "\t",row.names = F,col.names = F )




{
  
  
  setwd("/data/baix/bxDownload/homologene/")
  hfile=as.matrix(read.table(file = "./promoterWithAlu/human.orthologyPromoter_overlap0.9_Alu.txt",header = F))
  cfile=as.matrix(read.table(file = "./promoterWithAlu/chimp.orthologyPromoter_overlap0.9_alu.txt",header = F))
  
  #chr10 102887257 102890257 ENSG00000107807.8 TLX1 + protein_coding 3 chr10 102887913 102888068 FRAM SINE/Alu -
  #chr10	100751789	100754789	ENSPTRG00000002828	CNNM1	+	protein_coding	12085	chr10	100751656	100751927	AluJr	SINE	-
  h.midAlu=(as.numeric(hfile[,10])+as.numeric(hfile[,11]))/2
  c.midAlu=(as.numeric(cfile[,10])+as.numeric(cfile[,11]))/2
  h.aluP=apply(cbind(hfile[,c(2:3,6)],h.midAlu),1,function(x){if(x[3]=="+"){as.numeric(x[4])-as.numeric(x[1])}else{-(as.numeric(x[4])-as.numeric(x[2]))}})
  c.aluP=apply(cbind(cfile[,c(2:3,6)],c.midAlu),1,function(x){if(x[3]=="+"){as.numeric(x[4])-as.numeric(x[1])}else{-(as.numeric(x[4])-as.numeric(x[2]))}})
  
  {
    
    num=unique(c(hfile[,8],cfile[,8]))
    numL=length(num)
    hx.cy=c()
    for(i in 1:10){
      hflag=which(hfile[,8]==num[i])
      cflag=which(cfile[,8]==num[i])
      # hflag
      # cflag
      if(length(hflag)!=0 && length(cflag)!=0){
        hc.aluP=data.frame(c(h.aluP[hflag],c.aluP[cflag]),c(rep(1,length(hflag)),rep(2,length(cflag))),c(hfile[hflag,12],cfile[cflag,12]),c(hfile[hflag,4],cfile[cflag,4]),fix.empty.names  =F,stringsAsFactors = F)
        
        orders=order(hc.aluP[,1]) 
        hc.orders=hc.aluP[orders,]
        hflag1=which(hc.orders[,2]==1)
        cflag1=which(hc.orders[,2]==2)
        useflag=c()
        hc.dist=c()
        for(m in 1:length(hflag1)){
          if((hflag1[m]-1)>0&& !(hc.orders[hflag1[m]-1,2]%in%hc.orders[hflag1,2])){
            disUp=abs(hc.orders[hflag1[m]-1,1]-hc.orders[hflag1[m],1])
            hc.dist=rbind(hc.dist,c(hflag1[m],hflag1[m]-1,disUp))
          }
          if((hflag1[m]+1)<=length(orders)&& !(hc.orders[hflag1[m]+1,2]%in%hc.orders[hflag1[m],2])){
            disDown=abs(hc.orders[hflag1[m]+1,1]-hc.orders[hflag1[m],1])
            hc.dist=rbind(hc.dist,c(hflag1[m],hflag1[m]+1,disDown))
          }
        }
        hc.dist=matrix(hc.dist,ncol = 3)
        while(nrow(hc.dist)!=0){
          hc.flag=which.min(hc.dist[,3])
          hx.cy=rbind(hx.cy,c(hc.orders[hc.dist[hc.flag,1],c(1,3:4)],hc.orders[hc.dist[hc.flag,2],c(1,3:4)]))
          useflag=c(useflag,hc.dist[hc.flag,1:2])
          Left.hc.flag=which(hc.dist[,1]!=hc.dist[hc.flag,1]&hc.dist[,2]!=hc.dist[hc.flag,2])
          hc.dist=matrix(hc.dist[Left.hc.flag,],ncol=3)
          
        }
        
        hflag1.sect=setdiff(hflag1,useflag)
        cflag1.sect=setdiff(cflag1,useflag)
        if(length(hflag1.sect)!=0){
          for(k in 1:length(hflag1.sect)){
            hx.cy=rbind(hx.cy,c(hc.orders[hflag1.sect[k],c(1,3:4)],0,0,hc.orders[cflag1[1],4]))
          }
        }
        if(length(cflag1.sect)!=0){
          for(k in 1:length(cflag1.sect)){
            hx.cy=rbind(hx.cy,c(0,0,hc.orders[hflag1[1],4],hc.orders[cflag1.sect[k],c(1,3:4)]))
          }
        }
        
      }else if(length(hflag)==0 && length(cflag)!=0){
        for(k in 1:length(cflag)){
          hx.cy=rbind(hx.cy,c(0,0,0,c.aluP[cflag[k]],cfile[cflag[k],c(12,4)]))
        }
      }else{
        for(k in 1:length(hflag)){
          hx.cy=rbind(hx.cy,c(h.aluP[hflag[k]],hfile[hflag[k],c(12,4)],0,0,0))
        }
        
      }
      
      
    }
    
  }

}

#2021.04.28 ATAC-seq and FPKM
{
  setwd("~/GM12878/15ATAC-seq/")
  xfileName=c("GM12878eRNAandProteinTargetBothHaveAlu","GM12878eRNAHaveAluaveAndProteinTargetNoAlu","GM12878eRNAnoAluandProteinTargetHaveAlu","GM12878eRNAandProteinTargetNoneHaveAlu")
  efile=as.matrix(read.table(file="~/GM12878/2newGetEPpairsInTAD/CAGEhead/eRNA.protein.bed",header=F))
  tfile=as.matrix(read.table(file="~/GM12878/2newGetEPpairsInTAD/CAGEhead/proteinTarget.bed",header=F))
  rfile=file("eRNA_foldchange.bigwig.txt","r")
  ewfile=readLines(rfile)
  close(rfile)
  twfile=as.matrix(read.table(file="proteinTarget_foldchange.bigwig.txt",header=F))
  
  
  e.bw=list()
  t.bw=list()
  target=list()
  for(g in 1:4){
    xfile=as.matrix(read.table(file=(paste("~/eRNAprojectPlots/newPlot/partResultData/",xfileName[g],".txt",sep="")),header=FALSE))
    #chr10 101742138 101744233 GMIneRNA10- GMIneRNA10 -  chr10 101684767 101687767 ENSG00000227695.1 DNMBP-AS1 + antisense
    #gene=apply(xfile,1,function(x){strsplit(x[8],"~")[[1]][1]})
    eRNA=unique(xfile[,4])
    gene=unique(xfile[,10])
    e.flag=which(efile[,4]%in%eRNA)
    t.flag=which(tfile[,4]%in%gene)
    e.bw.g=apply(cbind(ewfile[e.flag],1), 1, function(x){median(as.numeric(strsplit(x[1],"\t")[[1]]))})
    t.bw.g=apply(twfile[t.flag,],1,mean)
    
    e.bw[[g]]=e.bw.g
    t.bw[[g]]=t.bw.g
    target[[g]]=tfile[t.flag,4]
  }
  
  wilcox.test(e.bw[[1]],e.bw[[4]],alternative = "greater")
  p=wilcox.test(log2(t.bw[[1]]),log2(t.bw[[4]]),alternative = "greater")$p.value
  
  #ATAC-seq_log.signalValue_EP.png
  pdf(file="~/eRNAprojectPlots/JGG.revise/ATAC-seq_log.signalValue_EP.pdf",family = "ArialMT")
  {
    par(mfrow=c(1,1),mar=c(2,3,1,2),mgp=c(1.8,0.5,0))
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    boxplot(log2(e.bw[[1]]),log2(e.bw[[3]]),log2(e.bw[[2]]),log2(e.bw[[4]]), 
            boxwex = 0.2, at = seq(0.6,1.5,0.3),
            col = cols,range=2.5,xaxt="n",
            main = "",
            xlab = "",
            ylab = "Median ATAC-seq signal values",
            xlim = c(0.5, 3), ylim = c(-4, 6))
    boxplot(log2(t.bw[[1]]),log2(t.bw[[3]]),log2(t.bw[[2]]),log2(t.bw[[4]]), add = TRUE,range=2.5,xaxt="n",
            boxwex = 0.2, at =  seq(2,3,0.3), col = cols,yaxt="n")
    axis(1,at=c(1.05,2.45),labels=c("Enhancer","Promoter"), col = 1, lty = 1, lwd = 1,cex.axis=1.2)
    legend("topright",inset=0.05,ncol = 2,legend = c("E+P+","E-P+","E+P-","E-P-"),col =cols[1:4],cex=1,bty = "n",fill = cols)
    
    lines(c(0.6,1.5),c(5,5))
    text(1.05,5.3,"ns")
    
    lines(c(2,2.9),c(3.3,3.3))
    lines(c(2,2),c(3.1,3.3))
    lines(c(2.9,2.9),c(3.1,3.3))
    text(2.45,3.6,format(p,scientific = T,digits = 3))
    
  }
  dev.off()  
  
 gfile=as.matrix(read.table(file=("/data/baix/bxDownload/hg19/GM12878/RNA-seq/rep1/genes.fpkm_tracking"),header=F))
 
  flag1=c()
 for(k in 1:length(t.bw[[4]])){
   a=which(t.bw[[1]]>(t.bw[[4]][k]-0.1) &t.bw[[1]]<(t.bw[[4]][k]+0.1))
   if(length(a)==0){
     a=which(t.bw[[1]]>(t.bw[[4]][k]-0.2) &t.bw[[1]]<(t.bw[[4]][k]+0.2))
   }
   if(length(a)>1){
     flag1=c(flag1,sample(a,1))
   }else if(length(a)==1){
     flag1=c(flag1,a)
   }
   
 }
 

  wilcox.test(t.bw[[1]][flag1],t.bw[[4]])

  FPKM1.flag=which(gfile[,4]%in%target[[1]][flag1])
  FPKM4.flag=which(gfile[,4]%in%target[[4]])
  FPKM1=as.numeric(gfile[FPKM1.flag,10])
  FPKM4=as.numeric(gfile[FPKM4.flag,10])
  p1=wilcox.test(FPKM1,FPKM4,alternative = "greater")$p.value
  
  pdf(file="~/eRNAprojectPlots/JGG.revise/ATAC-seq_log.signalValue_EP_E+P+_E-P-.pdf",family = "ArialMT")
  {
    par(mfrow=c(1,1),mar=c(2,3,1,3),mgp=c(1.8,0.5,0))
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    boxplot(log2(t.bw[[1]][flag1]),log2(t.bw[[4]]),
            boxwex = 0.2, at = seq(0.6,0.9,0.3),
            col = cols[c(1,4)],range=2.5,xaxt="n",
            main = "",
            xlab = "",las=1,
            ylab = "Median ATAC-seq signal values",
            xlim = c(0.5, 1.8), ylim = c(-4.5, 5))
    
    boxplot(log10(FPKM1),log10(FPKM4), add = TRUE,range=1.5,xaxt="n",yaxt="n",
            boxwex = 0.2, at =  seq(1.4,1.7,0.3), col = cols[c(1,4)])
    axis(1,at=c(0.75,1.55),labels=c("ATAC-seq","FPKM"), col = 1, lty = 1, lwd = 1,cex.axis=1.2)
    axis(4,at=seq(-4, 4,2),labels=seq(-4, 4,2), col = 1,las=1, lty = 1, lwd = 1,cex.axis=1.2)
    mtext("Log10 of Target genes' FPKM",side = 4,line =1.8,cex.lab=1)
    
    lines(c(0.6,0.9),c(3,3))
    lines(c(0.6,0.6),c(2.8,3))
    lines(c(0.9,0.9),c(2.8,3))
    text(0.75,3.3,"ns")
    
    lines(c(1.4,1.7),c(4.3,4.3))
    lines(c(1.4,1.4),c(4.1,4.3))
    lines(c(1.7,1.7),c(4.1,4.3))
    text(1.55,4.6,format(p1,scientific = T,digits = 3))
    
  }
  dev.off()
  
  
  
}

#2021.04.29 Hi-C EPI 
{
  setwd("~/GM12878/5newWithHic/")
  xfile=as.matrix(read.table(file="EP_listBy_up.down_up.down.txt",header = F))
  #chr1    916497  919497  ENSG00000187642.5       C1orf170        -       chr1    974026  974300  GMIneRNA+1      GMIneRNA+       + protein_coding   chr1    915000  920000  970000  975000  1.8473549
  
  yfile=as.matrix(read.table(file="~/eRNAprojectPlots/newPlot/partResultData/GM12878eRNAandProteinTargetBothHaveAlu.txt",header = F))
  #chr10   103663844       103664554       GMIneRNA+389    GMIneRNA+       +       chr10   103577696       103580696       ENSG00000198408.9 MGEA5   -       protein_coding
  zfile=as.matrix(read.table(file="~/eRNAprojectPlots/newPlot/partResultData/GM12878eRNAandProteinTargetNoneHaveAlu.txt",header = F))
  #chr10   103663844       103664554       GMIneRNA+389    GMIneRNA+       +       chr10   103577696       103580696       ENSG00000198408.9 MGEA5   -       protein_coding
  
  yfileL=nrow(yfile)
  flag=c()
  Alu.CI=c()
  for(i in 1:yfileL){
    flagi=which((xfile[,4]==yfile[i,4] & xfile[,10]==yfile[i,10])|(xfile[,4]==yfile[i,10] & xfile[,10]==yfile[i,4]))
    flag=c(flag,flagi)
    Alu.CI=c(Alu.CI,mean(as.numeric(xfile[flagi,19])))
  }
  
  zfileL=nrow(zfile)
  flag2=c()
  nonAlu.CI=c()
  for(i in 1:zfileL){
    flagi=which((xfile[,4]==zfile[i,4] & xfile[,10]==zfile[i,10])|(xfile[,4]==zfile[i,10] & xfile[,10]==zfile[i,4]))
    flag2=c(flag2,flagi)
    nonAlu.CI=c(nonAlu.CI,mean(as.numeric(xfile[flagi,19])))
  }
  
  
  Alu.CI=as.numeric(xfile[flag,19])
  
  
  nonAlu.CI=as.numeric(xfile[-flag,19])
  
  t.test(Alu.CI,nonAlu.CI,alternative = "less")
  
  
}

#2021.04.29 Hi-C EPI 
{
  setwd("~/GM12878/5newWithChIA-PET/")
  xfile=as.matrix(read.table(file="eTprotein_petCounts.txt",header = F))
  #chr2    88913448        88914355        GMIneRNA-1551   GMIneRNA-       -       chr2    88989162        88992162        ENSG00000153574.8 RPIA    +       protein_coding  0       0
  
  xfile=as.matrix(read.table(file="~/GM12878/5newWithGridseq/eRNA.targetProtein_ep.distal.Gridseq.txt",header = F))
  yfile=as.matrix(read.table(file="~/eRNAprojectPlots/newPlot/partResultData/GM12878eRNAandProteinTargetBothHaveAlu.txt",header = F))
  #chr10   103663844       103664554       GMIneRNA+389    GMIneRNA+       +       chr10   103577696       103580696       ENSG00000198408.9 MGEA5   -       protein_coding
  zfile=as.matrix(read.table(file="~/eRNAprojectPlots/newPlot/partResultData/GM12878eRNAandProteinTargetNoneHaveAlu.txt",header = F))
  #chr10   103663844       103664554       GMIneRNA+389    GMIneRNA+       +       chr10   103577696       103580696       ENSG00000198408.9 MGEA5   -       protein_coding
  
  yfileL=nrow(yfile)
  flag=c()
  Alu.CI=c()
  for(i in 1:yfileL){
    flagi=which((xfile[,4]==yfile[i,4] & xfile[,10]==yfile[i,10]))
    flag=c(flag,flagi)
    #Alu.CI=c(Alu.CI,mean(as.numeric(xfile[flagi,19])))
  }
  
  zfileL=nrow(zfile)
  flag2=c()
  nonAlu.CI=c()
  for(i in 1:zfileL){
    flagi=which((xfile[,4]==zfile[i,4] & xfile[,10]==zfile[i,10])|(xfile[,4]==zfile[i,10] & xfile[,10]==zfile[i,4]))
    flag2=c(flag2,flagi)
    #nonAlu.CI=c(nonAlu.CI,mean(as.numeric(xfile[flagi,19])))
  }
  
  
  Alu.CI=as.numeric(xfile[flag,14])+as.numeric(xfile[flag,15])
  nonAlu.CI=as.numeric(xfile[-flag,14])+as.numeric(xfile[-flag,15])

  
  Alu.len=nrow(unique(xfile[flag,c(4,10)]))
  nonAlu.len=nrow(unique(xfile[-flag,c(4,10)]))
  Alu.CI=table(apply(xfile[flag,c(4,10)],1,function(x){paste(x,collapse = ",")}))
  nonAlu.CI=table(apply(xfile[-flag,c(4,10)],1,function(x){paste(x,collapse = ",")}))
  
  nonAlu.CI=as.numeric(xfile[-flag,14])+as.numeric(xfile[-flag,15])
  
  wilcox.test(Alu.CI[Alu.CI>0],nonAlu.CI[nonAlu.CI>0],alternative = "greater")
  boxplot(log(Alu.CI+1),log(nonAlu.CI+1))
  
  wilcox.test(as.numeric(Alu.CI),as.numeric(nonAlu.CI))
  
  
}

a=rbind(c(794,6283),c(1560,12524))
fisher.test(a)
chiqs.test(a)



#2021.05.05 Alu motif provide better predictive power of E-P interaction ?
{
  setwd("~/GM12878/3newBlast/motif-blasthit/")
  etfile=as.matrix(read.table(file="~/GM12878/2newGetEPpairsInTAD/CAGEhead/eRNA.targetProtein.txt",header = F))
  #chr1    974026  974300  GMIneRNA+1      GMIneRNA+       +       chr1    916497  919497  ENSG00000187642.5       C1orf170        -protein_coding
  
  efile=as.matrix(read.table(file="eRNA.fimo.out/fimo.tsv",header = T))
  tfile=as.matrix(read.table(file="target.fimo.out/fimo.tsv",header = T))
  #motif_id        motif_alt_id    sequence_name   start   stop    strand  score   p-value q-value matched_sequence
  #RATCTYRGCTCACTGCAACCTCYRCCTCCYR MEME-5  GMExeRNA-444::chr19:6527309-6528220(-)  625     655     -       53.9808 6.31e-20  
  
  rfile=as.matrix(read.table(file="~/GM12878/5newWithChIA-PET/EPwithChIAPETloop/eRNA.targetProtein_evidentByChIApet.txt",header = F))
  #chr10   103647679       103649515       GMIneRNA-388    GMIneRNA-       -       chr10   103602677       103605677       ENSG00000120049.14        KCNIP2  -       protein_coding
  etfileL=nrow(etfile)
  retfile=cbind(etfile[,1:6],etfile[sample(etfileL,etfileL),7:13])
  Positive.flag=which((etfile[,4]%in%rfile[,4]) & (etfile[,10]%in%rfile[,10]))
  PN=rep(0,etfileL)
  PN[Positive.flag]=1
  motif.num_score=c()
  for(i in 1:etfileL){
    eflag=grep(etfile[i,4],efile[,3])
    tflag=grep(etfile[i,10],tfile[,3])
    if(length(eflag)!=0 & length(tflag)!=0){
      e.motif=matrix(efile[eflag,c(2,7)],ncol=2)
      t.motif=matrix(tfile[tflag,c(2,7)],ncol = 2)
      et.motif=intersect(e.motif[,1],t.motif[,1])
      
      if(length(et.motif)!=0){
        e.scorek=c()
        t.scorek=c()
        for(k in 1:length(et.motif)){
          e.flagk=which(e.motif[,1]==et.motif[k])
          t.flagk=which(t.motif[,1]==et.motif[k])
          
          e.scorek=c(e.scorek,max(as.numeric(e.motif[e.flagk,2])))
          t.scorek=c(t.scorek,max(as.numeric(t.motif[t.flagk,2])))
          
        }
        outi=c(length(et.motif),mean(e.scorek),mean(t.scorek),length(eflag),length(tflag),mean(as.numeric(efile[eflag,7])),mean(as.numeric(tfile[tflag,7])))
        
      }else{
        outi=c(length(et.motif),-20,-20,length(eflag),length(tflag),mean(as.numeric(efile[eflag,7])),mean(as.numeric(tfile[tflag,7])))
        
      }
            
      
    }else if(length(eflag)!=0 & length(tflag)==0){
      outi=outi=c(0,-20,-20,length(eflag),length(tflag),mean(as.numeric(efile[eflag,7])),-20)

    }else if(length(eflag)==0 & length(tflag)!=0){
      outi=outi=c(0,-20,-20,length(eflag),length(tflag),-20,mean(as.numeric(tfile[tflag,7])))
      
    }else{
      outi=outi=c(0,-20,-20,length(eflag),length(tflag),-20,-20)
      
    }
    motif.num_score=rbind(motif.num_score,outi)
    
  }
  
  
  seqs=seq(0,30)
  seqs=seq(-20,40,5)
  e_t.num=motif.num_score[,4]+motif.num_score[,5]
  
  par(mar=c(4,4,1,1),mgp=c(1.8,0.5,0),cex=1)
  
  len=c()
  mean.seqs=rowMeans(motif.num_score[,2:3])
  for(g in 1:length(seqs)){
    leng=which(mean.seqs>=mean.seqs[g])
    len=c(len,leng)
    
  }
  
  plen=len/length(mean.seqs)
  
  PP=data.frame("predictions"=plen,"labels"=PN)
  pred<-prediction(PP$predictions,PP$labels)
  perf<-performance(pred,"tpr","fpr")
  
  pdf(file="~/eRNAprojectPlots/JGG.revise/motif.AOC.pdf",family = "ArialMT")
  plot(perf,colorize=F)
  dev.off()
    
  if(F){
    
    TPR=c()
    FPR=c()
    API=c()
    precision=c()
    for(g in 1:length(seqs)){
      
      flag1=which(motif.num_score[,2]>=seqs[g]& motif.num_score[,3]>=seqs[g])
      TF=rep(0,etfileL)
      TF[flag1]=1
      
      TP=length(which(PN==1&TF==1))
      FN=length(which(PN==1&TF==0))
      FP=length(which(PN==0&TF==1))
      TN=length(which(PN==0&TF==0))
      
      TPR=c(TPR,TP/(TP+FN))
      FPR=c(FPR,FP/(TN+FP))
      API=c(API,(TP+TN)/(TP+FN+FP+TN))
      precision=c(precision,TP/(TP+FP))
    }
    
    plot(FPR,TPR,xlim=c(0,1),ylim=c(0,1))
    plot(TPR,precision,xlim=c(0,1),ylim=c(0,1))
    par(new=T)
    plot(c(0,1),c(0,1),type = "l",lty=2,col=2,xlim=c(0,1),ylim=c(0,1),ann = F)
    
    P_motif.num_score=motif.num_score[Positive.flag,]
    N_motif.num_score=motif.num_score[-Positive.flag,]
    
    wilcox.test(P_motif.num_score[,1],N_motif.num_score[,1],alternative = "greater")#p-value = 0.01086
    boxplot(P_motif.num_score[,1],N_motif.num_score[,1])
    
    wilcox.test(e_t.num[Positive.flag],e_t.num[-Positive.flag],alternative = "greater")#p-value = 0.01086
    boxplot(P_motif.num_score[,1],N_motif.num_score[,1])
    
    hist(P_motif.num_score[,4])
    hist(N_motif.num_score[,4])
    
    
    motif.num_score.sub=motif.num_score[Positive.flag,]
    Rmotif.num_score.sub=Rmotif.num_score[Positive.flag,]
    
    seqs=seq(0,5)
    TPR=c()
    FPR=c()
    API=c()
    precision=c()
    for(g in 1:length(seqs)){
      # flag1=which(motif.num_score.sub[,1]>=seqs[g]& motif.num_score.sub[,1]>=seqs[g])
      # flag2=which(Rmotif.num_score.sub[,1]>=seqs[g]& Rmotif.num_score.sub[,1]>=seqs[g])
      flag1=which(motif.num_score.sub[,1]>=seqs[g])
      flag2=which(Rmotif.num_score.sub[,1]>=seqs[g])
      TF=rep(0,etfileL)
      TF[flag1]=1
      
      # TP=length(which(PN==1&TF==1))
      # FN=length(which(PN==1&TF==0))
      # FP=length(which(PN==0&TF==1))
      # TN=length(which(PN==0&TF==0))
      TP=length(flag1)
      FN=nrow(motif.num_score.sub)-TP
      FP=length(flag2)
      TN=nrow(Rmotif.num_score.sub)-FP
      
      TPR=c(TPR,TP/(TP+FN))
      FPR=c(FPR,FP/(TN+FP))
      API=c(API,(TP+TN)/(TP+FN+FP+TN))
      precision=c(precision,TP/(TP+FP))
    }
    
    
    plot(FPR,TPR,xlim=c(0,1),ylim=c(0,1))
    plot(TPR,precision,xlim=c(0,1),ylim=c(0,1))
    par(new=T)
    plot(c(0,1),c(0,1),type = "l",lty=2,col=2,xlim=c(0,1),ylim=c(0,1),ann = F)
    
  }  
  
}

#2021.05.05 Alu motif provide better predictive power of E-P interaction ->TargetFinder EPI
{
  setwd("~/bxDownload/hg19/K562/TargetFinder/")
  etfile=as.matrix(read.table(file="TargetFinder_K562_interactionEPI_notInteractionEPI.bed",header = F))
  #chr1    540600  541000  K562-20332      chr1    2343060 2345336 K562-20332      0       0
  
  efile=as.matrix(read.table(file="enhancer.fimo/fimo.tsv",header = T))
  tfile=as.matrix(read.table(file="promoter.fimo/fimo.tsv",header = T))
  #motif_id        motif_alt_id    sequence_name   start   stop    strand  score   p-value q-value matched_sequence
  #RATCTYRGCTCACTGCAACCTCYRCCTCCYR MEME-5  GMExeRNA-444::chr19:6527309-6528220(-)  625     655     -       53.9808 6.31e-20  
  etfileL=nrow(etfile)
  motif.num_score=c()
  for(i in 1:etfileL){
    eflag=grep(etfile[i,4],efile[,3])
    tflag=grep(etfile[i,8],tfile[,3])
    if(length(eflag)!=0 & length(tflag)!=0){
      e.motif=matrix(efile[eflag,c(2,7)],ncol=2)
      t.motif=matrix(tfile[tflag,c(2,7)],ncol = 2)
      et.motif=intersect(e.motif[,1],t.motif[,1])
      
      if(length(et.motif)!=0){
        e.scorek=c()
        t.scorek=c()
        for(k in 1:length(et.motif)){
          e.flagk=which(e.motif[,1]==et.motif[k])
          t.flagk=which(t.motif[,1]==et.motif[k])
          
          e.scorek=c(e.scorek,max(as.numeric(e.motif[e.flagk,2])))
          t.scorek=c(t.scorek,max(as.numeric(t.motif[t.flagk,2])))
          
        }
        outi=c(length(et.motif),mean(e.scorek),mean(t.scorek),length(eflag),length(tflag),mean(as.numeric(efile[eflag,7])),mean(as.numeric(tfile[tflag,7])))
        
      }else{
        outi=c(length(et.motif),-20,-20,length(eflag),length(tflag),mean(as.numeric(efile[eflag,7])),mean(as.numeric(tfile[tflag,7])))
        
      }
      
      
    }else if(length(eflag)!=0 & length(tflag)==0){
      outi=outi=c(0,-20,-20,length(eflag),length(tflag),mean(as.numeric(efile[eflag,7])),-20)
      
    }else if(length(eflag)==0 & length(tflag)!=0){
      outi=outi=c(0,-20,-20,length(eflag),length(tflag),-20,mean(as.numeric(tfile[tflag,7])))
      
    }else{
      outi=outi=c(0,-20,-20,length(eflag),length(tflag),-20,-20)
      
    }
    motif.num_score=rbind(motif.num_score,outi)
    
  }

  seqs=seq(0,30)
  seqs=seq(-20,40,5)
  e_t.num=motif.num_score[,4]+motif.num_score[,5]
  
  PN=as.numeric(etfile[,9])
  TPR=c()
  FPR=c()
  API=c()
  precision=c()
  for(g in 1:length(seqs)){
    #flag1=which(motif.num_score[,2]>=seqs[g]& motif.num_score[,3]>=seqs[g])
    flag1=which(motif.num_score[,4]>=seqs[g]& motif.num_score[,5]>=seqs[g])
    TF=rep(0,etfileL)
    TF[flag1]=1
    
    TP=length(which(PN==1&TF==1))
    FN=length(which(PN==1&TF==0))
    FP=length(which(PN==0&TF==1))
    TN=length(which(PN==0&TF==0))
    
    TPR=c(TPR,TP/(TP+FN))
    FPR=c(FPR,FP/(TN+FP))
    API=c(API,(TP+TN)/(TP+FN+FP+TN))
    precision=c(precision,TP/(TP+FP))
  }
  
  plot(FPR,TPR,xlim=c(0,1),ylim=c(0,1))
  plot(TPR,precision,xlim=c(0,1),ylim=c(0,1))
  par(new=T)
  plot(c(0,1),c(0,1),type = "l",lty=2,col=2,xlim=c(0,1),ylim=c(0,1),ann = F)
  
  P_motif.num_score=motif.num_score[Positive.flag,]
  N_motif.num_score=motif.num_score[-Positive.flag,]
  
  wilcox.test(P_motif.num_score[,1],N_motif.num_score[,1],alternative = "greater")#p-value = 0.01086
  boxplot(P_motif.num_score[,1],N_motif.num_score[,1])
  
  wilcox.test(e_t.num[Positive.flag],e_t.num[-Positive.flag],alternative = "greater")#p-value = 0.01086
  boxplot(P_motif.num_score[,1],N_motif.num_score[,1])
  
  hist(P_motif.num_score[,4])
  hist(N_motif.num_score[,4])
  
  
  motif.num_score.sub=motif.num_score[Positive.flag,]
  Rmotif.num_score.sub=Rmotif.num_score[Positive.flag,]
  
  seqs=seq(0,5)
  TPR=c()
  FPR=c()
  API=c()
  precision=c()
  for(g in 1:length(seqs)){
    # flag1=which(motif.num_score.sub[,1]>=seqs[g]& motif.num_score.sub[,1]>=seqs[g])
    # flag2=which(Rmotif.num_score.sub[,1]>=seqs[g]& Rmotif.num_score.sub[,1]>=seqs[g])
    flag1=which(motif.num_score.sub[,1]>=seqs[g])
    flag2=which(Rmotif.num_score.sub[,1]>=seqs[g])
    TF=rep(0,etfileL)
    TF[flag1]=1
    
    # TP=length(which(PN==1&TF==1))
    # FN=length(which(PN==1&TF==0))
    # FP=length(which(PN==0&TF==1))
    # TN=length(which(PN==0&TF==0))
    TP=length(flag1)
    FN=nrow(motif.num_score.sub)-TP
    FP=length(flag2)
    TN=nrow(Rmotif.num_score.sub)-FP
    
    TPR=c(TPR,TP/(TP+FN))
    FPR=c(FPR,FP/(TN+FP))
    API=c(API,(TP+TN)/(TP+FN+FP+TN))
    precision=c(precision,TP/(TP+FP))
  }
  
  
  plot(FPR,TPR,xlim=c(0,1),ylim=c(0,1))
  plot(TPR,precision,xlim=c(0,1),ylim=c(0,1))
  par(new=T)
  plot(c(0,1),c(0,1),type = "l",lty=2,col=2,xlim=c(0,1),ylim=c(0,1),ann = F)
  
  
  
}

#2021.05.12 Alu density in TSS flanking 100kb 
{
  AluCount=function(TSS,alu.flag){
    seqs=seq(TSS-100000,TSS+100000,1000)
    SE=cbind(seqs-2500,seqs+2500)
    SEL=nrow(SE)
    outk=c()
    for(k in 1:SEL){
      cha1=SE[k,2]-alu.S[alu.flag]
      cha2=alu.E[alu.flag]-SE[k,1]
      countk=length(which((cha1>0) & (cha2>0)))
      outk=c(outk,countk)
    }
    return(outk)
  }
  
  
  
  setwd("~/GM12878/2newGetEPpairsInTAD/CAGEhead/")
  pfile=as.matrix(read.table(file="./proteinTarget.bed",header = F))
  #chr10	100162892	100163920	GMIneRNA-386	GMIneRNA-	-
  Alufile= afile <- as.matrix(read.table(file="/data/baix/bxDownload/hg19/Alu/AluInhg19_Repbase_chr1-23.bed",header=F,sep="\t"))
  #chr1 26791 27053  AluSp SINE/Alu +
  flagP=which(pfile[,6]=="+")
  flagM=which(pfile[,6]=="-")
  pfileM=pfile[flagM,]
  pfileP=pfile[flagP,]

  alu.S=as.numeric(Alufile[,2])
  alu.E=as.numeric(Alufile[,3])
  TSS.P=(as.numeric(pfile[flagP,2])+2000)
  TSS.M=(as.numeric(pfile[flagM,2])+1000)

  chr=unique(pfile[,1])
  outseq.p=c()
  outseq.m=c()
  for(i in 1:length(chr)){
    alu.flag=which(Alufile[,1]==chr[i])
    
    pflag=which(pfileP[,1]==chr[i])
    pflagL=length(pflag)
    for(j in 1:pflagL){
      TSSj=TSS.P[pflag[j]]
      outseq.pj=AluCount(TSSj,alu.flag)
      outseq.p=rbind(outseq.p,outseq.pj)
    }
    
    mflag=which(pfileM[,1]==chr[i])
    mflagL=length(mflag)
    for(j in 1:mflagL){
      TSSj=TSS.M[mflag[j]]
      outseq.mj=AluCount(TSSj,alu.flag)
      outseq.m=rbind(outseq.m,outseq.mj[length(outseq.mj):1])
    }
    
    
  }
  
  outseq=rbind(outseq.p,outseq.m) 
  outseq.mean=colMeans(outseq)
  par(mfrow=c(1,1))
  plot(-100:100,outseq.mean,type = "l")
  
}

#2021.05.13 k562 DRIPc-seq in Enhancer-> wether enriched in eRNA Alu
{
  setwd("~/newCell4Type/2getEP/K562/DRIPc-seq")
  e.file=as.matrix(read.table(file="../eRNA_withAlu.bed",header = F))
  ealu.file=as.matrix(read.table(file="../eRNA.Alu.txt",header = F))
  filenames=rbind(c("eRNA_withAlu_GSM4801435_K562_DRIPc_WT_LS60A_rep1_pos_bigwig.txt","eRNA_withAlu_GSM4801435_K562_DRIPc_WT_LS60A_rep1_neg_bigwig.txt"),
                  c("eRNA_withAlu_GSM4801436_K562_DRIPc_WT_LS60B_rep2_pos_bigwig.txt","eRNA_withAlu_GSM4801436_K562_DRIPc_WT_LS60B_rep2_neg_bigwig.txt"))
  counts=c()
  for(k in 1:2){
    DRIPc.fileP=file(filenames[k,1],"r")
    DRIPc.fileM=file(filenames[k,2],"r")
    
    DRIPc.lineP=readLines(DRIPc.fileP)
    DRIPc.lineM=readLines(DRIPc.fileM)
    close(DRIPc.fileP)
    close(DRIPc.fileM)
    
    value=c()
    len=c()
    e.fileL=nrow(e.file)
    for(i in 1: e.fileL){
      ealu.flag=which(ealu.file[,4]==e.file[i,4])
      ealu.S=as.numeric(ealu.file[ealu.flag,8])-as.numeric(ealu.file[ealu.flag,2])
      ealu.E=as.numeric(ealu.file[ealu.flag,9])-as.numeric(ealu.file[ealu.flag,2])
      
      if(e.file[i,6]=="+"){
        splitLineP=strsplit(DRIPc.lineP[i],"\t")[[1]]
        Alui=c()
        for(j in 1:length(ealu.flag)){
          ealu.Sj=max(1,ealu.S[j])
          ealu.Ej=min(ealu.E[j],length(splitLineP))
          Alui=c(Alui,ealu.Sj:ealu.Ej)
        }
        nonAlui=setdiff(1:length(splitLineP),Alui)
        value=rbind(value,c(median(as.numeric(splitLineP[Alui])),median(as.numeric(splitLineP[nonAlui]))))
        len=rbind(len,c(length(Alui),length(nonAlui)))
        
      }
      
      if(e.file[i,6]=="-"){
        splitLineM=strsplit(DRIPc.lineM[i],"\t")[[1]]
        Alui=c()
        for(j in 1:length(ealu.flag)){
          ealu.Sj=max(1,ealu.S[j])
          ealu.Ej=min(ealu.E[j],length(splitLineM))
          Alui=c(Alui,ealu.Sj:ealu.Ej)
        }
        nonAlui=setdiff(1:length(splitLineM),Alui)
        value=rbind(value,c(median(as.numeric(splitLineM[Alui])),median(as.numeric(splitLineM[nonAlui]))))
        len=rbind(len,c(length(Alui),length(nonAlui)))
        
      }
    }
    
    norm.alu.value=value[,1]*1000/len[,1]
    norm.nonalu.value=value[,2]*1000/len[,2]
    
    seqs=seq(0,14,0.5)
    a1=hist(log2(norm.alu.value+1),breaks = seqs,plot=F)
    a2=hist(log2(norm.nonalu.value+1),breaks = seqs,plot=F)
    a=cbind(a1$counts/sum(a1$counts),a2$counts/sum(a2$counts))
    counts=cbind(counts,a)    
  }
  
  b1=c()
  for(k in 1:4){
    b1=cbind(b1,cumsum(counts[,k]))
  }
  
  b=apply(b1,1,function(x){c(mean(x[c(1,3)]),mean(x[c(2,4)]))})
  sd1=apply(b1,1,function(x){c(sd(x[c(1,3)]),sd(x[c(2,4)]))})
  
  
  ymax=1
  ymin=0.2
  
  pdf(file="~/eRNAprojectPlots/figure6/K562_DRIPc_IneRNA.alu_nonAlu-errorbar.pdf",family = "ArialMT")
  {
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    par(mar=c(4,4,1,1),cex=1.1,mgp=c(1.8,0.5,0))
    plot(1:nrow(counts[-1,]),b[1,-nrow(counts)],xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts[-1,]),b[2,-nrow(counts)],xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    
    legend("bottomright",inset=0.05,legend = c("Alu elements","non-Alu region"),col =cols[1:2],lty = c(1, 1),cex=1,lwd = 2,bty = "n")
    axis(1,at=seq(1,nrow(counts[-1,]),4),labels=seq(0,13,2), col = "black", lty = 1, lwd = 1)
    
    for(i in 1:dim(counts)[1]){
      arrows(i,b[1,i]+sd1[1,i],i, b[1,i]-sd1[1,i], angle=90, code=3, length=0.03,lwd=1.5)
    }
    for(i in 1:dim(counts)[1]){
      arrows(i,b[2,i]+sd1[2,i],i, b[2,i]-sd1[2,i], angle=90, code=3, length=0.03,lwd=1.5)
    }
    text(17,0.5,paste("#eRNA=",nrow(e.file)))
  }
  dev.off()
  
  pdf(file="~/eRNAprojectPlots/figure6/K562_DRIPc_IneRNA.alu_nonAlu_rep1-2.pdf",family = "ArialMT")
  {
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    par(mar=c(4,4,1,1),cex=1.1,mgp=c(1.8,0.5,0))
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),1],lty=1,xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),2],lty=1,xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),3],lty=2,xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),4],lty=2,xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    
    legend("bottomright",inset=0.05,legend = rep(c("Alu elements","non-Alu region"),each=2),col =rep(cols[1:2],each=2),lty = c(1, 2,1,2),cex=1,lwd = 2,bty = "n")
    axis(1,at=seq(1,nrow(counts[-1,]),4),labels=seq(0,13,2), col = "black", lty = 1, lwd = 1)
    text(17,0.5,paste("#eRNA=",nrow(e.file)))
  }
  dev.off()

}

#2021.05.13 HeLa DRIPc-seq in Enhancer-> wether enriched in eRNA Alu
{
  setwd("~/newCell4Type/2getEP/HeLa-S3/DRIPc-seq/")
  e.file=as.matrix(read.table(file="../eRNA_withAlu.bed",header = F))
  ealu.file=as.matrix(read.table(file="../eRNA.Alu.txt",header = F))
  filenames=rbind(c("eRNA_withAlu_GSM4801428_HeLa_DRIPc_WT_LS61A_rep1_neg_bigwig.txt","eRNA_withAlu_GSM4801428_HeLa_DRIPc_WT_LS61A_rep1_pos_bigwig.txt"),
                  c("eRNA_withAlu_GSM4801429_HeLa_DRIPc_WT_LS61C_rep2_neg_bigwig.txt","eRNA_withAlu_GSM4801429_HeLa_DRIPc_WT_LS61C_rep2_pos_bigwig.txt"),
                  c("eRNA_withAlu_GSM4801430_HeLa_DRIPc_WT_LS61H_rep3_neg_bigwig.txt","eRNA_withAlu_GSM4801430_HeLa_DRIPc_WT_LS61H_rep3_pos_bigwig.txt"))
  counts=c()
  for(k in 1:3){
    DRIPc.fileP=file(filenames[k,1],"r")
    DRIPc.fileM=file(filenames[k,2],"r")
    
    DRIPc.lineP=readLines(DRIPc.fileP)
    DRIPc.lineM=readLines(DRIPc.fileM)
    close(DRIPc.fileP)
    close(DRIPc.fileM)
    
    value=c()
    len=c()
    e.fileL=nrow(e.file)
    for(i in 1: e.fileL){
      ealu.flag=which(ealu.file[,4]==e.file[i,4])
      ealu.S=as.numeric(ealu.file[ealu.flag,9])-as.numeric(ealu.file[ealu.flag,2])
      ealu.E=as.numeric(ealu.file[ealu.flag,10])-as.numeric(ealu.file[ealu.flag,2])
      
      if(e.file[i,6]=="+"){
        splitLineP=strsplit(DRIPc.lineP[i],"\t")[[1]]
        Alui=c()
        for(j in 1:length(ealu.flag)){
          ealu.Sj=max(1,ealu.S[j])
          ealu.Ej=min(ealu.E[j],length(splitLineP))
          Alui=c(Alui,ealu.Sj:ealu.Ej)
        }
        nonAlui=setdiff(1:length(splitLineP),Alui)
        value=rbind(value,c(median(as.numeric(splitLineP[Alui])),median(as.numeric(splitLineP[nonAlui]))))
        len=rbind(len,c(length(Alui),length(nonAlui)))
        
      }
      
      if(e.file[i,6]=="-"){
        splitLineM=strsplit(DRIPc.lineM[i],"\t")[[1]]
        Alui=c()
        for(j in 1:length(ealu.flag)){
          ealu.Sj=max(1,ealu.S[j])
          ealu.Ej=min(ealu.E[j],length(splitLineM))
          Alui=c(Alui,ealu.Sj:ealu.Ej)
        }
        nonAlui=setdiff(1:length(splitLineM),Alui)
        value=rbind(value,c(median(as.numeric(splitLineM[Alui])),median(as.numeric(splitLineM[nonAlui]))))
        len=rbind(len,c(length(Alui),length(nonAlui)))
        
      }
    }
    
    norm.alu.value=value[,1]*1000/len[,1]
    norm.nonalu.value=value[,2]*1000/len[,2]
    
    seqs=seq(0,10,0.5)
    a1=hist(log2(norm.alu.value+1),breaks = seqs,plot=F)
    a2=hist(log2(norm.nonalu.value+1),breaks = seqs,plot=F)
    a=cbind(a1$counts/sum(a1$counts),a2$counts/sum(a2$counts))
    counts=cbind(counts,a)    
  }
  
  b1=c()
  for(k in 1:6){
    b1=cbind(b1,cumsum(counts[,k]))
  }
  
  b=apply(b1,1,function(x){c(mean(x[c(1,3,5)]),mean(x[c(2,4,6)]))})
  sd1=apply(b1,1,function(x){c(sd(x[c(1,3,5)]),sd(x[c(2,4,6)]))})
  
  
  ymax=1
  ymin=0.8
  
  pdf(file="~/eRNAprojectPlots/figure6/HeLa_DRIPc_IneRNA.alu_nonAlu-errorbar.pdf",family = "ArialMT")
  {
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    par(mar=c(4,4,1,1),cex=1.1,mgp=c(1.8,0.5,0))
    plot(1:nrow(counts),b[1,],xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts),b[2,],xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    
    legend("bottomright",inset=0.05,legend = c("Alu elements","non-Alu region"),col =cols[1:2],lty = c(1, 1),cex=1,lwd = 2,bty = "n")
    axis(1,at=seq(1,nrow(counts)+1,6),labels=seq(0,10,3), col = "black", lty = 1, lwd = 1)
    
    for(i in 1:dim(counts)[1]){
      arrows(i,b[1,i]+sd1[1,i],i, b[1,i]-sd1[1,i], angle=90, code=3, length=0.03,lwd=1.5)
    }
    for(i in 1:dim(counts)[1]){
      arrows(i,b[2,i]+sd1[2,i],i, b[2,i]-sd1[2,i], angle=90, code=3, length=0.03,lwd=1.5)
    }
    
    text(15,0.88,paste("#eRNA=",nrow(e.file)))
  }
  dev.off()
  
  pdf(file="~/eRNAprojectPlots/figure6/HeLa_DRIPc_IneRNA.alu_nonAlu_rep1-2.pdf",family = "ArialMT")
  {
    cols=c("orangered","steelblue2","goldenrod1","darkseagreen3")
    par(mar=c(4,4,1,1),cex=1.1,mgp=c(1.8,0.5,0))
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),1],lty=1,xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),2],lty=1,xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),3],lty=2,xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),4],lty=2,xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),5],lty=2,xaxt="n",type="l",cex.lab=1.2,cex.axis=1,col=cols[1],ylim=c(ymin,ymax),lwd=2,xlab = "log2 of DRIPc-seq signal values per 1kb",ylab = "Cumulative frequency",main = "")
    par(new=T)
    plot(1:nrow(counts[-1,]),b1[-nrow(counts),6],lty=2,xaxt="n",yaxt="n",type="l",col=cols[2],ylim=c(ymin,ymax),lwd=2,xlab="",ylab ="",main="")
    
    legend("bottomright",inset=0.05,legend = rep(c("Alu elements","non-Alu region"),each=2),col =rep(cols[1:2],each=2),lty = c(1, 2,1,2),cex=1,lwd = 2,bty = "n")
    axis(1,at=seq(1,nrow(counts[-1,]),4),labels=seq(0,13,2), col = "black", lty = 1, lwd = 1)
    
  }
  dev.off()
  
}


#2021.07.05 read through result downstream gene expression increase
#bigwig TPM -no use
{
  setwd("~/LuChengjiang-dnCLIP-seq/RNA-seq/readThrough_downStreamGene")
  xfile=as.matrix(read.table(file="../bash/fileNames.txt"))
  plus.file=as.matrix(read.table(file="/data/baix/bxDownload/hg19/refSeq/mergedGene/upGene_intergenic_downGene_plus.bed"))
  minus.file=as.matrix(read.table(file="/data/baix/bxDownload/hg19/refSeq/mergedGene/upGene_intergenic_downGene_minus.bed"))
  plus.fileL=nrow(plus.file)
  minus.fileL=nrow(minus.file)
  Pbigwig=c()
  Mbigwig=c()
  
  for(i in 1:nrow(xfile)){
    Pbigwig.i=rep(0,plus.fileL)
    Mbigwig.i=rep(0,minus.fileL)
    
    mfile=as.matrix(read.table(file=paste(xfile[i,1],"_minus.mat",sep=""),header = F,skip = 1))
    pfile=as.matrix(read.table(file=paste(xfile[i,1],"_plus.mat",sep=""),header = F,skip=1))
    #chr1    89295   133566  ENSG00000238009.2+1     0.0     - 
    # m=apply(mfile[,1:6],1,function(x){as.numeric(strsplit(x[4],"\\+")[[1]][2])})
    # p=apply(pfile[,1:6],1,function(x){as.numeric(strsplit(x[4],"\\+")[[1]][2])})
    m=unlist(apply(mfile[,1:6],1,function(x){which(minus.file[,4]==x[4] & minus.file[,2]==x[2])}))
    p=unlist(apply(pfile[,1:6],1,function(x){which(plus.file[,4]==x[4] & plus.file[,2]==x[2])}))
    
    msum=apply(mfile[,-c(1:6)],1,function(x){sum(as.numeric(x))})
    psum=apply(pfile[,-c(1:6)],1,function(x){sum(as.numeric(x))})
    Pbigwig.i[p]=psum
    Mbigwig.i[m]=msum
    Pbigwig=cbind(Pbigwig,Pbigwig.i)
    Mbigwig=cbind(Mbigwig,Mbigwig.i)
    
  }
  
  plus.poi=as.numeric(apply(plus.file, 1, function(x){strsplit(x[4],"\\+")[[1]][2]}))
  minus.poi=as.numeric(apply(minus.file, 1, function(x){strsplit(x[4],"\\+")[[1]][2]}))
  poi=c(unique(plus.poi),unique(minus.poi))
  Con_rep1=cbind(rbind(matrix(Pbigwig[,1],ncol=3,byrow = T),matrix(Mbigwig[,1],ncol=3,byrow = T)[,3:1]),poi)      
  Con_rep2=cbind(rbind(matrix(Pbigwig[,2],ncol=3,byrow = T),matrix(Mbigwig[,2],ncol=3,byrow = T)[,3:1]),poi)  
  KD1_rep1=cbind(rbind(matrix(Pbigwig[,3],ncol=3,byrow = T),matrix(Mbigwig[,3],ncol=3,byrow = T)[,3:1]),poi)  
  KD1_rep2=cbind(rbind(matrix(Pbigwig[,4],ncol=3,byrow = T),matrix(Mbigwig[,4],ncol=3,byrow = T)[,3:1]),poi)
  KD2_rep1=cbind(rbind(matrix(Pbigwig[,5],ncol=3,byrow = T),matrix(Mbigwig[,5],ncol=3,byrow = T)[,3:1]),poi)  
  KD2_rep2=cbind(rbind(matrix(Pbigwig[,6],ncol=3,byrow = T),matrix(Mbigwig[,6],ncol=3,byrow = T)[,3:1]),poi)
  
  cutoff1=5
  cutoff2=0.1
  flag1=which(Con_rep1[,1]>cutoff1 & Con_rep2[,1]>cutoff1 & KD1_rep1[,1]>cutoff1
              &KD1_rep2[,1]>cutoff1 & KD2_rep1[,1]>cutoff1 & KD2_rep2[,1]>cutoff1)
  flag2=which(KD1_rep1[,2]>cutoff2 & KD1_rep2[,2]>cutoff2 & KD2_rep1[,2]>cutoff2 & KD2_rep2[,2]>cutoff2)
  flag3=which(KD1_rep1[,3]>cutoff2 & KD1_rep2[,3]>cutoff2 & KD2_rep1[,3]>cutoff2 & KD2_rep2[,3]>cutoff2)
  
  flag=intersect(intersect(flag1,flag2),flag3)
  
  Con=(Con_rep1[flag,]+Con_rep2[flag,])/2
  KD1=(KD1_rep1[flag,]+KD1_rep2[flag,])/2
  KD2=(KD2_rep1[flag,]+KD2_rep2[flag,])/2
  
  fd1=(KD1+1)/(Con+1)
  fd2=(KD2+1)/(Con+1)
  
  index.fd1=fd1[,2:3]/fd1[,1]
  index.fd2=fd2[,2:3]/fd2[,1]
  
  
  flag1=which(index.fd1[,1]>1.1 & index.fd1[,2]>1.1)
  flag2=which(index.fd2[,1]>1.1 & index.fd2[,2]>1.1)
  flag3=which(KD1[,2]>0 & KD1[,3]>0 &KD2[,2]>0 & KD2[,3]>0)
  
  flag=intersect(intersect(flag1,flag2),flag3)
  
  con.fd=(Con[flag,2:3]+1)/(Con[flag,1]+1)
  KD1.fd=(KD1[flag,2:3]+1)/(KD1[flag,1]+1)
  KD2.fd=(KD2[flag,2:3]+1)/(KD2[flag,1]+1)
  
  
  con.fd=(Con[flag[flag.l],2:3]+1)/(Con[flag[flag.l],1]+1)
  KD1.fd=(KD1[flag[flag.l],2:3]+1)/(KD1[flag[flag.l],1]+1)
  KD2.fd=(KD2[flag[flag.l],2:3]+1)/(KD2[flag[flag.l],1]+1)
  
  mycol=c("black","sienna2","skyblue1")
  #TPM Plot termination ratio
  pdf(file="~/LuChengjiang-dnCLIP-seq/Figures/RNAseq-readThrough_TerminationRatio_bigwig_neighborGene.pdf",family = "ArialMT")
  {
    par(mar=c(3,3,1,2),mgp=c(1.8,0.5,0))
    breaks=seq(-5,2,0.2)
    
    a1=hist(log10(con.fd[,1]),breaks = breaks,plot = F)
    a2=hist(log10(KD1.fd[,1]),breaks = breaks,plot = F)
    a3=hist(log10(KD2.fd[,1]),breaks = breaks,plot = F)
    
    b=cbind(cumsum(a1$counts)/sum(a1$counts),cumsum(a2$counts)/sum(a2$counts),
            cumsum(a3$counts)/sum(a3$counts))
    
    a12=hist(log10(con.fd[,2]),breaks = breaks,plot = F)
    a22=hist(log10(KD1.fd[,2]),breaks = breaks,plot = F)
    a32=hist(log10(KD2.fd[,2]),breaks = breaks,plot = F)
    
    b2=cbind(cumsum(a12$counts)/sum(a12$counts),cumsum(a22$counts)/sum(a22$counts),
             cumsum(a32$counts)/sum(a32$counts))
    
    xmax=max(breaks)-0.1
    xmin=-4
    seqs=5:(length(breaks)-1)
    plot(breaks[seqs],b[seqs,1],type="l",col=mycol[1],xlim=c(xmin,xmax),ylim=c(0,1),lwd=2,xlab="Termination ratio (log10)",ylab="Cumulative distribution")
    par(new=T)
    plot(breaks[seqs],b[seqs,2],type="l",col=mycol[2],xlim=c(xmin,xmax),ylim=c(0,1),lwd=2,ann=F,axes=F)
    par(new=T)
    plot(breaks[seqs],b[seqs,3],type="l",col=mycol[3],xlim=c(xmin,xmax),ylim=c(0,1),lwd=2,ann=F,axes=F)
    
    par(new=T)
    
    seqs=5:(length(breaks)-1)
    plot(breaks[seqs],b2[seqs,1],type="l",lty=2,col=mycol[1],xlim=c(xmin,xmax),ylim=c(0,1),lwd=2,ann=F,axes=F)
    par(new=T)
    plot(breaks[seqs],b2[seqs,2],type="l",lty=2,col=mycol[2],xlim=c(xmin,xmax),ylim=c(0,1),lwd=2,ann=F,axes=F)
    par(new=T)
    plot(breaks[seqs],b2[seqs,3],type="l",lty=2,col=mycol[3],xlim=c(xmin,xmax),ylim=c(0,1),lwd=2,ann=F,axes=F)
    
    legend("topleft",inset=0.05,legend = rep(c("Ctrl","KD1","KD2"),2),col = mycol,lty = c(rep(1,3),rep(2,3)),cex=1,lwd=2,bty = "n",ncol = 2)#title = "IntergenicR/UpstreamG DownstreamG/UpstreamG"
    text(0,0.4,paste("n=",length(flag)))
    text(0,0.3,"P< 2.2e-16")
    
    #t.test(log10(con.fd1[a]),log10(KD.fd1[a]))
    p1=wilcox.test(rep(log10(con.fd[,1]),2),c(log10(KD1.fd[,1]),log10(KD2.fd[,1])),paired = T,alternative = "less")$p.value
    p2=wilcox.test(rep(log10(con.fd[,2]),2),c(log10(KD1.fd[,2]),log10(KD2.fd[,2])),paired = T,alternative = "less")$p.value
    wilcox.test(log10(con.fd[,1]),log10(KD2.fd[,1]),paired = T,alternative = "less")$p.value
    
    #text(-1,0.4,paste("P=",format(p,scientific = T,digits = 2)))
    
    
  }
  dev.off()
  
  
  mycol=c("gray","sienna2","skyblue1")
  #read throght-> TPM boxplot
  pdf(file="~/LuChengjiang-dnCLIP-seq/Figures/RNAseq-readThrough_TPM.bigwig_neighborGene_boxplot.pdf",family = "ArialMT")
  {
    boxplot(log10(Con[flag,1:3]+1),
            boxwex = 0.25, at = 1:3 - 0.3,
            col = mycol[1],range=2.5,
            main = "",xlab = "",
            ylab = "TPM (log10)",#border="white",
            xlim = c(0.5, 3.5), ylim = c(-0.1, 5.1), yaxs = "i",xaxt="n")
    boxplot(log10(KD1[flag,1:3]+1), add = TRUE,range=2.5,
            boxwex = 0.25, at = 1:3 ,xaxt="n",yaxt="n",
            col = mycol[2])
    boxplot(log10(KD2[flag,1:3]+1), add = TRUE,range=2.5,
            boxwex = 0.25, at = 1:3 + 0.3,xaxt="n",yaxt="n",
            col =mycol[3] )
    legend("topright",inset = 0.05, c("Ctrl", "KD1","KD2"),fill = mycol,bty="n",ncol = 3)
    axis(1,at=1:3,labels=c("upstream gene","intergenic","downstream gene"), col = 1, lty = 1, lwd = 1,cex.axis=1.2)
    
    # p1=wilcox.test(Con[flag,1],c(KD1[flag,1],KD2[flag,1]),alternative = "less")$p.value
    # p2=wilcox.test(Con[flag,2],c(KD1[flag,2],KD2[flag,2]),alternative = "less")$p.value
    # p3=wilcox.test(Con[flag,3],c(KD1[flag,3],KD2[flag,3]),alternative = "less")$p.value
    p1=wilcox.test(rep(Con[flag,1],2),c(KD1[flag,1],KD2[flag,1]),paired = T,alternative = "less")$p.value
    p2=wilcox.test(rep(Con[flag,2],2),c(KD1[flag,2],KD2[flag,2]),paired = T,alternative = "less")$p.value
    p3=wilcox.test(rep(Con[flag,3],2),c(KD1[flag,3],KD2[flag,3]),paired = T,alternative = "less")$p.value
    
    p=format(c(p1,p2,p3),scientific = T,digits = 3)
    
    for(i in 1:3){
      y=max(log10(KD1[flag,i]+1))
      lines(c(i,i+0.3),c(y+0.15,y+0.15))
      lines(c(i-0.3,i+0.15),c(y+0.35,y+0.35))
      lines(c(i-0.3,i-0.3),c(y+0.15,y+0.35))
      lines(c(i+0.15,i+0.15),c(y+0.15,y+0.35))
      text(i,y+0.5,paste("P=",p[i],sep=""))
    }
    
  }
  dev.off()
  
  
  pout.name=plus.file[which(plus.poi%in%Con[flag,4]),]
  mout.name=minus.file[which(minus.poi%in%Con[flag,4]),]
  
  #a=c(plus.poi[which(plus.poi%in%Con[flag,4])],minus.poi[which(minus.poi%in%Con[flag,4])])
  #plot.png
  up_inter_down.bwget=function(input){
    mfile=as.matrix(read.table(file=paste(input,"_minus.mat",sep=""),header = F,skip = 1))
    pfile=as.matrix(read.table(file=paste(input,"_plus.mat",sep=""),header = F,skip=1))
    #chr1    89295   133566  ENSG00000238009.2+1     0.0     - 
    outm=c()
    position=1
    
    for(m in 1:nrow(mout.name)){
      flagi=which(mfile[,2]==mout.name[m,2] & mfile[,4]==mout.name[m,4])
      
      if(length(flagi)!=0){
        subout=as.numeric(mfile[flagi,-c(1:6)])
      }else{
        subout=rep(0,50)
      }
      if(position==1){
        outi=subout
        position=position+1
      }else if(position==2){
        outi=c(outi,subout)
        position=position+1
      }else{
        outi=c(outi,subout)
        outm=rbind(outm,outi)
        position=1
      }
      
      
    }
    
    outp=c()
    position=1
    for(p in 1:nrow(pout.name)){
      flagi=which(pfile[,2]==pout.name[p,2] & pfile[,4]==pout.name[p,4])
      
      if(length(flagi)!=0){
        subout=as.numeric(pfile[flagi,-c(1:6)])
      }else{
        subout=rep(0,50)
      }
      if(position==1){
        outi=subout
        position=position+1
      }else if(position==2){
        outi=c(outi,subout)
        position=position+1
      }else{
        outi=c(outi,subout)
        outp=rbind(outp,outi)
        position=1
      }
      
      
    }
    
    
    return(rbind(outp,outm[,ncol(outm):1]))
    
  }
  
  Con_rep1_bw=up_inter_down.bwget(xfile[1])
  Con_rep2_bw=up_inter_down.bwget(xfile[2])
  KD1_rep1_bw=up_inter_down.bwget(xfile[3])
  KD1_rep2_bw=up_inter_down.bwget(xfile[4])
  KD2_rep1_bw=up_inter_down.bwget(xfile[5])
  KD2_rep2_bw=up_inter_down.bwget(xfile[6])
  
  Con_bw=(Con_rep1_bw+Con_rep2_bw)/2
  KD1_bw=(KD1_rep1_bw+KD1_rep2_bw)/2
  KD2_bw=(KD2_rep1_bw+KD2_rep2_bw)/2
  cutoff3=0
  KD1.TF=apply(KD1_bw[,51:100],1,function(x){length(which(x>cutoff3))>10})
  KD2.TF=apply(KD2_bw[,51:100],1,function(x){length(which(x>cutoff3))>10})
  flag.l=which(KD1.TF & KD2.TF)
  
  pdf(file="~/LuChengjiang-dnCLIP-seq/Figures/RNAseq-readThrough_TPM.bigwig_neighborGene_plot.pdf",family = "ArialMT")
  {
    
    Con.median=apply(Con_bw[flag.l,],2,median)
    KD1.median=apply(KD1_bw[flag.l,],2,median)
    KD2.median=apply(KD2_bw[flag.l,],2,median)
    
    mycol=c("black","sienna2","skyblue1")
    ymin=-6
    ymax=0
    xmin=1
    xmax=length(Con.median)
    plot(1:xmax,log2(Con.median),ylim=c(ymin,ymax),xlim=c(xmin,xmax),lwd=2,col=mycol[1],type = "l",xaxt="n",xlab="",ylab = "Median TPM (log2)")
    par(new=T)
    plot(1:xmax,log2(KD1.median),ylim=c(ymin,ymax),xlim=c(xmin,xmax),lwd=2,col=mycol[2],type = "l",axes = F,xlab="",ylab = "")
    par(new=T)
    plot(1:xmax,log2(KD2.median),ylim=c(ymin,ymax),xlim=c(xmin,xmax),lwd=2,col=mycol[3],type = "l",axes = F,xlab="",ylab = "")
    axis(1,at=c(1,50,101,150),labels=c("TSS","TTS","TSS","TTS"), las=1,lty = 1, lwd = 1)
    legend("top",inset=0.05,legend = c("Con","KD1","KD2"),col = mycol,lty = 1,cex=1,lwd=2,bty = "n")
    text(70,-3,paste("n=",length(flag),sep=""))
    abline(v=50,lty=2,col="gray")
    abline(v=101,lty=2,col="gray")
    
  }
  dev.off()
  
  write.table(rbind(pout.name,mout.name),file="merged.upGene_intergenic_downGene_indexLarger1.1.bed",quote = F,col.names = F,row.names = F,sep="\t")
  
}


{
  a=c(2068086,
  2001318,
  1343926,
  644958,
  1343737,
  1192952)
  b=c(7390847,
    5216526,
    6998487,
    2014635,
    3752781,
    3332059)
  cc=c(6056270,
      4492482,
      5596717,
      1608632,
      3094971,
      2790609)
  paste("(",round(a/b,4)*100,"%)",sep="")
  
  
}

aa=cbind(a,b,cc)
bb=as.data.frame(aa)
colnames(bb)=c("a","b","c")

group_da=data.frame(type=rownames(group),group$cluster)
meta_da=data.frame(type=rownames(meta),meta)
merge_da=left_join(meta_da,group_da,by=type)


{
  ## boxplot on a formula:
  boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
  # *add* notches (somewhat funny here <--> warning "notches .. outside hinges"):
  boxplot(count ~ spray, data = InsectSprays,
          notch = TRUE, add = TRUE, col = "blue")
  
  boxplot(decrease ~ treatment, data = OrchardSprays, col = "bisque",
          log = "y")
  ## horizontal=TRUE, switching  y <--> x :
  boxplot(decrease ~ treatment, data = OrchardSprays, col = "bisque",
          log = "x", horizontal=TRUE)
  
  rb <- boxplot(decrease ~ treatment, data = OrchardSprays, col = "bisque")
  title("Comparing boxplot()s and non-robust mean +/- SD")
  mn.t <- tapply(OrchardSprays$decrease, OrchardSprays$treatment, mean)
  sd.t <- tapply(OrchardSprays$decrease, OrchardSprays$treatment, sd)
  xi <- 0.3 + seq(rb$n)
  points(xi, mn.t, col = "orange", pch = 18)
  arrows(xi, mn.t - sd.t, xi, mn.t + sd.t,
         code = 3, col = "pink", angle = 75, length = .1)
  
  ## boxplot on a matrix:
  mat <- cbind(Uni05 = (1:100)/21, Norm = rnorm(100),
               `5T` = rt(100, df = 5), Gam2 = rgamma(100, shape = 2))
  boxplot(mat) # directly, calling boxplot.matrix()
  
  ## boxplot on a data frame:
  df. <- as.data.frame(mat)
  par(las = 1) # all axis labels horizontal
  boxplot(df., main = "boxplot(*, horizontal = TRUE)", horizontal = TRUE)
  
  ## Using 'at = ' and adding boxplots -- example idea by Roger Bivand :
  boxplot(len ~ dose, data = ToothGrowth,
          boxwex = 0.25, at = 1:3 - 0.2,
          subset = supp == "VC", col = "yellow",
          main = "Guinea Pigs' Tooth Growth",
          xlab = "Vitamin C dose mg",
          ylab = "tooth length",
          xlim = c(0.5, 3.5), ylim = c(0, 35), yaxs = "i")
  boxplot(len ~ dose, data = ToothGrowth, add = TRUE,
          boxwex = 0.25, at = 1:3 + 0.2,
          subset = supp == "OJ", col = "orange")
  legend(2, 9, c("Ascorbic acid", "Orange juice"),
         fill = c("yellow", "orange"))
  
  ## With less effort (slightly different) using factor *interaction*:
  boxplot(len ~ dose:supp, data = ToothGrowth,
          boxwex = 0.5, col = c("orange", "yellow"),
          main = "Guinea Pigs' Tooth Growth",
          xlab = "Vitamin C dose mg", ylab = "tooth length",
          sep = ":", lex.order = TRUE, ylim = c(0, 35), yaxs = "i")
  
  ## more examples in  help(bxp)
}


rfile=as.matrix(read.table(file="~/LuChengjiang-dnCLIP-seq/RNA-seq/1test_TPM.bed"),header=F)
len=as.numeric(rfile[,3])-as.numeric(rfile[,2])
bw=apply(rfile[,6:8],1,function(x){if(x[1]=="+"){as.numeric(x[2])}else{as.numeric(x[3])}})
a=bw/500

N=c(21.0946405006086,  0 ,0)
O=c(23.4892058470132,33.5785049251802,54.7316083710994)

(mean(N)-mean(O))/sqrt((sd(N)/3)^2+(sd(O)/3)^2)






