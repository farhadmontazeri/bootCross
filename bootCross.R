library(parallel)
library(qgraph)
library(igraph)
library(colorspace)
#community detection and cluster overlay
overlayClusters<-function(graph, method=c("spinglass","walktrap"),label.cex=0.5){
  igraphQgRes<-NetworkToolbox::convert2igraph(abs(getWmat(graph)))
  if(method=="spinglass"){
    clustering<-cluster_spinglass(igraphQgRes,gamma=0.5,spins = 20)
    clusters<-clustering$membership
  }else{
    clustering<-cluster_walktrap(igraphQgRes)  
    clusters<-clustering$membership  
  }
  igraphQgRes<-as.igraph(graph,attributes = TRUE)
  clusters[is.na(clusters)]<-0
  clusters_numbers <- unique(clusters)
  n_clusters <- length(clusters_numbers)
  cols <- rainbow_hcl(n_clusters)
  Graph<-data.frame("id"=1:graph$graphAttributes$Graph$nNodes,
                    "node"=graph$graphAttributes$Nodes$labels)
  Graph$clusterID<-clusters
  if(length(which(Graph$clusterID==0)))Graph<-Graph[-which(Graph$clusterID==0),]
  for(i in 1:length(clusters_numbers)){
    Graph$borderColor[Graph$clusterID==clusters_numbers[[i]]]<- cols[[i]] 
  }
  valborderColor<-Graph$borderColor
  valborderWidth<-7
  group_ids<-list()
  for(i in 1:length(unique(Graph$clusterID))){
    group_ids[[i]]<-as.vector(Graph$id[Graph$clusterID==unique(Graph$clusterID)[[i]]])
    names(group_ids)[[i]]<-unique(Graph$clusterID)[[i]]
  }
  colsMatched<-vector()
  for(i in 1:n_clusters){
    
    colsMatched[[i]]<-cols[match(as.integer(names(group_ids)[i]),clusters_numbers)]
  }
  
  par(mar = rep(0.1, 4))
  group_color <- colsMatched
  # the fill gets an additional alpha value for transparency:
  group_color_fill <- paste0(group_color, '30')
  valqgresult2<-plot(igraphQgRes,
                     mark.groups = group_ids,
                     mark.col = group_color_fill,
                     mark.border = group_color,
                     layout=graph$layout,
                     mark.expand=15,
                     vertex.label.cex=label.cex)#input$icex)#,
  result<-list()
  result$Graph<-Graph
  return(Graph)
}
bootCross<-function(...,prefixes=c("S1","S2"),keyCol,nboots=1000,
                  ncores){
  
  if (missing(ncores)) {
    ncores <- parallel::detectCores()-1
  }else {
    ncores<-ncores
  }
  #sorting and merging cross-sections into "merged dataset"
    cross.sections<-list(...)
    cross.sections<-lapply(cross.sections,function(x){
      colnames(x)<-tolower(colnames(x))
      x
      })
    if (length(cross.sections)!=length(prefixes)){
      stop("number of cross.sectios must equal the number of prefixes")
    }else{
      commonColumns<-Reduce(intersect, lapply(cross.sections, names))
      for (i in (1:length(cross.sections))){
        cross.sections[[i]]<-cross.sections[[i]][,commonColumns]
        colnames(cross.sections[[i]])[-which(colnames(cross.sections[[i]])==keyCol)]<-paste0(prefixes[i],".", colnames(cross.sections[[i]])[-which(colnames(cross.sections[[i]])==keyCol)])
        
      }
    }
  data <- Reduce(function(...) merge(...,by=keyCol, all=T), cross.sections)
  
  data<-data[,-which(colnames(data)==keyCol)]
  for (i in (1:length(cross.sections))){
    cross.sections[[i]]<-cross.sections[[i]][,-which(colnames(cross.sections[[i]])==keyCol)]  
    }
  message("\nCross-sections Merged.")
  n<-nboots
  cases <- nrow(data)
  datalist <- list()
  #generating bootstrapped data
  message("\nBootstrapping ...", appendLF = FALSE)
  for(i in 1:n) {
    testdata<- data[sample(1:cases, replace = TRUE),]
    varCount<-apply(testdata,2,function(x){
      length(unique(na.omit(x)))
    })
    if(min(varCount,na.rm=TRUE)<2){
      next
    }else{
      datalist[[i]]<-testdata
    }
  }
  message("done", appendLF = TRUE)
  #averaging the bootstrapped data matrices to calculate node size
  nulls<-unlist(lapply(datalist,function(x)is.null(x)))
  if(length(which(nulls))){
    datalist<-datalist[-which(nulls)]
  }
  n<-length(datalist)
  arr.data <- array( unlist(datalist) , c(cases,ncol(data),n) )
  arr.data<-apply( arr.data , c(1,2) ,function(x){mean(x,na.rm=TRUE)}) #function(x){mean(x,na.rm = TRUE)/sd(x,na.rm = TRUE)} )
  nodesize<-as.vector(apply(arr.data,2,function(x){round(mean(x,na.rm=TRUE),2)}))#/sd(x,na.rm = TRUE)}))#function(x)mean(x,na.rm=TRUE)))
  #generating correlation matrices for bootstrapped data
  corlist <- list()
  message("\nComputing Correlation Matrices...\n", appendLF = FALSE)
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist","corlist", "cases"),
                          envir = environment())
  corlist <- pbapply::pblapply(X = datalist, cl = cl, FUN = qgraph::cor_auto)
  parallel::stopCluster(cl)
  message("done", appendLF = TRUE)
  #generating Representative Correlation Matrix(RCM)
  arr <- array( unlist(corlist) , c(ncol(data),ncol(data),n) )
  typical.core<-apply( arr , c(1,2) , function(x){median(x,na.rm = TRUE)} )
  #making RCM positive definite
  if(any(eigen(typical.core)$values < 0))
  {
    warning("RCM is not positive definite.\nFinding the nearest positive definite matrix using Matrix::nearPD()")
    typical.core<-as.matrix(Matrix::nearPD(typical.core, corr = TRUE, ensureSymmetry = TRUE)$mat)
  }
  #naming rows and columns of RCM
  colnames(typical.core)<-rownames(typical.core)<-colnames(data)
  #output
  result<-list()
  result$mergedCross<-data
  result$crossSections<-cross.sections
  result$longCore<-typical.core
  result$nodeSize<-nodesize
  result$bootDat<-datalist
  result$prefixes<-prefixes
  return(result)
}
#calling bootCross function
bootCrossed<-bootCross(crossSectionS1, crossSectionS3,prefixes=c("S1","S3"),keyCol="pid",nboots=1000)
#estimating the longitudinal network
groupsLong
#$RRB
#[1]  1  2  3  4  5  6  7 33 34 35 36 37 38 39

#$COM.NV
#[1]  8  9 10 11 12 13 14 40 41 42 43 44 45 46

#$COM.V
#[1] 15 16 17 47 48 49

#$SOC
#[1] 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 50 51 52 53 54 55 56 57 58 59
#[26] 60 61 62 63 64
groupColors<-c("#5CC96D","#9198FF","#DEB962","#FF9391")
longNodeSizes<-bootCrossed$nodeSize
longCores<-bootCrossed$longCore
graphLong<- qgraph::EBICglasso(longCores, nrow(bootCrossed$mergedCross),0.25)
qgraphLong<- qgraph::qgraph(graphLong,DoNotPlot=FALSE, legend=FALSE,
                                   arrows=F,rescale=TRUE,
                                   edge.labels=FALSE,
                                   fade=TRUE,
                                   edge.width=1,
                                   diag=FALSE,vTrans=250,
                                   label.cex=4.6,
                                   theme="colorblind",mar=rep(1,4),
                                   sampleSize=nrow(bootCrossed$mergedCross),
                                   vsize=longNodeSizes*5.3,
                                   aspect=TRUE,label.scale.equal=TRUE,
                                   plot=TRUE,groups=groupsLong,
                                   legend.mode = "style2",
                                   legend.cex = 0.32,
                                   color=as.character(groupColors),
                                   labels=colnames(graphLong),cut=0.2,
                                   layout="spring")

qgraphLongOverlay<-overlayClusters(qgraphLong,"spinglass")
#extracting stage correlation matrices and estimating cross-sectional stage sub-networks
groupsStage
#$RRB
#[1] 1 2 3 4 5 6 7

#$COM.NV
#[1]  8  9 10 11 12 13 14

#$COM.V
#[1] 15 16 17

#$SOC
#[1] 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32

graphStages<-list()
stageCores<-list()
qgraphStages<-list()
stageNodeSizes<-list()
dfSplitRowNames<-as.data.frame(strsplit(rownames(bootCrossed$longCore),".",fixed=TRUE))
for (i in (1:length(bootCrossed$prefixes))){
  stageColNums<-which(dfSplitRowNames[1,] ==bootCrossed$prefixes[[i]])
  stageNodeSizes[[i]]<-bootCrossed$nodeSize[stageColNums]
  stageCores[[i]]<-bootCrossed$longCore[stageColNums,stageColNums]
  colnames(stageCores[[i]])<-rownames(stageCores[[i]])<-substring(colnames(stageCores[[i]]),nchar(prefixes[[i]])+2)
  graphStages[[i]]<- qgraph::EBICglasso(stageCores[[i]], nrow(bootCrossed$mergedCross),0.5)
  qgraphStages[[i]]<- qgraph::qgraph(graphStages[[i]],DoNotPlot=FALSE , legend=FALSE,
                               arrows=F,rescale=TRUE,
                               edge.labels=FALSE,
                               fade=TRUE,
                               edge.width=1,
                               diag=FALSE,vTrans=250,
                               label.cex=4.6,
                               theme="colorblind",mar=rep(1,4),
                               sampleSize=nrow(bootCrossed$mergedCross),
                               vsize=stageNodeSizes[[i]]*5.3,
                               aspect=TRUE,label.scale.equal=TRUE,
                               plot=TRUE,groups=groupsStage,
                               legend.mode = "style2",
                               legend.cex = 0.32,
                               color=as.character(groupColors),
                               labels=colnames(graphStages[[i]]),cut=0.2,
                               layout="spring")
}
 layoutAvg<-averageLayout(qgraphStages)
 graphStagesAvg<-list()
 qgraphStageAvgOverlay<-list()
 for (i in (1:length(qgraphStages))){
  graphStagesAvg[[i]]<-qgraph(qgraphStages[[i]],layout=layoutAvg,DoNotPlot=TRUE)
  qgraphStageAvgOverlay[[i]]<-overlayClusters(graphStagesAvg[[i]],"spinglass")
  }
#Plotting Centrality Measures
plot.bootCross<-function (..., labels,
                     scale = c("z-scores", "raw", "raw0", 
                               "relative"),
                     include = c( "Strength", 
                                  "Closeness", "Betweenness", 
                                  "ExpectedInfluence"),
                     print = TRUE, verbose = TRUE,
                     standardized, relative, weighted = TRUE, 
                     signed = TRUE, orderBy="default" ,
                     decreasing = FALSE) 
 {
   if (any(include == "all") | any(include == "All")) {
     include <- c("Degree", "Strength", "OutDegree", "InDegree", 
                  "OutStrength", "InStrength", "Closeness", "Betweenness", 
                  "ExpectedInfluence", "OutExpectedInfluence", "InExpectedInfluence")
   }
   
   graphs<-list(...)
   if(all(unlist(lapply(graphs[[1]],function(x)class(x)))=="qgraph")){
     colorLabelsBy<-graphs[[1]][[1]]
   }else{
     stop("input must be a list of qgraph objects")
   }
   
   scale <- match.arg(scale)
   if (!missing(standardized)) {
     warning("'standardized' argument is deprecated and will be removed.")
   }
   else {
     standardized <- scale== "z-scores"
   }
   if (!missing(relative)) {
     warning("'relative' argument is deprecated and will be removed.")
   }
   else {
     relative <- scale == "relative"
   }
   if (scale == "z-scores") {
     if (verbose) 
       message("Note: z-scores are shown on x-axis rather than raw centrality indices.")
   }
   if (scale == "relative") {
     if (verbose) 
       message("Note: relative centrality indices are shown on x-axis rather than raw centrality indices.")
   }
   measure <- NULL
   value <- NULL
   node <- NULL
   type <- NULL
   Long <- centralityTable(..., standardized = standardized, 
                           labels = labels, relative = relative, weighted = weighted, 
                           signed = signed)
   Long <- subset(Long, measure %in% include)
   
   if (length(unique(Long$type))==2 & orderBy=="d") {
     Long$node<-factor(Long$node)
     wide = Long %>% 
       spread(measure, value)
     wide.scq<-wide
     r.scq<-nrow(wide.scq)/2
     wide.scq$d<-(wide.scq[(r.scq+1):nrow(wide.scq),6])-wide.scq[1:r.scq,6]
     wide.scq<- wide.scq[order(wide.scq$d),]
     wide.scq %>%
       mutate(node = fct_reorder(node, d))
     wide.scq$type<-factor(wide.scq$type)
     keycol <- "measure"
     valuecol <- "value"
     gathercols <- c(include,"d")
     Long<-gather_(wide.scq, keycol, valuecol, gathercols)
   }
   if (orderBy == "default") {
     nodeLevels <- unique(gtools::mixedsort(as.character( Long$node), 
                                           decreasing = decreasing))
   }else {
     nodeLevels <- names(sort(tapply(Long$value[Long$measure ==orderBy], Long$node[Long$measure == orderBy], mean),decreasing = decreasing))
   }
   Long$node <- factor(as.character(Long$node), levels = nodeLevels)
   Long <- Long[gtools::mixedorder(Long$node), ]
   wide.scq <- spread(Long, measure, value)
   
      match.color<-match(wide.scq$node,colorLabelsBy$graphAttributes$Nodes$labels)
   lblCol<-colorLabelsBy$graphAttributes$Nodes$color[match.color]
   lblCol<-lblCol[1:length(colorLabelsBy$graphAttributes$Nodes$color)]
   if (length(unique(Long$type)) > 1) {
     Long2<-Long[which(Long$measure != "d" ),]
     g <- ggplot(Long2, aes(x = factor(node), y = value, group = type, 
                            colour = type))
     g <- g + geom_line(size=1) + xlab("") + ylab("") + geom_point(size=2,fill="white",stroke=0.5,pch=24, colour="slateblue1")+
       facet_grid(measure ~ ., scales = "free", space = "free") +
       theme_minimal() + 
       theme(
         panel.grid.major.x = element_line(colour = "gray"),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.y = element_line(colour = "gray"),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x=element_text(angle=90, vjust = 0.5,colour = lblCol,face="bold"),
         strip.text.x = element_text(size = 2),
         strip.placement = "outside",
         strip.background = element_rect(fill="gray90", color = "lightblue2"), 
         panel.background = element_rect(fill = "lightcyan1",colour = "lightblue2"),
         legend.position = "top")
     g
   } 
}
plot.bootCross(qgraphStages,include =c("Strength","Closeness","Betweenness","ExpectedInfluence"),orderBy="Strength")
plot.bootCross(qgraphStages,include =c("Strength","Closeness","Betweenness","ExpectedInfluence"), orderBy="d")
