aggregateTaxa<-function(x,lineages,taxon.level="class", unknown=c(NA,"unclassified"), keepUnknown = FALSE, keepSumUnknown=FALSE, returnLineages=FALSE){
  # this is needed, so lineages.higher.x can be filled without error caused by levels
  if(is.data.frame(lineages)){
    lineages=as.matrix(lineages)
  }
  if(nrow(lineages)!=nrow(x)){
    stop("The lineage matrix should have as many rows as the abundance matrix!")
  }
  if(keepSumUnknown && keepUnknown){
    stop("Please choose either keepUnknown or keepSumUnknown to keep unknown taxa.")
  }
  taxon.level=tolower(taxon.level) # conversion to lower case
  unknown.lineage=c("Dummy_phylum","Dummy_kingdom","Dummy_class","Dummy_order","Dummy_family","Dummy_genus","Dummy_species")
  unknown=c(unknown,unknown.lineage)
  indicesUnknown=c()
  levelIndex=NA
  if(taxon.level=="species"){
    levelIndex=7
  }else if(taxon.level=="genus"){
    levelIndex=6
  }else if(taxon.level=="family"){
    levelIndex=5
  }else if(taxon.level=="order"){
    levelIndex=4
  }else if(taxon.level=="class"){
    levelIndex=3
  }else if(taxon.level=="phylum"){
    levelIndex=2
  }else if(taxon.level=="kingdom"){
    levelIndex=1
  }else{
    stop("Requested taxon level not supported")
  }
  levelMembers=unique(lineages[,levelIndex])
  if(!keepSumUnknown && !keepUnknown){
    # discard unclassified or dummy higher-level taxa
    levelMembers=setdiff(levelMembers,unknown)
    rowNames=levelMembers
    rowNum=length(levelMembers)
  }else{
    # keep unclassified higher-level taxa, but do not count them as rows in the result matrix
    temp=setdiff(levelMembers,unknown)
    rowNames=temp
    rowNum=length(temp)
  }
  print(paste("Number of higher-level taxa:",length(levelMembers)))
  higher.x=matrix(NA,nrow=rowNum,ncol=ncol(x))
  lineages.higher.x=matrix(NA,nrow=rowNum,ncol=7)
  print(dim(lineages.higher.x))
  rownames(higher.x)=rowNames
  colnames(higher.x)=colnames(x)
  rowCounter=1
  for(levelMember in levelMembers){
    #print(levelMember)
    # get indices of member taxa
    member.indices=which(lineages[,levelIndex]==levelMember)
    if(length(member.indices)>0){
      #print(paste("level index:",levelIndex))
      #print(paste("row counter:",rowCounter))
      #print("first member index:")
      #print(member.indices[1])
      #print(lineages[member.indices[1],1:levelIndex])
      # pick the lineage of the first member, but only up to the selected level
      lineages.higher.x[rowCounter,1:levelIndex]=lineages[member.indices[1],1:levelIndex]
      if((keepUnknown || keepSumUnknown) && (levelMember %in% unknown)){
        indicesUnknown=c(indicesUnknown,member.indices)
      }else{
        # sum member taxa in each column
        vec=c()
        if(length(member.indices)==1){
          vec=colSums(t(as.matrix(x[member.indices,])))
          #print(dim(t(as.matrix(x[member.indices,]))))
        }else{
          vec=colSums(x[member.indices,])
          #print(dim(x[member.indices,]))
        }
        higher.x[rowCounter,]=vec
        rowCounter=rowCounter+1
      }
      # deal with missing values in taxonomic lineages
    }else{
      if(is.na(levelMember) && (keepSumUnknown || keepUnknown)){
        #print(paste("Keeping taxa with assignment NA for level",taxon.level))
        indicesUnknown=c(indicesUnknown,which(is.na(lineages[,levelIndex])))
      }else{
        warning(paste("No members found for higher-level taxon",levelMember))
      }
    }
  }
  indicesUnknown=unique(indicesUnknown)
  if(keepSumUnknown){
    print(paste("Adding sum of",length(indicesUnknown),"taxa with unknown classification."))
    unknownSum=colSums(x[indicesUnknown,])
    higher.x=rbind(higher.x,"sumUnknown"=unknownSum)
    lineages.higher.x=rbind(lineages.higher.x,"sumUnknown"=unknown.lineage)
  }else if(keepUnknown){
    print(paste("Appending",length(indicesUnknown),"taxa with unknown classification."))
    unknown.taxa.matrix=x[indicesUnknown,]
    higher.x=rbind(higher.x,unknown.taxa.matrix)
    lineages.higher.x=rbind(lineages.higher.x,lineages[indicesUnknown,])
  }
  if(returnLineages){
    rownames(lineages.higher.x)=rownames(higher.x)
    res=list(higher.x,lineages.higher.x)
    names(res)=c("abundances","lineages")
    return(res)
  }else{
    return(higher.x)
  }
}

#' @title PCoA for microbial sequencing data
#'
#' @description A wrapper around various PCoA-based analyses implemented in vegan. The wrapper can handle groups and
#' metadata. PCoA is carried out sample-wise. The na.action is set to na.omit, however envfit cannot deal with
#' missing values, therefore if metadata are provided, they should be free of missing values.
#'
#' @details When a reference and groups are provided and the number of group memberships does not equal the number of samples in the combined abundance table,
#' groups are automatically extended such that reference samples are assigned to a single group with name refName, which is colored in gray.
#' The color vector is likewise extended if provided. If a clusters, time and/or labels vector is provided together with a referene, it has to refer to both data sets.
#' Samples in the abundance matrix are appended after the reference samples, so cluster memberships, time points and/or labels have to be provided in the same order.
#' Different total counts in abundances and reference samples may bias the result, therefore rarefyRef allows to rarefy both to the same total count after matching.
#' RarefyRef will rarefy all samples to the lowest total count found in any sample. Rows with zero counts after rarefaction are removed.
#' When topTaxa is set larger zero, significant top-varying taxa are shown. The permutation test is carried out by shuffling the selected number of top-covarying taxa.
#' Multiple testing correction on parameter-free p-values is then only applied to these top-covarying taxa. The strength of covariance is determined as
#' the norm of the vectors resulting from multiplying the standardized eigen vectors with the taxa.
#' In contrast, envfit p-values are computed for all metadata and multiple-testing correction is consequently applied to all metadata provided, though
#' only the selected number of most significant metadata are shown. Thus, topTaxa ranks taxa by covariance with significance as a filter, whereas
#' topMetadata ranks metadata by significance.
#' The number in the axis label brackets refers to the proportion of variance explained as computed with vegan's eigenvals function.
#' Note that ordination values in plots are not scaled (equivalent to plot.cca(scaling=0)).
#' If groups are provided, drawEllipse can be enabled to draw ellipses with vegan's ordiellipse and to run vegan's adonis (PERMANOVA) to assess whether group composition
#' differs significantly. Adonis is sensitive to differences in within-group variation. To assess whether within-group variation differs, betadisper can be carried out by
#' enabling groupDispersion.
#' In addition, if groups are provided, seqPCoA computes cluster quality indices to assess group separation with package clusterCrit.
#' Concerning interpretation of cluster quality indices: The silhouette index ranges between -1 and 1, a large Dunn's, Silhouette or Calinski-Harabasz index and a small
#' Davies-Bouldin or C-index indicate well-defined clusters, respectively.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param reference an optional reference data set on which abundances are mapped; data are merged by matching row names (nonmatching ones are kept as sum); cannot be combined with rda or topMetadata (topMetadata needs to be set to zero)
#' @param rarefyRef rarefy abundance and reference samples to the minimum total count found in any of the samples; recommended when the total counts differ
#' @param addToRefStepwise compute ordination coordinates for each sample added to the reference separately; cannot be combined with metadata, rda, errorbars or drawEllipse
#' @param refName group name for reference samples
#' @param metadata an optional data frame with metadata items as columns, where samples are in the same order as in abundances and data types (factor vs numeric) are supposed to be correct; if provided and rda is FALSE, envfit is carried out
#' @param groupAttrib optional: the name of a metadata item that refers to a vector that provides for each sample its group membership
#' @param groups an optional vector that provides for each sample its group membership in the same order as samples are given in abundances; overridden by groupAttrib
#' @param groupColors an optional map of predefined colors for groups that matches names in groups (which should be strings); adds a color legend to the plot; if reference is provided, refName is added if absent
#' @param colors an optional vector of colors with entries in the same order as samples to be colored; it overrides groupColors if provided
#' @param timeColors an optional vector of colors for arrows connecting samples interpreted as time points (only relevant when a time vector is supplied)
#' @param clusters an optional vector that provides for each sample its cluster membership (cluster membership is visualized through shape, up to 10 different shapes are possible, not supported for errorbars)
#' @param labels an optional vector that provides for each sample a label to display
#' @param sizes an optional vector that provides for each sample a numeric value that will be displayed as dot size (sizes will be shifted into positive range if necessary and scaled between 0.5 and 2.5)
#' @param size.legend a string displayed as a legend for size, FALSE suppresses size legend
#' @param errorbars an optional vector with as many entries as samples, defines groups for drawing error bars with ordibar (sd, conf set via ellipseConf); time connect centroids; color/size/group of first sample in error group assigned
#' @param errorbarLabels an optional vector with as many entries as groups in errorbars to label centroids; defaults to unique(errorbars); overrides labels (only applies when errorbars vector is supplied)
#' @param time an optional vector with as many time points as samples, adds arrows between consecutive time points (time points should be in ascending order; if groups are provided, time should be in ascending order per group)
#' @param hiddenTaxa an optional vector with names of taxa to be hidden (they will be taken into account for PCoA/RDA/envfit, but are not displayed among top co-varying taxa)
#' @param hiddenSamples an optional vector with indices of samples to be hidden (they will be taken into account for PCoA/RDA/envfit, but are not displayed)
#' @param dis dissimilarity or distance supported by vegan's vegdist function (if set to cor, a PCA is carried out using vegan's function rda with scale set to true)
#' @param rda carry out an RDA instead of a PCoA using vegan's capscale function
#' @param scale scale numeric metadata (subtract the mean and divide by standard deviation)
#' @param doScree do a Scree plot
#' @param topTaxa if larger than zero: show the top N taxa most strongly covarying with principal components as arrows in the PCoA if they are significant according to a permutation test
#' @param topMetadata if larger than zero, metadata provided and rda false: show the top N most significant numeric metadata as arrows and the top N most significant factor metadata as text in the PCoA
#' @param arrowFactor the length of taxon arrows (determined by scaled covariance) is multiplied with this factor
#' @param metadataFactor the length of numeric metadata arrows (determined by Pearson correlation) is multiplied with this factor
#' @param centroidFactor centroid positions (representing categoric metadata) are multiplied with this factor
#' @param taxonColor the color of the taxon arrows and text
#' @param metadataColor the color of the metadata arrows and text
#' @param hideGroupLegend do not show the group legend
#' @param drawEllipse if groups or groupAttrib given, draw polygons encapsulating groups using vegan's ordiellipse function (kind is sd, conf given via ellipseConf); print adonis R2 and p-value (permutation number given via env.permut)
#' @param ellipseOnClusters draw ellipse around cluster members (indicated by shape) instead of group members (indicated by color)
#' @param ellipseColorMap color ellipses according to this color list; entries have to be group names or, if ellipseOnClusters is true, cluster names; overrides groupColors for ellipses
#' @param clusterQualityIndex if groups or groupAttrib given, report cluster quality according to silhouette function or any criterium supported by package clusterCrit (default: silhouette, set to none to disable computation of cluster quality)
#' @param groupDispersion if groups or groupAttrib given, report Tukey's HSD test on differences in group dispersions (avg distance of group members to centroid) and do boxplot (wraps vegan's betadisper)
#' @param xlim range shown on the x axis, by default the minimum and maximum of the first selected component
#' @param ylim range shown on the y axis, by default the minimum and maximum of the second selected component
#' @param permut number of permutations for top-covarying taxa; if NA, NULL or smaller than 1, no permutation test is carried out
#' @param env.permut number of permutations for envfit, if drawEllipse is true, for adonis
#' @param pAdjMethod method for multiple testing correction supported by p.adjust for top-covarying taxon and envfit p-values
#' @param qvalThreshold threshold on multiple-testing corrected top-covarying taxon and envfit p-values
#' @param ellipseConf confidence limit for drawEllipse and errorbars
#' @param dimensions the principal components used for plotting, by default the first and second
#' @param verbose print out more information
#' @param \\dots Additional arguments passed to plot()
#' @examples
#' data("ibd_taxa")
#' data("ibd_metadata")
#' ibd_metadata=assignMetadataTypes(ibd_metadata,categoric=c("SRA_metagenome_name","Diagnosis"))
#' seqPCoA(ibd_taxa,groups=as.vector(ibd_metadata$Diagnosis),topTaxa=30, drawEllipse=TRUE)
#' # remove 65 samples with missing calprotectin measurements or other missing values in the metadata
#' na.indices=unique(which(is.na(ibd_metadata),arr.ind=TRUE)[,1])
#' indices.to.keep=setdiff(1:nrow(ibd_metadata),na.indices)
#' ibd_metadata=ibd_metadata[indices.to.keep,]
#' ibd_taxa=ibd_taxa[,indices.to.keep]
#' seqPCoA(ibd_taxa,metadata=ibd_metadata,groups=as.vector(ibd_metadata$Diagnosis),topTaxa=30)
#' @export
#'

seqPCoA<-function(abundances, reference=NULL, rarefyRef=FALSE, addToRefStepwise=FALSE, refName="ref", metadata=NULL, groupAttrib="", groups=c(), groupColors=NULL, colors=c(), timeColors=c(), clusters=c(), labels=c(), sizes=c(), size.legend="", errorbars = c(), errorbarLabels=unique(errorbars), time=c(), hiddenTaxa=c(), hiddenSamples=c(), dis="bray", rda=FALSE, scale=FALSE, doScree=FALSE, topTaxa=round(nrow(abundances)/2), topMetadata=round(nrow(metadata)/2), arrowFactor=0.5, metadataFactor=1, centroidFactor=1, taxonColor="brown", metadataColor="blue", hideGroupLegend=FALSE, drawEllipse=FALSE, ellipseOnClusters=FALSE, ellipseColorMap=NULL, clusterQualityIndex="silhouette", groupDispersion=FALSE, xlim=NULL, ylim=NULL, permut=1000, env.permut=1000, pAdjMethod="BH", qvalThreshold=0.05, ellipseConf=0.95, dimensions=c(1,2), verbose=FALSE, ...){
  
  # Test
  # path.vdp="/Users/u0097353/Documents/Documents_Karoline/MSysBio_Lab/Results/Nephrology/Data/vdp_genera.txt"
  # vdp=read.table(path.vdp,row.names=1,header=TRUE,sep="\t")
  # data("ibd_lineages")
  # ibd.genera=aggregateTaxa(ibd_taxa, lineages=ibd_lineages,taxon.level="genus")
  # ibd.genera.counts = round(ibd.genera*10000) # scale to counts (vdp already sums to 10K)
  # seqPCoA(ibd.genera.counts,reference=vdp, rarefyRef=TRUE, groups=groups, drawEllipse = TRUE, topMetadata = 0) # 53 matching genera
  # # beautified:
  # seqPCoA(ibd.genera.counts,reference=vdp, rarefyRef=TRUE, groups=groups, drawEllipse = TRUE, topMetadata = 0, xlim=c(-0.1,0.15), ylim=c(-0.1,0.1), arrowFactor=0.0001, clusterQualityIndex = "silhouette")
  # test step-wise
  # seqPCoA(ibd.genera.counts[,1:20],reference=vdp, rarefyRef=TRUE, groups=groups,xlim=c(-0.1,0.15), ylim=c(-0.1,0.1), addToRefStepwise = TRUE)
  
  ###### Internal constants ######
  env.locator=FALSE # label numeric env variables interactively (does not work with zoom)
  locator.width=0.35 # default is 0.25
  min.cex=0.5
  max.cex=2.5
  pch.value=16
  errorbarTransparency=0.3
  # 15=square, 16=circle, 17=triangle point up, 18=diamond, 25=triangle point down,
  # 3=plus sign, 4=multiplier sign, 8=star sign, 9=diamond with plus sign, 7=square with plus sign
  clus.pch.values=c(15,16,17,18,25,3,4,8,9,7)
  shift=0.05 # for taxa and numeric env data arrow labeling
  defaultColor="gray"
  
  ###### Initialization #########
  lastIndexRef=0
  cluster.colors=c()
  error.centroid.x=c()
  error.centroid.y=c()
  errorgroup.vs.group=list()
  errorgroupColors=c()
  errorgroupSizes=c()
  errorbarTransparentColors=c()
  errorgroupTimeColors=c()
  sample.coords=matrix(NA, nrow=ncol(abundances),ncol=2) # samples x axes
  display.size.legend=FALSE
  metadata.to.plot=c()
  #metadata.to.plot=c("IS_TOTAL","PCSG")
  
  ####### Checks #########
  if(length(clusters)>0){
    clusters=as.character(clusters)
  }
  
  if(!is.character(groups)){
    groups=as.character(groups)
  }
  
  if(length(errorbarLabels)>0 && length(errorbars)==0){
    stop("Please provide errorbar groups.")
  }
  
  if(rda && is.null(metadata)){
    stop("Metadata are needed for RDA!")
  }
  
  if(rda && !is.null(reference)){
    stop("RDA is not supported when a reference is provided.")
  }
  
  if(addToRefStepwise && is.null(reference)){
    stop("A reference is needed to add samples to reference step-wise.")
  }
  
  if(doScree && addToRefStepwise){
    stop("Scree plot is not available for step-wise addition of samples to reference.")
  }
  
  if(rda && addToRefStepwise){
    stop("RDA is not available for step-wise addition of samples to reference.")
  }
  
  # define default time vector colors
  if(length(timeColors)==0){
    timeColors=rep(defaultColor,ncol(abundances))
  }
  
  if(length(errorbars)>0){
    num.error.groups=length(unique(errorbars))
    if(length(errorbarLabels)>0 && length(errorbarLabels)!=num.error.groups){
      stop("Please provide as many error bar labels as there are error bar groups.")
    }
  }
  
  # make it easier for users
  if(addToRefStepwise){
    if(topTaxa>0){
      print("Biplot with taxa is not available for step-wise addition of samples to reference. TopTaxa is set to 0.")
    }
    topTaxa=0
    if(!is.null(metadata)){
      print("Metadata cannot be handled when samples are added step-wise to reference. Metadata are ignored.")
    }
    metadata=NULL
  }
  
  # if no metadata are provided, we cannot look for top-significant metadata items
  if(is.null(metadata)){
    topMetadata=0
  }
  
  # TODO: there is no reason why envfit is not supported in case metadata object provides data for both - to be modified
  if(topMetadata > 0 && !is.null(reference)){
    stop("Envfit is not supported when a reference is provided. Please set topMetadata to 0.")
  }
  
  ##### Preprocessing and more checks ##########
  
  if(!is.null(reference)){
    lastIndexRef=ncol(reference)
    # match abundances and reference by their row names
    res=intersectTables(reference,abundances,byRow = TRUE, keepSumNonMatched = TRUE)
    # append matched abundances to reference
    abundances=cbind(res$table1,res$table2)
    if(rarefyRef){
      # rarefy and discard taxa with zero abundance after rarefaction
      filtered=rarefyFilter(abundances)
      abundances=filtered$rar
      # update lastIndexRef to deal with columns removed during rarefaction
      ref.indices=intersect(1:lastIndexRef, filtered$colindices)
      lastIndexRef=length(ref.indices)
      # check
      print(paste("Name of last reference sample:",colnames(abundances))[lastIndexRef])
    }
    if(length(groups)>0){
      if(length(groups)!=ncol(abundances)){
        # extend groups to reference
        groups=c(rep(refName,ncol(res$table1)),groups)
      }
    }
    if(!is.null(groupColors)){
      if(!(refName %in% names(groupColors))){
        groupColors[[refName]]=defaultColor
      }
    }
    if(length(colors)>0){
      if(length(colors)!=ncol(abundances)){
        colors=c(rep(defaultColor,ncol(res$table1)),colors)
      }
    }
    if(length(clusters)>0){
      if(length(clusters)!=ncol(abundances)){
        stop("Please provide cluster memberships of reference and abundances combined, with reference first.")
      }
    }
    if(length(labels)>0){
      if(length(labels)!=ncol(abundances)){
        stop("Please provide labels of reference and abundances combined, with reference first.")
      }
    }
    if(length(time)>0){
      if(length(time)!=ncol(abundances)){
        stop("Please provide time points for reference and abundances combined, with reference first.")
      }
    }
  } # reference provided
  
  if(!is.null(metadata) && scale){
    # print("Scaling metadata")
    numeric.metadata=getMetadataSubset(metadata,type="numeric")
    catbin.metadata=getMetadataSubset(metadata,type="catbin")
    scaled.numeric.metadata=scale(numeric.metadata,center=TRUE, scale=TRUE)
    metadata=cbind(scaled.numeric.metadata,catbin.metadata)
  }
  
  if(verbose){
    print("Preparations done.")
  }
  
  ####### PCoA ###############
  
  if(rda){
    
    # it is not possible to give dissimilarities to capscale, so they need to be recomputed in case betadisper is enabled
    pcoa.res=capscale(data.frame(t(abundances))~.,metadata,distance=dis, na.action = "na.omit")
    
  }else if(addToRefStepwise){
    
    # loop over samples, carry out PCoA and extract coordinates of joint ordination with reference
    background=abundances[,1:lastIndexRef]
    samples=abundances[,(lastIndexRef+1):ncol(abundances)]
    print(paste("Carrying out",ncol(samples),"ordinations to match samples to reference one by one"))
    # loop samples to match to the background
    for(ord.index in 0:ncol(samples)){
      if(ord.index==0){
        matched=background
      }else{
        matched=cbind(background,samples[,ord.index])
      }
      if(dis=="cor"){
        # carry out standard PCA
        pcoa.res=rda(data.frame(t(matched)), scale=TRUE, na.action = "na.omit")
      }else{
        pcoa.res=capscale(data.frame(t(matched))~1,distance=dis, na.action = "na.omit")
      }
      # fill background sample positions from PCoA without extra sample
      if(ord.index==0){
        sample.coords[1:lastIndexRef,]=pcoa.res$CA$u[, dimensions]
      }else{
        print(paste("Completed ordination",ord.index))
        # extract coordinates of sample ordinated together with background
        sample.coords[(lastIndexRef+ord.index),]=pcoa.res$CA$u[(lastIndexRef+1), dimensions]
      }
    }
    # check
    #print(sample.coords[1:100,])
    
  }else{
    # carry out PCoA
    if(dis=="cor"){
      # carry out standard PCA
      # from vegan's doc: using correlation coefficients instead of covariances will give a more balanced ordination (scale=TRUE)
      print("Carrying out PCA...")
      pcoa.res=rda(data.frame(t(abundances)), scale=TRUE, na.action = "na.omit")
    }else{
      pcoa.res=capscale(data.frame(t(abundances))~1,distance=dis, na.action = "na.omit")
    }
    if(!is.null(metadata) && topMetadata>0){
      
      # carry out envfit
      # Lisa: >> betadisper computes eigenvectors that differ slightly from those in capscale
      # to be consistent, give the vector slot of the output of eigen to envfit
      # from betadisper code:
      # n <- attr(dis, "Size")
      # x <- matrix(0, ncol = n, nrow = n)
      # x[row(x) > col(x)] <- dis^2x <- x + t(x)
      # x <- dblcen(x)
      # e <- eigen(-x/2, symmetric = TRUE) <<
      # for the moment ignored, since betadisper/adonis is only used for assessment of cluster variability & quality and not for plotting
      ef=envfit(pcoa.res,metadata,perm=env.permut, choices=dimensions)
      #print(ef$vectors)
      #print(names(ef$vectors$r))
      if(length(metadata.to.plot)>0){
        for(metadatum.to.plot in metadata.to.plot){
          index.metadatum.to.plot=which(names(ef$vectors$r)==metadatum.to.plot)
          print(paste(sep="","R2 of ",metadatum.to.plot,": ",ef$vectors$r[index.metadatum.to.plot]))
          print(paste(sep="","P-value (not corrected) of ",metadatum.to.plot,": ",ef$vectors$pvals[index.metadatum.to.plot]))
        }
      }
      # correct for multiple testing using code in http://www.davidzeleny.net/anadat-r/doku.php/en:indirect_ordination_suppl
      pvals.vectors=p.adjust(ef$vectors$pvals, method=pAdjMethod)
      pvals.factors=p.adjust(ef$factors$pvals, method=pAdjMethod)
      # sort p-values in ascending orders, keep sorting result for other vector output
      pvals.vectors.sorted=sort(pvals.vectors, index.return=TRUE)
      pvals.vectors=pvals.vectors.sorted$x
      # factors make only use of names, but not of any other result
      pvals.factors=sort(pvals.factors)
      indices.vectors=which(pvals.vectors<qvalThreshold)
      indices.factors=which(pvals.factors<qvalThreshold)
      print(paste(length(indices.vectors),"significant numeric metadata found, in order of significance:"))
      #print(pvals.vectors)
      for(index in indices.vectors){
        #print(names(ef$vectors$r[index]))
        print(names(ef$vectors$r[pvals.vectors.sorted$ix[index]]))
      }
      print(paste(length(indices.factors),"significant categoric metadata found, in order of significance:"))
      for(index in indices.factors){
        print(names(pvals.factors)[index])
      }
      if(length(indices.vectors)>topMetadata){
        indices.vectors=indices.vectors[1:topMetadata]
      }
      if(length(indices.factors)>topMetadata){
        indices.factors=indices.factors[1:topMetadata]
      }
    } # do envfit
  } # not rda
  
  if(doScree){
    barplot(pcoa.res$CA$eig, names.arg=c(1:length(pcoa.res$CA$eig)), main="Scree plot", xlab="Eigenvalue index", ylab="Eigenvalue")
  }
  
  # assign the right labels
  redun.method="PCoA"
  if(dis=="cor"){
    redun.method="PCA"
  }
  
  if(!addToRefStepwise){
    # proportion of variance explained, see documentation of eigenvals function in vegan
    eig.sum=summary(eigenvals(pcoa.res))
    var.explained.1=eig.sum[2,dimensions[[1]]] # second row of eigenvalue summary: proportion of variance explained
    var.explained.2=eig.sum[2,dimensions[[2]]]
    xlab=paste(redun.method,dimensions[1]," [",round(var.explained.1,2),"]",sep="")
    ylab=paste(redun.method,dimensions[2]," [",round(var.explained.2,2),"]",sep="")
  }else{
    xlab=paste(redun.method,dimensions[1],sep="")
    ylab=paste(redun.method,dimensions[2],sep="")
  }
  
  # select indices of samples to show (if hiddenSamples is empty, this will be all)
  selected.sample.indices=setdiff(1:ncol(abundances),hiddenSamples)
  
  if(length(colors)>0){
    # use pre-assigned colors
    if(length(colors)!=ncol(abundances)){
      stop("There should be as many colors as samples in the color vector!")
    }
  }else{
    colors=rep(defaultColor,ncol(abundances))
    if(groupAttrib!=""){
      groups=metadata[[groupAttrib]]
      colors=assignColorsToGroups(groups,refName = refName, myColors = groupColors)
      if(length(errorbars)>0){
        errorbarTransparentColors=assignColorsToGroups(groups,refName = refName, myColors = groupColors, alpha=errorbarTransparency)
      }
    }else if(length(groups)>0){
      colors=assignColorsToGroups(groups, refName = refName, myColors = groupColors)
      if(length(errorbars)>0){
        errorbarTransparentColors=assignColorsToGroups(groups,refName = refName, myColors = groupColors, alpha=errorbarTransparency)
      }
    }
  }
  
  
  # assign cluster memberships as shapes
  if(length(clusters)>0){
    clusnum=length(unique(clusters))
    if(clusnum>length(clus.pch.values)){
      stop(paste("No more than",length(clus.pch.values),"clusters can be visualized."))
    }
    pch.value=c()
    clustersymbol.lookup=list()
    clus.values.counter=1
    for(clus.index in 1:length(clusters)){
      current.clus=clusters[clus.index]
      if(!(current.clus %in% names(clustersymbol.lookup))){
        clustersymbol.lookup[[current.clus]]=clus.pch.values[clus.values.counter]
        clus.values.counter=clus.values.counter+1
      }
    }
    for(clus.index in 1:length(clusters)){
      pch.value=c(pch.value,clustersymbol.lookup[[clusters[clus.index]]])
    }
    pch.value=pch.value[selected.sample.indices]
    #print(pch.value[20:50])
  }
  
  #print(selected.sample.indices)
  #print(colors)
  # determine range
  if(is.null(xlim)){
    xlim=range(pcoa.res$CA$u[selected.sample.indices,dimensions[1]])
  }
  if(is.null(ylim)){
    ylim=range(pcoa.res$CA$u[selected.sample.indices,dimensions[2]])
  }
  
  if(verbose){
    print("Size info:")
    print(sizes)
  }
  
  if(!is.null(sizes) || length(sizes)>0){
    display.size.legend=TRUE
    #print(paste("legend: ",size.legend))
    if(size.legend==FALSE){
      display.size.legend=FALSE
    }
    # shift into positive range
    if(min(sizes)<0){
      sizes=sizes-min(sizes)
    }
    # scale between 0.5 and 2.5
    sizes=(sizes-min(sizes))/(0.5*(max(sizes)-min(sizes)))
    sizes=sizes+min.cex # make sure there is no dot of size zero
    print(range(sizes))
  }else{
    #sizes=rep(par()$cex,length(selected.sample.indices))
    sizes=rep(par()$cex,ncol(abundances)) # default size
  }
  
  #print(length(sizes[selected.sample.indices]))
  #print(sizes[selected.sample.indices])
  
  if(length(errorbars)>0){
    # draw empty plot with the right dimensions and labels
    plot(pcoa.res$CA$u[selected.sample.indices,dimensions], type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    # calculate centroids but do not display ellipses
    res.ordi=ordiellipse(pcoa.res,groups = errorbars[selected.sample.indices],kind = "sd",draw="none", scaling=0)
    # assign group and color to error groups
    errorTemp=errorbarTransparentColors
    errorbarTransparentColors=c()
    for(i in 1:length(res.ordi)){
      error.sample.name=names(res.ordi)[i]
      error.sample.index=which(errorbars==error.sample.name)[1] # take first sample in error group to assign group and color
      if(length(groups)>0){
        group.of.error.sample=groups[error.sample.index]
        errorgroup.vs.group[[error.sample.name]]=group.of.error.sample
      }
      errorgroupColors=c(errorgroupColors,colors[error.sample.index])
      errorgroupSizes=c(errorgroupSizes,sizes[error.sample.index])
      errorbarTransparentColors=c(errorbarTransparentColors,errorTemp[error.sample.index])
      if(length(timeColors)>0){
        errorgroupTimeColors=c(errorgroupTimeColors,timeColors[error.sample.index])
      }
    }
    # draw centroids
    for(i in 1:length(res.ordi)){
      #print(paste("Error group size: ",errorgroupSizes[i]))
      points(x=res.ordi[[i]]$center[1], y=res.ordi[[i]]$center[2], col=errorgroupColors[i],cex=errorgroupSizes[i], pch=pch.value)
      error.centroid.x=c(error.centroid.x, res.ordi[[i]]$center[1])
      error.centroid.y=c(error.centroid.y, res.ordi[[i]]$center[2])
    }
    # add error bars
    #print(paste("Confidence interval for error bars",ellipseConf))
    # ordibar suddenly adds arrow heads...
    # length=0 should switch them off, but doesn't work. Cannot disable the plot to draw bars myself, plot=FALSE not an accepted graphical parameter. Could not pass on code=0.
    #res.ordibar=ordibar(pcoa.res,groups = errorbars[selected.sample.indices], length=0, kind = "sd", scaling=0, col=errorbarTransparentColors, conf=ellipseConf)
    # draw bars myself using ordibar output
    #for(i in 1:length(res.ordibar)){
    # print(paste("x0=",res.ordibar[[i]]$center[1]))
    # print(paste("y0=",res.ordibar[[i]]$center[2]))
    # print((res.ordibar[[i]]$center[1]+res.ordibar[[i]]$cov[1,1]))
    # print((res.ordibar[[i]]$center[2]+res.ordibar[[i]]$cov[1,2]))
    #arrows(x0=res.ordibar[[i]]$center[1], y0=res.ordibar[[i]]$center[2], x1=(res.ordibar[[i]]$center[1]+res.ordibar[[i]]$cov[1,1]), y1=(res.ordibar[[i]]$center[2]+res.ordibar[[i]]$cov[1,2]), col=errorbarTransparentColors[i], length=0.1, lty=2, code=0) # code = 0: no arrow head
    #}
    # ordiellipse as best alternative solution (respects conf)
    ordiellipse(pcoa.res,groups = errorbars[selected.sample.indices],kind = "sd", scaling=0, col=errorbarTransparentColors, conf=ellipseConf, draw="lines")
    #ordispider(pcoa.res,groups = errorbars[selected.sample.indices],scaling=0, col=errorbarTransparentColors)
  }else{
    if(!addToRefStepwise){
      plot(pcoa.res$CA$u[selected.sample.indices,dimensions],cex=sizes[selected.sample.indices], col=colors[selected.sample.indices], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, bg=colors[selected.sample.indices], pch=pch.value, ...)
    }else{
      plot(sample.coords[selected.sample.indices,dimensions], cex=sizes[selected.sample.indices], col=colors[selected.sample.indices], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, bg=colors[selected.sample.indices], pch=pch.value, ...)
    }
  }
  
  # add sample labels if requested
  if(length(labels)>0 && length(errorbarLabels)==0){
    #print(paste("label num:", length(labels[selected.sample.indices])))
    #print(labels)
    text(pcoa.res$CA$u[selected.sample.indices,dimensions[1]],pcoa.res$CA$u[selected.sample.indices,dimensions[2]],labels=labels[selected.sample.indices], pos=3, cex=0.9)
    if(verbose){
      print("Labels added.")
    }
  }
  if(length(errorbarLabels)>0){
    text(error.centroid.x,error.centroid.y,labels=errorbarLabels, pos=3, cex=0.9)
    if(verbose){
      print("Error bar labels added.")
    }
  }
  
  # add arrows between consecutive time points if requested
  if(length(time)>0){
    if(length(errorbars)>0){
      errorArrowColors=rep(defaultColor,length(unique(errorbars)))
      if(length(timeColors)>0){
        errorArrowColors=errorgroupTimeColors
      }
      for(i in 1:(length(res.ordi)-1)){
        arrows(x0=res.ordi[[i]]$center[1], y0=res.ordi[[i]]$center[2], x1=res.ordi[[i+1]]$center[1], y1=res.ordi[[i+1]]$center[2], col=errorArrowColors, length=0.1, lty=2) # angle=20
        #print(paste("arrow start coordinates",res.ordi[[i]]$center[1], res.ordi[[i]]$center[2]))
        #print(paste("arrow end coordinates",res.ordi[[i+1]]$center[1], res.ordi[[i+1]]$center[2]))
      }
    }else{
      for(i in 1:(nrow(pcoa.res$CA$u)-1)){
        # both samples to be connected are visible
        if(!(i %in% hiddenSamples) && !((i+1) %in% hiddenSamples)){
          if(length(groups)==0 || groups[i]==groups[i+1]){
            arrows(x0=pcoa.res$CA$u[i,dimensions[1]],y0=pcoa.res$CA$u[i,dimensions[2]],x1=pcoa.res$CA$u[i+1,dimensions[1]],y1=pcoa.res$CA$u[i+1,dimensions[2]], length=0.1, col=timeColors[i], lty=2, angle=20) # length=length of the edges of the arrow head in inches, lty=2 is dashed, angle refers to shaft vs edge of arrow head
          } # within the same group
        } # visible
      } # loop samples
    } # no errorbars
    
    if(verbose){
      print(paste("Number of time points:",length(time)))
      print(paste("Number of hidden samples:",length(hiddenSamples)))
      print("Arrows added.")
    }
  } # time provided
  
  
  groups.factor=factor(groups,levels=unique(as.character(groups)))
  
  if(verbose){
    #print("group factors:")
    #print(groups.factor)
  }
  
  # add ellipses if requested
  if(drawEllipse==TRUE && !addToRefStepwise){
    if(length(groups)==0){
      warning("Please specify groups in order to draw ellipses!")
    }else{
      show.groups=unique(groups[selected.sample.indices])
      #print(show.groups)
      # draw transparent ellipses around samples of the same groups with given confidence
      # by default, color order is not correct, because ordiellipse calls factor on the groups, which sorts entries alphabetically
      # for this reason, factor with correct ordering is reassigned to groups
      # color ellipse
      #ordiellipse(pcoa.res,scaling=0,groups=groups.factor,draw="polygon",col=unique(colors),alpha=0.25,conf=ellipseConf,lwd=0.5, kind="sd", border=0)
      # color border rather than ellipse itself
      if(ellipseOnClusters && length(clusters)>0){
        show.groups=unique(clusters[selected.sample.indices])
        cluster.factor=factor(clusters,levels=unique(clusters))
        if(!is.null(ellipseColorMap)){
          cluster.colors=assignColorsToGroups(clusters, refName = refName, myColors = ellipseColorMap)
          #print(unique(cluster.colors))
        }else{
          cluster.colors=assignColorsToGroups(clusters, refName = refName)
        }
        ordiellipse(pcoa.res,scaling=0,groups=cluster.factor, show.groups = show.groups ,draw="polygon",alpha=0.25,conf=ellipseConf,lwd=0.5, kind="sd", border=unique(cluster.colors))
      }else{
        ordiellipse(pcoa.res,scaling=0,groups=groups.factor, show.groups = show.groups, draw="polygon",alpha=0.25,conf=ellipseConf,lwd=0.5, kind="sd", border=unique(colors))
      }
      # adonis fails for step-wise addition of samples to reference
      adonis_results = adonis(data.frame(t(abundances)) ~ groups.factor, permutations = env.permut, method=dis)
      # variance explained through groups
      adonis.r2=adonis_results$aov.tab[1,5]
      adonis.pval=adonis_results$aov.tab[1,6]
      print("Adonis to test for significant difference in group compositions")
      print(paste("Adonis R2: ",round(adonis.r2,4),", p-value: ",round(adonis.pval,4),sep=""))
    }
  }
  
  # report cluster quality of groups
  if(length(groups)>0){
    if(clusterQualityIndex != "" && clusterQualityIndex != "none"){
      if(clusterQualityIndex=="CH"){
        clusterQualityIndex="Calinski_Harabasz"
      }
      print(paste("Cluster quality index", clusterQualityIndex))
      #print(is.numeric(abundances))
      if(clusterQualityIndex=="silhouette"){
        print(silhouette(abundances,groups=groups,method=dis))
      }else{
        # note that clusterCrit gives an error for large input matrices, but not for small ones; this error is not yet fixed
        print(intCriteria(t(abundances), part=groupsToNumeric(groups), crit=clusterQualityIndex)[[1]])
      }
    }
  }
  
  if(topTaxa>0){
    # taken from biplot.pcoa in ape
    # columns are taxa
    Y=t(abundances)
    n <- nrow(Y) # n = sample number
    # standardize eigen vectors (subtract mean and divide by standard deviation)
    # eigen vector dimensions: samples x selected dimensions
    ev.stand <- scale(pcoa.res$CA$u[,dimensions])
    # covariance between taxa (as columns) and selected principal components (standardized eigen vectors)
    S <- cov(Y, ev.stand)
    # scale S by the eigen values
    U <- S %*% diag((pcoa.res$CA$eig[dimensions]/(n-1))^(-0.5))
    colnames(U) <- colnames(pcoa.res$CA$u[,dimensions])
    # U dimensions: taxa x selected dimensions
    
    # select top covarying taxa in U
    # arrrow length codes strength of covariance with eigen vectors
    norms=apply(U,1,myNorm)
    sorted=sort(norms,index.return=TRUE,decreasing=TRUE)
    sorted.top.indices=sorted$ix[1:topTaxa]
    pvalues=c()
    U.sub=U[sorted.top.indices,]
    Y.sub=Y[,sorted.top.indices]
    was.permuted=FALSE
    # if requested, carry out permutations
    if(!is.null(permut) && !is.na(permut) && permut>0){
      was.permuted=TRUE
      # taxa vs permutations
      permuted.norms=matrix(NA,nrow=ncol(Y.sub),ncol=permut)
      # carry out permutation test
      for(iteration in 1:permut){
        # shuffle abundances separately per column (in Y, taxa are columns)
        Y.rand=apply(Y.sub,2,base::sample)
        # only look at top covarying taxa
        S.rand=cov(Y.rand,ev.stand)
        # scale the shuffled S by eigen values
        U.rand <- S.rand %*% diag((pcoa.res$CA$eig[dimensions]/(n-1))^(-0.5))
        # taxon arrows are defined by U, which has as many rows as taxa and as many columns as selected eigen vectors
        # compute arrow length as vector norm for each row (each row represents one arrow)
        permuted.norms[,iteration]=apply(U.rand,1,myNorm)
      }
      # compute signficance of vector norms
      for(taxon.index in 1:ncol(Y.sub)){
        obs.norm=norms[taxon.index]
        # check how many permuted norms are larger than observed norm
        r=length(which(permuted.norms[taxon.index,]>obs.norm))
        # compute parameter-free p-value
        pvalues=c(pvalues,(r+1)/(permut+1))
      }
      # adjust p-values for multiple testing and discard corrected p-values below selected significance level
      pvalues=p.adjust(pvalues,method=pAdjMethod)
      sig.pvalue.indices=which(pvalues<qvalThreshold)
      # only keep significant top covarying taxa
      U.selected=U.sub[sig.pvalue.indices,]
    } # end permutation test
    else{
      U.selected=U.sub
      sig.pvalue.indices=c()
    }
    if(length(sig.pvalue.indices)>0 || !was.permuted){
      if(was.permuted){
        print(paste("Among the top ",topTaxa," covarying taxa, ",length(sig.pvalue.indices)," are significant.",sep=""))
        for(sig.taxon.name in rownames(U.selected)){
          print(sig.taxon.name)
        }
      }
      if(length(hiddenTaxa)>0){
        visibleTaxa.indices=c()
        # loop selected taxa
        for(U.taxon.index in 1:nrow(U.selected)){
          if(!(rownames(U.selected)[U.taxon.index] %in% hiddenTaxa)){
            visibleTaxa.indices=c(visibleTaxa.indices,U.taxon.index)
          } # check whether selected taxon should be hidden
        }
        U.selected = U.selected[visibleTaxa.indices,]
      }
      arrows(0, 0, U.selected[, 1] * arrowFactor, U.selected[, 2] * arrowFactor, col = taxonColor,length = 0.1, lty=2)
      
      for(i in 1:nrow(U.selected)){
        for(j in 1:ncol(U.selected)){
          if(U.selected[i,j]>0){
            U.selected[i,j]=U.selected[i,j]+shift
          }else if(U.selected[i,j]<0){
            U.selected[i,j]=U.selected[i,j]-shift
          }
        }
      }
      text(U.selected*arrowFactor, rownames(U.selected), cex = 0.9, col = taxonColor)
    }else{
      if(was.permuted){
        warning("None of the top-covarying taxa is significant.")
      }
    }
  }
  
  if(!is.null(metadata) && topMetadata>0){
    if(length(indices.vectors)>0 || length(metadata.to.plot)>0){
      if(length(metadata.to.plot) > 0){
        #print(names(ef$vectors$r[pvals.vectors.sorted$ix]))
        for(metadatum.to.plot in metadata.to.plot){
          index.metadatum.to.plot=which(names(ef$vectors$r[pvals.vectors.sorted$ix])==metadatum.to.plot)
          indices.vectors=c(indices.vectors,index.metadatum.to.plot)
        }
      }
      # add arrows for numeric metadata
      ef.arrows=as.matrix(ef$vectors$arrows[pvals.vectors.sorted$ix[indices.vectors],dimensions])
      if(length(indices.vectors)==1){
        ef.arrows=t(ef.arrows)
      }
      lengths=ef$vectors$r[pvals.vectors.sorted$ix[indices.vectors]]
      ef.names=names(ef$vectors$r[pvals.vectors.sorted$ix[indices.vectors]])
      arrows(0, 0, ef.arrows[, 1]*lengths*metadataFactor, ef.arrows[, 2]*lengths*metadataFactor, col = metadataColor, length=0.1)
      
      # move metadata labels to nicer places
      for(i in 1:nrow(ef.arrows)){
        for(j in 1:ncol(ef.arrows)){
          if(ef.arrows[i,j]>0){
            ef.arrows[i,j]=ef.arrows[i,j]+shift
          }else if(ef.arrows[i,j]<0){
            ef.arrows[i,j]=ef.arrows[i,j]-shift
          }
        }
      }
      # pos=1, las = 1 means horizontal but is not effective, cex=0.8
      if(env.locator){
        # tolerance defaults to 0.25
        identify(ef.arrows*lengths*metadataFactor, tolerance=locator.width, labels=ef.names, cex = 0.8, col = metadataColor)
      }else{
        text(ef.arrows*lengths*metadataFactor, ef.names, cex = 0.8, col = metadataColor)
      }
    }
    if(length(indices.factors)>0){
      # the order of factors was altered, but since we match by name, this is not problematic
      sig.factors=names(indices.factors)
      # find centroids belonging to significant factors
      # centroid: mean or median dissimilarity of samples belonging to given level of a factor in PCoA
      # significance of centroid separation: now assessed with adonis
      for(sig.factor in sig.factors){
        sig.factor.centroid.indices=which(ef$factors$var.id==sig.factor)
        # several values of the categoric variable can be significant
        for(sig.factor.centroid.index in sig.factor.centroid.indices){
          sig.factor.value=rownames(ef$factors$centroids)[sig.factor.centroid.index]
          centroid.pos=ef$factors$centroids[sig.factor.centroid.index,]
          text(x=centroid.pos[1]*centroidFactor,y=centroid.pos[2]*centroidFactor,sig.factor.value, cex=0.8,col=metadataColor)
        } # end loop values of significant factors
      } # end loop significant factors
    } # significant factors found
  }
  
  if(length(groups)>0 && !hideGroupLegend){
    # cex=0.9
    legend("topright",legend=unique(groups[selected.sample.indices]),cex=0.9, pch = rep("*",length(unique(groups[selected.sample.indices]))), col = unique(colors[selected.sample.indices]), bg = "white", text.col="black")
  }
  if(length(clusters)>0){
    if(ellipseOnClusters){
      legend("topleft",legend=unique(clusters[selected.sample.indices]),cex=0.9, pch = unique(pch.value), col = "black", bg = "white", text.col=unique(cluster.colors))
    }else{
      legend("topleft",legend=unique(clusters[selected.sample.indices]),cex=0.9, pch = unique(pch.value), col = "black", bg = "white", text.col="black")
    }
  }
  # 1 = default size
  if(display.size.legend==TRUE){
    min.legend="min"
    max.legend="max"
    if(size.legend!=""){
      min.legend=paste(min.legend," ",size.legend,sep="")
      max.legend=paste(max.legend," ",size.legend,sep="")
    }
    legend("bottomleft",legend=c(min.legend,max.legend),pt.cex=c(min.cex,max.cex), pch = c(1,1), col = "black", bg = "white", text.col="black")
  }
  
  # color map legend
  if(!is.null(groupColors)){
    legend.colors=c()
    for(name in names(groupColors)){
      legend.colors=c(legend.colors,groupColors[[name]])
    }
    # invisible
    trans.col= rgb(255, 255, 255, maxColorValue = 255, alpha = 125)
    legend("bottomright",legend=names(groupColors), box.lty=0, box.col=trans.col, bg=trans.col, cex=0.9, text.col=legend.colors)
  }
  
  # A non-significant result in betadisper is not necessarily related to a significant/non-significant result in adonis.
  # Betadisper tests homogeneity of dispersion among groups, which is a condition (assumption) for adonis.
  # Betadisper can be done to see if one group has more compositional variance than another.
  # Adonis tests whether composition among groups is similar or not.
  if(groupDispersion){
    beta.out=betadisper(vegdist(t(abundances),method=dis), groups.factor, type='centroid')
    print("Tukey's HSD to test for significant differences in group variance (betadisper)")
    t.hsd=TukeyHSD(beta.out)
    print(t.hsd$group)
    # fourth column: p-values
    if(length(which(t.hsd$group[,4]<0.05))){
      print("Variance differs signficantly for at least 2 groups, so adonis results may be biased!")
    }
    boxplot(beta.out,xlab="Group")
  }
  
}

#' @title Assign colors to groups
#'
#' @description Given a group membership vector, each group receives
#' a unique color such that the resulting color vector is as long as
#' the group membership vector.
#'
#' @param groups a vector of group memberships
#' @param refName the name of the reference group; receives gray as a color
#' @param myColors an optional map of predefined colors for groups
#' @param returnMap whether to return the color map together with the colors
#' @param hiddenSamples optional vector of indices of samples to be hidden; these are assigned an invisible color
#' @param alpha transparency level, 0 for most transparent and 1 for most opaque
#' @return a vector of colors; if returnMap is true, a vector of colors and the color map
#' @export
#'
assignColorsToGroups<-function(groups, refName="ref", myColors = NULL, returnMap=FALSE, hiddenSamples=c(), alpha=1){
  # make a fully transparent color to hide samples
  trans.col=rgb(0, 1, 0, alpha=0)
  
  refContained=FALSE
  if(refName %in% groups){
    refContained=TRUE
  }
  groupNum=length(unique(groups))
  if(refContained){
    groupNum=groupNum-1 # reference counts extra
  }
  col.vec = seq(0,1,1/groupNum)
  hues = hsv(col.vec,alpha=alpha)
  #print(hues)
  hueCounter=1
  colors=c()
  # fill the color map
  if(is.null(myColors)){
    myColors=list()
    for(group.index in 1:length(groups)){
      group=as.character(groups[group.index])
      if(!(group %in% names(myColors))){
        if(group==refName){
          myColors[[group]]="gray"
        }else{
          myColors[[group]]=hues[hueCounter]
          hueCounter=hueCounter+1
        }
      }
    } # loop samples
  } # color map already filled
  
  #print(myColors)
  # assign colors from the color map
  for(group.index in 1:length(groups)){
    #print(groups[group.index])
    #print(myColors[[as.character(groups[group.index])]])
    if(group.index %in% hiddenSamples){
      colors=c(colors,trans.col)
    }else{
      colors=c(colors,myColors[[as.character(groups[group.index])]])
    }
  }
  
  if(returnMap){
    res=list(colors,myColors)
    names(res)=c("colors","colorMap")
    return(res)
  }
  return(colors)
}

# convert a group vector with strings in a numeric vector of cluster membership integers
groupsToNumeric<-function(groups=c()){
  clus.mem=groups
  groups.unique=unique(groups)
  for(i in 1:length(groups.unique)){
    current.group=groups.unique[i]
    group.indices=which(groups==current.group)
    clus.mem[group.indices]=i
  }
  return(as.integer(as.numeric(clus.mem)))
}


# compute the norm of a vector
myNorm<-function(x){
  if(!is.numeric(x)){
    stop("x should be a numeric vector.")
  }
  return(sqrt(sum(x^2)))
}

# Rarefaction combined with sample filtering
#
# Rarefy a matrix to the given minimum count number column-wise
# using vegan's rrarefy function. If columns have less than the minimum count number,
# they are discarded. Rows that have a sum of zero after rarefaction are also discarded.
# x a matrix
# min minimum count to which x is to be rarefied (if equal to zero, the minimum column sum is taken as min)
# a list with the rarefied matrix (rar) and the indices of the rows and the columns that were kept (rowindices and colindices)
rarefyFilter<-function(x,min = 0){
  keep=c()
  if(min < 0){
    stop("Min should be either 0 or positive.")
  }
  if(min == 0){
    min=min(colsums=apply(x,2,sum))
    print(paste("Rarefy to minimum count",min))
    keep=c(1:ncol(x))
  }else{
    colsums=apply(x,2,sum)
    # there are columns below the minimum
    if(min(colsums) < min){
      # loop column sums
      for(j in 1:ncol(x)){
        if(colsums[j] >= min){
          keep=c(keep,j)
        }
      }
      print(paste("Number of columns",ncol(x)))
      print(paste("Keeping ",length(keep)," columns with column sums equal or above",min))
      x=x[,keep]
    }
  }
  rar=t(vegan::rrarefy(t(x),min))
  zero.indices=which(rowSums(rar)==0)
  # discard taxa with zero sums
  if(length(zero.indices)>0){
    keep.nonzero=setdiff(1:nrow(rar),zero.indices)
    rar=rar[keep.nonzero,]
  }
  res=list(rar,keep.nonzero,keep)
  names(res)=c("rar","rowindices","colindices")
  return(res)
}

#' @title Convert metadata into numeric form
#'
#' @description Binary metadata items are converted into binary numeric metadata items (with 0/1).
#' Categoric metadata with more than two categories can be binarized, such that
#' each category is represented by a separate binary metadata item. Numeric metadata are kept as is.
#' Dates, when specified, are converted into days since the reference date.
#' Note that constant metadata are removed.
#'
#' @param metadata a dataframe
#' @param yes the symbol used for the first value in a binary metadata item (e.g. "Y")
#' @param no the symbol used for the second value in a binary metadata item (e.g. "N")
#' @param na.threshold remove metadata with more than the maximum allowed percentage of missing values
#' @param to.skip names of metadata to skip from binarization
#' @param binarize convert categoric metadata items into as many binary metadata items as there are categories (if false, metadata with more than 2 categories are removed)
#' @param date.items names of metadata items that represent dates
#' @param format the format used for the date items (an example date fitting the default format is 26/1/80)
#' @param referenceDate reference date used for conversion of dates into numbers (the reference date format is always d/m/Y)
#' @param remove.neg remove metadata items with negative values
#' @return a purely numeric dataframe
#' @export
metadataToNumeric<-function(metadata, yes="Y", no="N", na.threshold=100, to.skip=c(), binarize=TRUE, date.items=c(),format="%d/%m/%y", referenceDate="1/1/1900", remove.neg=TRUE){
  binarizedMetadata=list()
  metadata.to.remove=c()
  skip.from.conversion=c()
  # loop metadata items
  for(name in names(metadata)){
    #print(paste("Processing",name))
    # cannot remove leading or trailing white spaces, since this alters the content of factors
    # convert date into numeric
    if(name %in% date.items){
      print(paste("Processing date item",name))
      metadata[[name]]=as.numeric(as.Date(metadata[[name]],format=format)-as.Date(referenceDate,format="%d/%m/%Y"))
    }else{
      if(is.factor(metadata[[name]])){
        levels=levels(metadata[[name]])
        levels=setdiff(levels,NA) # remove NA
        # process binary metadata
        if(length(levels)==2){
          skip.from.conversion=c(skip.from.conversion,name)
          numLevels=as.numeric(levels) # gives warnings (invalid factor level, NA generated) in case levels are non-numeric
          if(0 %in% numLevels && 1 %in% numLevels){
            print(paste("Binary metadata item",name,"is already numeric."))
            metadata[[name]]=as.character(metadata[[name]])
          }else{
            metadata[[name]]=as.character(metadata[[name]])
            yes.indices=which(metadata[[name]]==yes)
            no.indices=which(metadata[[name]]==no)
            # convert into a numeric metadata item, keep NA values
            replacement=rep(NA,length(metadata[[name]]))
            replacement[yes.indices]=1
            replacement[no.indices]=0
            metadata[[name]]=replacement
          }
        }else if(length(levels)>2){
          #print("Categoric factor")
          metadata.to.remove=c(metadata.to.remove,name)
          # binarize categoric metadata
          if(binarize && !(name %in% to.skip)){
            print(paste("Binarizing metadata", name))
            #print(levels)
            # there are less categories than samples
            if(length(levels)<nrow(metadata)){
              categoryNames=c()
              # initialize category-specific metadata items as absent
              # levels are NA-free
              for(level in levels){
                categoryName=paste(name,trimws(level),sep="-")
                categoryNames=c(categoryNames,categoryName)
                binarizedMetadata[[categoryName]]=rep(0,nrow(metadata))
              }
              # set presence values
              for(level in levels){
                for(j in 1:nrow(metadata)){
                  value=metadata[[name]][j]
                  if(!is.na(value)){
                    categoryName=paste(name,trimws(value),sep="-")
                    binarizedMetadata[[categoryName]][j]=1
                  }else{
                    # set all categories to NA
                    for(categoryName in categoryNames){
                      binarizedMetadata[[categoryName]][j]=NA
                    }
                  }
                }
              }
              #print(binarizedMetadata)
            }else{
              warning(paste("Cannot binarize categoric metadata item",name,"because it has as many categories as samples."))
            }
          }
        } # only 1 category - will be removed below
      } # end factor metadata
    } # not a date
  } # loop metadata
  # append binarized metadata
  if(binarize && length(binarizedMetadata)>0){
    for(append in names(binarizedMetadata)){
      #print(append)
      #print(length(binarizedMetadata[[append]]))
      metadata[[append]]=binarizedMetadata[[append]]
    }
  }
  metadata.filtered=filterMetadata(metadata,toFilter = metadata.to.remove,na.threshold = na.threshold, remove.neg=remove.neg)
  return(metadata.filtered)
}

#' @title Filter metadata
#'
#' @description Remove metadata with more than the given percentage
#' of missing values, constant metadata and metadata with names as provided
#' in the toFilter vector.
#'
#' @param data a data frame with metadata items
#' @param toFilter names of metadata to remove
#' @param na.threshold percentage of NA values allowed, between 0 and 100
#' @param remove.neg remove metadata items with negative values
#' @return a data frame of filtered metadata items
#' @export
filterMetadata<-function(data, toFilter=c(), na.threshold=0, remove.neg=FALSE){
  if(!is.data.frame(data)){
    stop("Please provide a data frame.")
  }
  onePerc=nrow(data)/100
  na.number.threshold=round(onePerc*na.threshold)
  print(paste("Allowed number of NAs:",na.number.threshold))
  indices.tokeep=c()
  index.counter=1
  # loop metadadata items
  for(name in names(data)){
    values=unique(data[[name]])
    values=setdiff(values,NA)
    numberNA=length(which(is.na(data[[name]])))
    #print(paste(name,":",numberNA))
    if(length(values)>1 && numberNA<=na.number.threshold){
      if(name %in% toFilter){
        print(paste("Skipping metadata item: ",name,sep=""))
      }else{
        if(!is.factor(data[[name]]) && remove.neg){
          if(length(which(as.numeric(values)<0))==0){
            indices.tokeep=c(indices.tokeep,index.counter)
          }else{
            print(paste("Skipping metadata item with negative values: ",name,sep=""))
          }
        }else{
          indices.tokeep=c(indices.tokeep,index.counter)
        }
      }
    }else{
      if(numberNA>na.number.threshold){
        print(paste("Filtering metadata",name,"with index",index.counter,"because it has more missing values (namely",numberNA,") than the allowed percentage"))
      }else{
        print(paste("Filtering metadata",name,"with index",index.counter,"because it is constant"))
      }
    }
    index.counter=index.counter+1
  } # end loop
  return(data[,indices.tokeep])
}

#' @title Assign data types to metadata
#'
#' @description Treat all metadata with more than two values
#' as numeric, unless they are in the vector with categoric values.
#' Metadata with two values (without counting missing values) are treated as
#' binary (so categoric).
#' Metadata are supposed to have been filtered previously (so constant rows were removed).
#' Note that a message is printed to indicate which metadata type was assigned and to
#' warn about numeric metadata with non-numeric values.
#'
#' @param data a metadata matrix with rows as samples and columns as items or a dataframe with metadata items as columns
#' @param categoric metadata items to be treated as categoric data
#' @export
assignMetadataTypes<-function(data, categoric=c()){
  if(is.data.frame(data)==FALSE){
    data=as.data.frame(data)
  }
  for(i in colnames(data)){
    values=unique(data[[i]])
    values=setdiff(values,NA) # remove NA
    if(i %in% categoric){
      data[[i]]=factor(data[[i]],ordered=FALSE)
      print(paste("Metadata",i,"is categoric"))
    }else if(length(values)==2){
      data[[i]]=factor(data[[i]],ordered=FALSE)
      print(paste("Metadata",i,"is binary"))
    }else{
      print(paste("Metadata",i,"is numeric"))
      tryCatch(as.numeric(as.character(data[[i]])), warning=function(w) print(paste(i,"has non-numeric values!")))
      data[[i]]=as.numeric(as.character(data[[i]]))
    }
  } # loop column names
  return(data)
}

#' @title Remove groups that have metadata items with missing values
#' @param data a metadata matrix with rows as samples and columns as items or a dataframe with metadata items as columns
#' @param groups group membership vector with as many entries as samples
#' @return a vector with indices of samples without problematic groups
#' @export
removeGroupsWithMissingValues<-function(data, groups=c()){
  if(is.data.frame(data)==FALSE){
    data=as.data.frame(data)
  }
  unique.groups=unique(groups)
  indices.samples.tokeep=c()
  for(group in unique.groups){
    group.member.indices=which(groups==group)
    group.metadata=data[group.member.indices,]
    na.indices=which(is.na(group.metadata)==TRUE)
    # no missing metadata: keep the group
    if(length(na.indices)==0){
      indices.samples.tokeep=c(indices.samples.tokeep,group.member.indices)
    }else{
      print(paste("Removing group",group))
    }
  } # end group loop
  return(indices.samples.tokeep)
}

#' @title Replace missing values by group mean.
#'
#' @description Missing values are replaced either by the group mean
#' for numeric metadata or the most frequent group value for categoric
#' metadata. Missing values are only replaced if there are less
#' missing values than the threshold in the group, else the metadata item
#' concerned is removed. The output are metadata without missing values.
#'
#' @param metadata.df a data frame with metadata
#' @param groups a vector that specifies the group for each sample
#' @param na.threshold number of allowed missing values per group
#' @param metadata.to.skip metadata for which filling with the group mean does not make sense (e.g. date)
#' @export
setNAToGroupMean<-function(metadata.df, groups=c(), na.threshold=4, metadata.to.skip=c()){
  filter.indices=c()
  patientmean=NA
  # loop metadata items (which are columns)
  for(i in 1:length(names(metadata.df))){
    numberNA=length(which(is.na(metadata.df[[i]])))
    metadata.item=names(metadata.df)[i]
    if(numberNA==0){
      # nothing to do
    }else{
      if(numberNA<=na.threshold && length(which(metadata.to.skip==metadata.item))<1){
        # loop over samples (which are rows)
        for(j in 1:nrow(metadata.df)){
          if(is.na(metadata.df[j,i])){
            patient=groups[j]
            patient.indices=which(groups==patient)
            if(length(patient.indices)==1){
              # cannot replace NA, because only a single value is available
              print(paste("Metadata item",names(metadata.df)[i],"cannot be filled, because group",patient,"has only a single value!"))
              filter.indices=append(filter.indices,i)
            }else{
              # separate treatment for numeric and factor
              if(is.numeric(metadata.df[,i])){
                patientmean=mean(na.omit(metadata.df[patient.indices,i]))
              }else if(is.factor(metadata.df[,i])){
                values=unique(metadata.df[patient.indices,i])
                values=setdiff(values,c(NA))
                if(length(values)==1){
                  patientmean=values[1]
                  #print(paste("filled with",patientmean))
                }else{
                  print(paste("Metadata item",names(metadata.df)[i],"has more than one factor for group",patient,"!"))
                  patientmean=NA
                  filter.indices=append(filter.indices,i)
                }
              }
              metadata.df[j,i]=patientmean
              print(paste("Metadata item",names(metadata.df)[i],"is filled with group mean",patientmean,"for sample",rownames(metadata.df)[j]))
            }
          }
        } # loop over samples
      }# not too many missing values
      else{
        filter.indices=append(filter.indices,i)
        print(paste("Metadata item",names(metadata.df)[i],"has too many NA or is in the list of metadata to be skipped!"))
      }
    } # missing values present
  } # loop over metadata items
  keep.indices=setdiff(c(1:length(names(metadata.df))),filter.indices)
  return(metadata.df[,keep.indices])
}

# type can be either numeric, categoric, binary or catbin
# categoric does not include binary, but catbin includes both
getMetadataSubset<-function(metadata.df, type="numeric"){
  keep=c()
  for(i in 1:length(names(metadata.df))){
    if(type=="numeric" && is.numeric(metadata.df[,i])){
      keep=append(keep,i)
    }else if(is.factor(metadata.df[,i]) && (type=="categoric" || type=="binary" || type=="catbin")){
      if(type=="catbin"){
        keep=append(keep,i)
      }else{
        values=setdiff(unique(metadata.df[,i]),NA)
        if(type=="binary" && length(values)==2){
          keep=append(keep,i)
        }else if(type=="categoric" && length(values)>2){
          keep=append(keep,i)
        }
      }
    }
  }
  return(metadata.df[,keep])
}

# remove samples with indicated metadata item value
# if random is true, remove as many random samples as there are samples that have the particular metadata item value
# function returns filtered taxa and metadata objects in a list
filterSamplesWithMetadataValues<-function(taxa, metadata.df, metadata.item="", metadata.value=1, random=FALSE){
  metadata.values=metadata.df[[metadata.item]]
  to.filter=which(metadata.values==metadata.value)
  to.keep=setdiff(c(1:ncol(taxa)),to.filter)
  if(random){
    rand.indices=sample(c(1:ncol(taxa)))[1:length(to.keep)]
    to.keep=rand.indices
  }
  taxa=taxa[,to.keep]
  metadata.df=metadata.df[to.keep,]
  res=list(taxa, metadata.df)
  names(res)=c("taxa","metadata")
  return(res)
}
#' @title Get silhouette score for a clustering result
#'
#' @description The silhouette score is a cluster qualiy index that assessses both cluster cohesion
#' (how close points within the same cluster are to each other) and cluster separation
#' (how well different clusters are separated) for a given clustering result (encoded in the group membership vector).
#' The silhouette score is bounded by -1 (worst value) and 1 (best value). The silhouette score s of a point i is defined as:
#' s(i)=b(i)-a(i)/max(b(i),max(i)), where a(i) is the mean distance of this point to all other points within the same cluster
#' and b(i) is the smallest mean distance to data points in the nearest neighboring cluster.
#' The silhouette score of the entire clustering is computed as the mean of silhouette scores across all points.
#' This function wraps vegdist to compute distances/dissimilarities.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param groups vector with group memberships
#' @param method name of dissimilarity or distance to use. See vegdist for options
#' @export
#'
silhouette <- function(abundances, groups, method='bray'){
  if (length(groups) == length(colnames(abundances))){
    abundances <- t(abundances)
  } else if (length(groups) != length(rownames(abundances))){
    stop('Group vector is not the same length as the abundance table dimensions.')
  }
  #print(dim(abundances))
  dis <- as.matrix(vegdist(abundances, method=method))
  # for every sample in a cluster, calculate silhoutte coefficient
  scores <- c()
  empty=FALSE
  for (i in 1:nrow(dis)){
    row <- dis[i,]
    cluster <- groups[i]
    intra_cluster_distance <- 0
    smallest_cluster_distance <- 1
    for (clus in unique(groups)){
      ids <- which(groups == clus)
      if (clus == cluster){
        ids <- ids[ids != i]  # distance to self should not be included
        # Manual merge of treatment of zero-member clusters with alternative
        # <<<<<<< HEAD
        #        if(length(ids)==0){
        #          warning(paste("Cluster",cluster," of",length(unique(groups)),"clusters is empty! Silhouette is set to NA"))
        #          empty=TRUE
        #        }else{
        #          #print(paste("ids=",ids))
        #          intra_cluster_distance <- mean(row[ids])
        #        }
        #        # fix: replaced meandist below by intra_cluster_distance
        #      } else if (intra_cluster_distance < smallest_cluster_distance){
        #        smallest_cluster_distance <- mean(row[ids])
        # =======
        # only done if there are members in the cluster: intra_cluster_distance <- mean(row[ids])
        if(length(ids)==0){
          warning(paste("Cluster",cluster," of",length(unique(groups)),"clusters is empty! Silhouette is set to NA"))
          empty=TRUE
        }else{
          #print(paste("ids=",ids))
          intra_cluster_distance <- mean(row[ids])
        }
      } else {
        inter_cluster_distance <- mean(row[ids])
        if (inter_cluster_distance < smallest_cluster_distance){
          smallest_cluster_distance <- inter_cluster_distance
        }
        # >>>>>>> 0d6c5931d0584caae022c3a1319236ef0c91ed4c
      }
    } # loop clusters
    if (intra_cluster_distance < smallest_cluster_distance){
      silhouette <- 1-(intra_cluster_distance / smallest_cluster_distance)
      
    } else if (intra_cluster_distance > smallest_cluster_distance){
      silhouette <- (smallest_cluster_distance / intra_cluster_distance) - 1
    } else {
      silhouette <- 0
    }
    if(empty){
      silhouette=NA
    }
    #print(paste("silhouette",silhouette))
    scores <- c(scores, silhouette)
  } # loop rows of dissimilarity matrix
  #print(scores)
  return(mean(scores))
}

#Filtering out taxa
prevFilter = function(phy, prev){
  #' @title Prevalence filter
  #'
  #' @description Filters taxa present in fewer samples than the specified fraction. 
  #' @details Filtered taxa are retained in the Bin taxon to preserve sample sums. 
  #'
  #' @param phy Phyloseq object
  #' @param prev Minimum fraction of samples
  #'
  prev0 = apply(X = otu_table(phy),
                MARGIN = ifelse(taxa_are_rows(phy), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
  prevalenceThreshold = prev * nsamples(phy)
  nonrares = prune_taxa((prev0 > prevalenceThreshold), phy)
  rares = prune_taxa((prev0 < prevalenceThreshold), phy)
  rares = merge_taxa(rares, taxa_names(rares))
  otus = data.frame(otu_table(nonrares))
  otus = cbind(otus, data.frame(otu_table(rares)))
  tax = data.frame(tax_table(nonrares), stringsAsFactors = FALSE)
  tax = rbind(tax, rep("Other", 7))
  otus <- t(otus)
  rownames(otus) <- c(taxa_names(nonrares), 'Bin')
  rownames(tax) <- c(taxa_names(nonrares), 'Bin')
  newphy = phyloseq(otu_table(otus, taxa_are_rows = TRUE), sample_data(phy), tax_table(as.matrix(tax)))
  return(newphy)
}



compareGroups<-function(abundances, property="beta", method="dissim", groups=c(), colors=c(), colorMap=NULL, plot.type="pergroup", avg="none", all=FALSE, noSameGroup=TRUE, rowNorm=FALSE, subsample=TRUE, xlab="", pvalViz=FALSE, pAdjMethod="BH", returnValues=FALSE){
  supported.properties=c("richness","evenness","alpha","beta")
  if(plot.type=="intravsinter" && property!="beta"){
    stop("The intravsinter plot type is only supported for beta diversity!")
  }
  if(plot.type=="intravsinter" && method=="DM"){
    stop("The intravsinter plot type is not supported for method DM.")
  }
  if(length(groups)==0){
    if(plot.type=="intravsinter"){
      stop("For plot type intravsinter, a vector of group memberships is required.")
    }else{
      warning("Empty group vector provided. All samples are assigned to the same group.")
      groups=rep("all",ncol(abundances))
    }
  }
  if(plot.type=="intravsinter"){
    all=TRUE
    property="beta"
  }
  if(method %in% supported.properties){
    stop(paste("The method parameter specifies the way in which beta diversity is computed. Please use the property parameter."))
  }
  scaleMatrix=FALSE
  if(property=="richness"){
    if(!is.Count.Matrix(abundances)){
      warning("Richness estimation requires counts. Data were scaled to counts.")
      scaleMatrix=TRUE
    }
  }
  if(property=="beta" && method=="DM"){
    if(!is.Count.Matrix(abundances)){
      warning("The Dirichlet-Multinomial distribution can only be estimated from count data. Data were scaled to counts.")
      scaleMatrix=TRUE
    }
  }
  if(scaleMatrix){
    vec=as.vector(abundances)
    scaling.factor=1/min(vec[vec!=0])
    if(scaling.factor>10000){
      scaling.factor=10000
    }
    print(paste("Abundances have been scaled with factor:",scaling.factor))
    abundances=round(abundances*scaling.factor)
  }
  if(property!="beta"){
    method=""
  }
  
  unique.groups=unique(groups)
  if(is.numeric(unique.groups)){
    unique.groups=sort(unique.groups)
  }
  
  las=1
  if(length(unique(groups))>10){
    las=2 # vertical x axis labels
  }
  
  if(length(colors)>0 && length(groups)>0){
    if(length(colors)!=length(unique(groups))){
      #print(length(unique(groups)))
      #print(length(colors))
      stop("Please provide as many colors as there are groups!")
    }
  }
  
  if(length(colors)==0 && length(groups)>0){
    colors=rep("white",length(groups))
  }else if(length(colors)==0 && length(groups)==0){
    # default color if no groups are provided
    colors=c("black")
  }
  
  pseudo=min(abundances[abundances>0])/100 # to take logarithm for Taylor law
  #print(paste("pseudo:",pseudo))
  # computing intra-group variability
  groups.with.theta=c()
  interdissim=c()
  intradissim=c()
  pergroupdissim=list()
  groupspecdissim=c()
  groupspecavgprops=c()
  main=""
  ylab=""
  
  # return
  res=NULL
  
  if(xlab==""){
    if(plot.type=="intravsinter"){
      xlab="Bray-Curtis dissimilarity"
    }else{
      xlab="Groups"
    }
  }
  
  if(rowNorm == TRUE){
    # normalize row-wise
    abundances=abundances/rowSums(abundances)
  }
  
  constrainSampleNum=FALSE
  sampleNum=NA
  
  if(subsample==TRUE){
    constrainSampleNum=TRUE
    sampleNum=min(table(groups))
    print(paste("Constraining sample number randomly to the same minimal group sample number of",sampleNum))
  }
  
  # loop groups
  for(group in unique.groups){
    groupspecdissim=c()
    #print(paste("Processing group",group))
    group.indices=which(groups==group)
    if(constrainSampleNum){
      group.indices=sample(group.indices)[1:sampleNum]
    }
    group.data=abundances[,group.indices]
    if(length(group.indices)==1){
      group.data=as.matrix(group.data)
    }
    if(property=="beta"){
      if(is.null(group.data)==FALSE && ncol(group.data) > 0){
        groups.with.theta=c(groups.with.theta,group)
        if(method == "dissim"){
          main="Beta diversity (Bray Curtis dissimilarities)"
          ylab="Bray Curtis dissimilarity"
          dissimMat=as.matrix(vegdist(t(group.data),method="bray"))
          dist=as.matrix(vegdist(t(group.data),method="bray"))
          pergroupdissim[[as.character(group)]]=dist[lower.tri(dist)]
          # average beta diversity
          if(avg!="none"){
            if(avg=="median"){
              groupspecavgprops=c(groupspecavgprops,median(groupspecdissim))
            }else if(avg=="mean"){
              groupspecavgprops=c(groupspecavgprops,mean(groupspecdissim))
            }
          }
        }else if(method == "DM"){
          main="Beta diversity (overdispersion)"
          ylab="Overdispersion"
          # dirmult: counts of alleles X (items=columns) vary across sub-populations (observations=rows), so transpose is necessary
          dm.fit=dirmult(t(group.data))
          intradissim=c(intradissim, dm.fit$theta)
        } # end beta diversity methods
      } # avoid patients with a single column
      else{
        warning(paste("Group",group,"has less than 2 samples: cannot compute beta diversity."))
      }
    }else if(property=="taylor"){
      main="Taylor law slope"
      ylab="Slope"
      means=apply(group.data,1, mean, na.rm=TRUE) # compute taxon-wise
      vars=apply(group.data,1, var, na.rm=TRUE)
      logvars=log(vars+pseudo)
      logmeans=log(means+pseudo)
      reg.data=data.frame(logvars,logmeans)
      linreg = lm(formula = logvars~logmeans)
      # print(paste("Intercept:",linreg$coefficients[1]))
      # print(paste("Slope:",linreg$coefficients[2]))
      slope=linreg$coefficients[2]
      intradissim=c(intradissim, slope)
    }else if(property=="alpha" || property=="richness" || property=="evenness"){
      for(sample.index in 1:ncol(group.data)){
        if(property=="alpha"){
          main="Alpha diversity (Shannon index)"
          ylab="Shannon"
          # Shannon
          value=vegan::diversity(group.data[,sample.index])
        }else if(property=="richness"){
          main="Richness (Chao1 index)"
          ylab="Chao1"
          # Chao1
          value=estimateR(group.data[,sample.index])[2]
        }else if(property=="evenness"){
          main="Evenness (Sheldon index)"
          ylab="Sheldon"
          value=sheldon(group.data[,sample.index])
        }
        groupspecdissim=c(groupspecdissim,value)
        # average values over group members or keep all values
        if(avg=="mean"){
          groupspecavgprops=c(groupspecavgprops,mean(groupspecdissim))
        }else if(avg=="median"){
          groupspecavgprops=c(groupspecavgprops,median(groupspecdissim))
        }else if(avg=="none"){
          pergroupdissim[[as.character(group)]]=groupspecdissim
        }
      }
    }else{
      stop(paste("Property",property,"is not supported."))
    }
  } # end group loop
  
  # compute global property
  if(all){
    if(property=="beta"){
      if(method == "dissim"){
        dissimMat=as.matrix(vegdist(t(abundances),method="bray"))
        for(index1 in 1:(ncol(abundances)-1)){
          for(index2 in (index1+1):ncol(abundances)){
            group1=groups[index1]
            group2=groups[2]
            dissimVal = dissimMat[index1,index2]
            if(!noSameGroup || group1 != group2){
              interdissim=c(interdissim,dissimVal)
            }
          }
        }
      }else if(method == "DM"){
        dm.fit=dirmult(t(abundances))
        interdissim=dm.fit$theta
      }
    }else{
      warning("All is only supported for beta diversity.")
      all=FALSE
    }
  } # end all
  
  if(las==2){
    xlab=""
  }
  #if(!is.null(colorMap)){
  # TODO to make space in par for the legend and shift it there, no longer needed for violin plots, but still required for barplots
  #}
  
  # plot
  if(plot.type=="pergroup"){
    if(method=="DM"){
      thetas=intradissim
      print(interdissim)
      names.arg=groups.with.theta
      if(all==TRUE){
        names.arg=c(groups.with.theta,"all")
        thetas=c(intradissim,interdissim)
      }
      range=c(0,max(thetas))
      res=thetas
      compareGroupsPlot(mat=thetas, labels=names.arg, xlab=xlab, ylab="Theta", main="Estimated overdispersion", range=range, plotType = "bar",colors=colors, las=las, colorMap=colorMap)
    }else if(property=="taylor"){
      range=c(min(intradissim),max(intradissim)) # slope can be negative
      if(min(intradissim)>0){
        range=c(0,max(intradissim))
      }
      compareGroupsPlot(mat=intradissim, labels=unique.groups, xlab=xlab, ylab=ylab, main=main, range=range, plotType = "bar",colors=colors, las=las, colorMap=colorMap)
    }else{
      if(avg!="none"){
        names(groupspecavgprops)=unique.groups
        if(all && property=="beta"){
          avgall=NA
          if(avg=="mean"){
            avgall=mean(interdissim)
          }else if(avg=="median"){
            avgall=median(interdissim)
          }
          groupspecavgprops=c(groupspecavgprops,avgall)
          names(groupspecavgprops)=c(unique.groups,"all")
        }
        range=c(0,max(groupspecavgprops))
        if(property=="beta"){
          range=c(0,1)
        }
        ylab=paste(avg,ylab)
        res=groupspecavgprops
        compareGroupsPlot(mat=groupspecavgprops, xlab=xlab, ylab=ylab, main=main, range=range, plotType = "bar",colors=colors, las=las, colorMap=colorMap)
      }else{
        # richness, evenness and alpha-diversity that are not averaged as well as beta-diversity with method dissim are displayed as box plots
        if(all && property=="beta"){
          pergroupdissim[["all"]]=interdissim
        }
        mat=listToMat(pergroupdissim)
        ylim=c(0,max(mat,na.rm = TRUE))
        res=mat
        # display box plot with p-values
        if(pvalViz){
          combinations=list()
          units=colnames(mat)
          # go through upper triangle
          for(index1 in 1:(length(units)-1)){
            for(index2 in (index1+1):length(units)){
              unit1=units[index1]
              unit2=units[index2]
              val1=unique(mat[,index1])
              val2=unique(mat[,index2])
              combi=paste(unit1,unit2,sep="")
              # avoid result only consisting of missing values (e.g. if group has only 1 sample)
              if((length(val1)==1 && is.na(val1)) || (length(val2)==1 && is.na(val2))){
                print("Skipping group combination",combi,"because one group contains missing values or has only 1 sample")
              }else{
                combinations[[combi]]=c(unit1,unit2)
              }
            }
          }
          compareGroupsPlot(mat=mat, groups=groups, colors=colors, combinations=combinations, pAdjMethod = pAdjMethod, xlab=xlab, ylab=ylab, main=main, plotType = "pvalues")
        }else{
          # alternative plot type: box, violin
          compareGroupsPlot(mat=mat, xlab=xlab, ylab=ylab, groups=groups, main=main, colors=colors, las=las, colorMap=colorMap, plotType = "violin")
        }
      }
    }
  }else if(plot.type=="intravsinter"){
    intradissim=as.vector(listToMat(pergroupdissim))
    print(paste("Number of inter-group values:",length(interdissim)))
    print(paste("Number of intra-group values:",length(intradissim)))
    if(subsample){
      print("Number of inter and intra-group values is equalised.")
      if(length(intradissim)<length(interdissim)){
        mat=cbind(intradissim,interdissim[sample(1:length(intradissim))])
      }else{
        mat=cbind(intradissim[sample(1:length(interdissim))],interdissim)
      }
    }else{
      dissimList=list(intradissim, interdissim)
      names(dissimList)=c("intra","inter") # names are required here
      mat=listToMat(dissimList)
    }
    colnames(mat)=c("intra","inter")
    out=wilcox.test(mat[,1],mat[,2],na.rm=TRUE) # unpaired, two-sided
    pval=round(out$p.value,4)
    min=min(mat,na.rm =TRUE)
    max=max(mat,na.rm =TRUE)
    maxD=max(max(hist(mat[,1],plot=F)$counts, na.rm=TRUE),max(hist(mat[,2],plot=F)$counts, na.rm=TRUE))
    #print(min)
    #print(max)
    #print(maxD)
    cols=c(rgb(0,1,0,0.5),rgb(1,0,0,0.5))
    res=mat
    hist(mat[,1],xlim=c(min,max),ylim=c(0,maxD), main=paste("P-value Wilcoxon",pval,sep=": "), xlab=xlab, col=cols[1], breaks="FD")
    hist(mat[,2],col=cols[2], breaks="FD", add=T)
    legend("topleft", colnames(mat), lty = rep(1,length(cols)), col = cols, merge = TRUE, bg = "white", text.col="black")
  }else{
    stop(paste("Plot type",plot.type,"is not supported."))
  }
  
  if(returnValues){
    return(res)
  }
}

# Internal function for group-wise property plotting.
# mat: values to be plotted
# xlab: label for x axis
# ylab: label for y axis
# range: y axis range
# main: title of the plot
# groups: group membership vector
# labels: bar plot labels
# colors: bar/box color
# combinations: for pvalue plots
# pAdjMethod: multiple testing correction method for pvalue plots
# las: orientation of axis labels
# plotType: box, bar, pvalues or violin
# colorMap: color map for a metadata item, adds color legend
compareGroupsPlot<-function(mat, xlab="", ylab="", range=c(), main="", groups=c(), labels=c(), colors=c(), combinations=c(), pAdjMethod="BH", las=1, plotType="box", colorMap=NULL){
  
  legend.colors=c()
  if(!is.null(colorMap)){
    for(entry in names(colorMap)){
      legend.colors=c(legend.colors,colorMap[[entry]])
    }
  }
  
  if(plotType=="box"){
    # varwidth does not work, notch generates a warning message
    boxplot(mat,ylab=ylab, main=main, xlab=xlab, notch=FALSE,ylim=range, col = colors, cex.names = 0.9, las=las) # border=colors
    for(i in 1:ncol(mat)){
      points(rep(i,length(mat[,i])),mat[,i])
    }
    # make a plot with ggplot2
  }else if(plotType=="violin" || plotType=="pvalues"){
    # ggplot2 and reshape2 are imported by phyloseq, which is imported by seqgroup
    # stat_compare_means requires package ggpubr, also imported by seqgroup
    # to avoid error messages in package built
    variable=""
    value=""
    color=""
    df=as.data.frame(mat)
    df_melt <- melt(df)
    # if a color map is given, add color as third factor
    if(!is.null(colorMap)){
      # melt color map into two vectors
      legend.names=names(colorMap)
      legend.colors=c()
      for(name in names(colorMap)){
        legend.colors=c(legend.colors,colorMap[[name]])
      }
      #print(colors)
      #print(legend.names)
      #print(legend.colors)
      df_melt_extended=matrix(nrow=nrow(df_melt),ncol=3)
      unique.groups=unique(groups)
      #print(unique.groups)
      # loop rows of df_melt object and add color
      for(i in 1:nrow(df_melt)){
        group=as.character(df_melt[i,1])
        group.index=which(unique.groups==group)
        df_melt_extended[i,1]=group
        df_melt_extended[i,2]=as.numeric(as.character(df_melt[i,2]))
        # find the name of the metadata item in the color map
        color.name.index=which(legend.colors==as.character(colors[group.index]))
        df_melt_extended[i,3]=as.character(legend.names[color.name.index])
      }
      colnames(df_melt_extended)=c("variable","value","color")
      df_melt_extended=as.data.frame(df_melt_extended)
      # make sure value is not treated as a factor
      df_melt_extended$value=as.numeric(as.character(df_melt_extended$value))
    }
    if(plotType=="pvalues"){
      # cannot set ylim, else p-values are not plotted correctly
      if(!is.null(colorMap)){
        p <- ggplot(df_melt_extended, aes(variable, value, fill=color))+geom_violin()+ ggpubr::stat_compare_means(comparisons=combinations, p.adjust.method = pAdjMethod)+xlab(xlab)+ylab(ylab)+ggtitle(main) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # 90 degree x axis test
        plot(p)
      }else{
        # geom_jitter(position = position_jitter(0.2))
        p <- ggplot(df_melt, aes(variable, value, fill=variable))+ guides(fill=FALSE) +geom_violin(show.legend=FALSE)+ scale_fill_manual(values=colors) + ggpubr::stat_compare_means(comparisons=combinations, p.adjust.method = pAdjMethod)+xlab(xlab)+ylab(ylab)+ggtitle(main) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # 90 degree x axis test
        plot(p)
      }
    }else{
      if(!is.null(colorMap)){
        p <- ggplot(df_melt_extended, aes(variable,value, fill=color)) + geom_violin()  +  xlab(xlab) + ylab(ylab) + ggtitle(main) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        plot(p)
      }else{
        # to remove legend: http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
        p <- ggplot(df_melt, aes(variable,value, fill=variable)) + guides(fill=FALSE) + geom_violin()  + scale_fill_manual(values=colors) + xlab(xlab) + ylab(ylab) + ggtitle(main) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        plot(p)
      }
    }
  }else if(plotType=="bar"){
    if(length(labels)>0){
      barplot(mat,names.arg=labels, cex.names = 0.9, xlab=xlab,ylab=ylab,main=main,ylim=range, col = colors, las=las)
    }else{
      barplot(mat,xlab=xlab,ylab=ylab,main=main, ylim=range, col=colors,cex.names = 0.9, las=las)
    }
  }else{
    stop(paste("Plot type",plotType,"is not supported."))
  }
  if(!is.null(colorMap) && (plotType!="pvalues") && (plotType!="violin")){
    legend("topright",legend=names(colorMap), cex=0.9, bg = "white", text.col=legend.colors)
  }
}


#' @title Assign colors such that they encode a metadata item
#' @description For a binary or categoric metadata item, as many colors as categories exist are assigned.
#' For a numeric metadata item, the item is binned according to the equalfreq or equalwidth
#' strategy in as many bins as requested. Alternatively, the user can also provide thresholds to define bins.
#' Each bin is assigned its color. If the groups vector is empty, a color vector is returned with as many entries as samples.
#' If groups are given, a group color vector is returned with as many entries as groups. The color map is returned along with
#' the color vector.
#' @details Notation: a round bracket means that the end point is excluded, a square bracket that the end point is included.
#' For groupColorStrategy maxfreq: if there are several equally frequent bins in the group, the first one is selected.
#' Missing values receive color "gray".
#' @param metadataItem a vector with binary, categoric or numeric values
#' @param groups an optional group vector with as many entries as samples in the target data set
#' @param binning binning strategy, only required for numeric metadata, either a string (equalfreq: the same number of samples in each group or equaldist: the same range of values) or a vector with thresholds
#' @param numBins the number of bins into which to bin a numeric metadata item, if zero: as many bins as samples (only required for numeric metadata and equalfreq/equaldist binning)
#' @param binLabels optional names for bins
#' @param groupColorStrategy maxfreq (assign the most frequent bin color to the group), binnum (color encodes number of bins per group)
#' @param returnBinned return the binned metadataItem
#' @return a list with color vector and a color map, if returnBinned is true, in addition bins
#' @export
makeColorsGivenMetadata<-function(metadataItem=c(), groups=c(), binning="equalfreq", numBins=0, binLabels=c(), groupColorStrategy="maxfreq",returnBinned=FALSE){
  
  bins=numBins
  if(bins==0){
    bins=length(metadataItem)
  }
  
  # vector with sample membership in binned metadata item
  binned=c()
  # vector with sample-wise or group-wise colors
  colors=c()
  # map to look up color assigned to bin
  colorMap=list()
  # result list
  res=list()
  # indices of missing values
  na.indices=which(is.na(metadataItem))
  
  # numeric metadata item and binning requested
  if(is.numeric(metadataItem) && bins<length(metadataItem)){
    print("Metadata item is numeric.")
    # binning strategy: thresholds were provided
    if(is.numeric(binning)){
      # check if minimum and maximum are provided as thresholds
      if(!(min(metadataItem) %in% binning)){
        binning=c(min(metadataItem),binning)
        print(paste("Adding minimum (",min(metadataItem),") as bin."))
      }
      if(!(max(metadataItem) %in% binning)){
        binning=c(binning,max(metadataItem))
        print(paste("Adding maximum (",max(metadataItem),") as bin."))
      }
      if(length(binLabels)>0){
        if(length(binLabels)==length(binning) || length(binLabels)>length(binning)){
          stop("There should be one label less than binning thresholds!")
        }
        binned=cut(metadataItem,breaks=binning, include.lowest=TRUE,labels=binLabels)
      }else{
        binned=cut(metadataItem,breaks=binning, include.lowest=TRUE)
      }
      # equaldist strategy
    }else if(binning == "equaldist"){
      # bin vector (equaldist)
      cut.out=cut(metadataItem,breaks=bins)
      print("Thresholds equaldist binning:")
      print(levels(cut.out))
      # rename bins (previously named after borders)
      if(length(binLabels)==bins){
        levels(cut.out)=binLabels
      }
      binned=cut.out
      # equalfreq strategy
    }else if(binning == "equalfreq"){
      # divide into requested number of intervals
      intervals=seq(0,1,(1/bins))
      thresholds=quantile(metadataItem,probs=intervals, na.rm=TRUE)
      print("Thresholds equalfreq binning:")
      print(thresholds[2:(length(thresholds)-1)])
      binned=metadataItem
      for(i in 2:length(thresholds)){
        set1=which(metadataItem>=thresholds[(i-1)])
        set2=which(metadataItem<thresholds[i])
        indices=intersect(set1,set2)
        if(length(binLabels)>0){
          binned[indices]=binLabels[(i-1)]
        }else{
          if(i==2){
            binned[indices]=paste("[",thresholds[(i-1)],",",thresholds[i],"]",sep="")
          }else if(i==length(thresholds)){
            binned[indices]=paste("[",thresholds[(i-1)],",",thresholds[i],"]",sep="")
          }else{
            binned[indices]=paste("(",thresholds[(i-1)],",",thresholds[i],")",sep="")
          }
        }
      }
      maxindex=which(metadataItem==max(metadataItem, na.rm=TRUE))
      if(length(binLabels)>0){
        binned[maxindex]=binLabels[length(binLabels)]
        levels(binned)=binLabels
      }else{
        binned[maxindex]=paste("[",thresholds[(length(thresholds)-1)],",",thresholds[length(thresholds)],"]",sep="")
      }
    }else{
      stop("Binning strategy ",binning," not supported!")
    }
    # do not bin numeric metadata item
  }else if(bins == length(metadataItem)){
    print("Each sample is assigned into its own bin.")
    binned = metadataItem
  }else{
    print("Metadata item is binary/categoric")
    bins=length(unique(metadataItem))
    print(paste("Number of categories:",bins))
    binned = metadataItem
  }
  
  # drop levels efficiently
  binned=as.character(binned)
  
  # treat missing values
  if(length(na.indices)>0){
    print(paste(length(na.indices),"missing values found."))
    binned[na.indices]="NA" # will be assigned a gray color
    colorMap[["NA"]]="gray"
  }
  
  #print(binned)
  
  # get a vector with as many colors as samples
  colors=assignColorsToGroups(groups=binned, refName="NA")
  
  # assemble color map
  for(index in 1:length(binned)){
    abin=binned[index]
    if(!(abin %in% names(colorMap))){
      colorMap[[as.character(abin)]]=colors[index]
    }
  }
  
  #print(names(colorMap))
  
  if(length(groups)>0){
    group.colors=c()
    unique.groups=unique(groups)
    hues=c()
    # count number of bins in each group to identify the range
    bins.per.group=c()
    group.counter=1
    if(groupColorStrategy=="binnum"){
      colorMap=list() # reset color map
      for(group in unique.groups){
        group.member.indices=which(groups==group)
        bins.per.group=c(bins.per.group,length(unique(binned[group.member.indices])))
      }
      col.vec = seq(0,1,1/max(bins.per.group))
      hues = hsv(col.vec)
    }
    # loop groups
    for(group in unique.groups){
      group.member.indices=which(groups==group)
      if(groupColorStrategy=="maxfreq"){
        # get the name of the most frequent bin in the group
        maxfreq.bin.label=names(sort(table(binned[group.member.indices]),decreasing = TRUE))[1]
        #print(paste("Most frequent bin in group",group,":",maxfreq.bin.label))
        # get an index of the bin
        maxfreq.bin.index=which(binned==maxfreq.bin.label)[1]
        # get the corresponding color
        maxfreq.bin.color=colors[maxfreq.bin.index]
        group.colors=c(group.colors,maxfreq.bin.color)
      }else if(groupColorStrategy=="binnum"){
        bin.number=bins.per.group[group.counter]
        group.colors=c(group.colors,hues[bin.number])
        if(!(bin.number %in% names(colorMap))){
          colorMap[[bin.number]]=hues[bin.number]
        }
      }else{
        stop(paste("Group color strategy",groupColorStrategy,"not supported."))
      }
      group.counter=group.counter+1
    } # end group loop
    colors=group.colors
  }
  
  if(returnBinned){
    res=list(colors, colorMap, binned)
    names(res)=c("colors","colormap","bins")
  }else{
    res=list(colors,colorMap)
    names(res)=c("colors","colormap")
  }
  return(res)
}

# Convert a list into a matrix. When a list entry
# has less values than the entry with the largest
# number of values, complete it with missing values.
listToMat<-function(groupprops=list()){
  # get group with the largest number of entries
  maxnum=0
  #print(names(groupprops))
  for(name in names(groupprops)){
    if(length(groupprops[[name]])>maxnum){
      maxnum=length(groupprops[[name]])
    }
  }
  mat=matrix(NA,nrow=maxnum,ncol=length(names(groupprops)))
  colnames(mat)=names(groupprops)
  counter=1
  for(name in names(groupprops)){
    for(index in 1:length(groupprops[[name]])){
      mat[index,counter]=groupprops[[name]][index]
    }
    counter=counter+1
  }
  return(mat)
}

# Compute evenness using Sheldon's index
#
# Sheldon's index is defined as \eqn{S=\frac{e^H}{N}}, where H is the Shannon diversity and N the species number.
# It ranges from 0 to 1, where 1 signifies a perfectly even abundance distribution.
#
#  A.L. Sheldon 1969. Equitability indices: dependence on the species count. Ecology, 50, 466-467.
# C Heip 1974. A new index measuring evenness. J. mar. biol. Ass. UK 54, 555-557.
#
#  Note that the N2N1 mode results in evenness smaller than 1 for equal taxon probabilities.
#
# x a vector of species abundances
# correction whether or not to apply the correction described in Alatalo, Oikos 37, 199-204, 1981
# N2N1 whether to compute Sheldon's evenness as the ratio of e raised to the power of H (H = Shannon diversity) and Simpson's diversity
# Sheldon's evenness
sheldon<-function(x, correction = TRUE, N2N1 = FALSE){
  H = diversity(x, index="shannon")
  if(N2N1){
    simpson = diversity(x, index="simpson")
    numerator = 1/simpson
    denominator = exp(1)^H
  }else{
    numerator = exp(1)^H
    denominator = specnumber(x)
  }
  if(correction){
    numerator = numerator - 1
    denominator = denominator - 1
  }
  S = numerator/denominator
  S
}


















