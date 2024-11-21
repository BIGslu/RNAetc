#'  Calculate WGCNA module coherence
#'
#' Calculate within module pairwise gene correlations from WGCNA modules and gene expression data.
#' Test the coherence of modules built from Dataset A in data from Dataset B.
#'
#'
#' @param module_gene_sets Dataframe giving module membership. Should contain columns specifying ensembl gene IDs & module definitions.
#' @param voom_obj Voom object containing expression dataset in which module coherence is to be tested. Should either me a voom EList with "E" expression matrix, or an expression matrix.
#' @param gene_set_column Character string specifying name of the column in module_gene_sets df which defines modules. Defaults to "mod.color" which is output by make_modules()
#' @param ensembl_ID_column Character string specifying name of the column in module_gene_sets df which gives gene IDs. Defaults to "geneName" which is output by make_modules()
#' @param module_set Character string specifying name the study/dataset from which modules were built. This is used simply for labeling of outputs. Default = "STUDY"
#' @param sample_set Character string specifying name of the study from which the data come. This is used simply for labeling of outputs. Default = "STUDY"
#' @param remove_sets Vector of character strings naming modules which you want removed from the analysis, eg. "0" or "Grey".
#' @param return_plot Logical indicating whether plot should be printed when function runs. Defaults = TRUE.
#' @param r_cutoff Vector or single numeric value at which to draw red cutoff lines in correlation plot. Default = 0.3
#' @param p_cutoff Vector or single numeric value at which to draw red cutoff lines in significance plot. Default = 0.01
#'
#' @return List including:
#' \itemize{
#'   \item{coherence_boxplot_combined} A ggplot plot object visualizing distributions of correlation and p values.
#'   \item{coherence_boxplot_cor} A ggplot plot object visualizing distributions of correlation values.
#'   \item{coherence_boxplot_p} A ggplot plot object visualizing distributions of correlation p values.
#'   \item{subgene_correlation_df} Dataframe containing pairwise gene correlations within modules.
#'   \item{GSabsent} A list of missing genes enumerating genes defined in modules but missing from expression data.
#' }#'
#' @export
#'
#' @examples
#' moduleCoherence <- calculate_module_coherence(module_gene_sets = example.mods,
#' voom_obj = example.voom, module_set = "Example WGCNA Modules",
#' sample_set = "Example Data", r_cutoff = 0.3, p_cutoff = 0.01)
calculate_module_coherence<-
  function(module_gene_sets,
           voom_obj,
           gene_set_column = "mod.color",
           ensembl_ID_column="geneName",
           module_set="STUDY",
           sample_set="STUDY",
           remove_sets = NULL,
           return_plot=TRUE,
           r_cutoff = 0.3,
           p_cutoff = 0.01){

    # set undefined global variables to null
    X1 <- X2 <- allCorVals<- allPVals <- V3<- Cor<- P <- Set <- gene1<- gene2<- P.est<- h_label<- h_pos<- negLogP<- value<- h_line<- NULL

    moduleSets<- module_gene_sets
    if(class(voom_obj)[1]== "EList"){voomData<-voom_obj$E}else{voomData<-voom_obj}
    if(!is.numeric(voomData)){
      print("Warning, count/expression matrix is non-numeric.
        This usually arises when a non-expression column has been retained from a metadata object (eg due to grouping).
        Be advised to check for stray metadata in the expression matrix.")}


    #---------------------------------------------------------
    #  Calculate the gene set means (module expression)
    #---------------------------------------------------------
    uniGS<-moduleSets%>%
      dplyr::pull(gene_set_column)%>%
      unique()%>%as.character()%>%gtools::mixedsort() # Pull unique modules & order correctly

    GSsub<-c()
    GSabsent<-list()

    for (i in uniGS){
      curGS<-i
      curIDs<-as.character(moduleSets[which(moduleSets[, gene_set_column]==curGS), ensembl_ID_column])

      matchIndex<-match(curIDs, rownames(voomData))
      if (any(is.na(matchIndex))){
        matchIndex<-matchIndex[-which(is.na(matchIndex))]}
      GSabsent[[i]]<-curIDs[which(!(curIDs %in% rownames(voomData)))]
      if (length(curIDs[which(curIDs %in% rownames(voomData))]) > 1){
        GSsub<-rbind(GSsub, apply(voomData[matchIndex,], 2, mean))
      }else{
        GSsub<-rbind(GSsub, voomData[matchIndex,])}
    }

    #rownames(GSsub)<-uniGS

    if(length(unlist(GSabsent))>0){
      print(paste0(length(unique(unlist(GSabsent))), " Module Genes absent in voom object. Missing genes accesible in .$GSabsent"))}

    #----------------------------
    #  Set function for correlation p values
    #----------------------------

    # Correlation p value function
    cor2pvalue = function(r, n) {
      t <- (r*sqrt(n-2))/sqrt(1-r^2)
      p <- 2*stats::pt(abs(t),(n-2), lower.tail=FALSE)
      se <- sqrt((1-r*r)/(n-2))
      out <- list(r, n, t, p, se)
      names(out) <- c("r", "n", "t", "p", "se")
      return(out)
    }


    # remove gene sets that are zero
    if(!is.null(remove_sets)){uniGS<-uniGS[-which(uniGS %in% remove_sets)]}

    # Setup storage
    SubGeneCorDF<-list()
    SubGeneCorMedian<-c()
    SubGenePMedian<-c()

    for (i in uniGS){
      curGS<-i
      curIDs<-as.character(moduleSets[which(moduleSets[, gene_set_column]==curGS), ensembl_ID_column])
      matchIndex<-match(curIDs, rownames(voomData))
      if (any(is.na(matchIndex))){
        matchIndex<-matchIndex[-which(is.na(matchIndex))]}

      if(length(matchIndex)>1){
        setCor<-stats::cor(t(voomData[matchIndex,])) # Calculate pairwise gene correlations from expression data
        n<-t(!is.na(t(voomData[matchIndex,]))) %*% (!is.na(t(voomData[matchIndex,]))) # create matrix with sample size to calculate correlation p (dim=curIDs*curIDs, value = sample size)

        allCorInfo<-cor2pvalue(setCor, n)
        setP<-allCorInfo$p

        SubGeneCorDF[[i]]<- # Assemble df of gene pairs, correlations, and p values
          data.frame(t(utils::combn(colnames(setCor), 2)),
                     allCorVals=setCor[t(utils::combn(colnames(setCor), 2))],
                     allPVals = setP[t(utils::combn(colnames(setP), 2))])%>%
          dplyr::rename(gene1 = X1, gene2 = X2)%>%
          dplyr::mutate(V3=i)} else {print(paste0("Faulty gene-set:", curGS, " only contains 1 gene present in expression data. Skipping."))}
    }


    # Assemble and format results dataframe
    SubGeneCorDF<-
      dplyr::bind_rows(SubGeneCorDF)%>%
      dplyr::rename(Cor = allCorVals, P = allPVals, Set = V3)%>%
      dplyr::mutate(Cor = as.numeric(as.character(Cor)))%>%
      dplyr::mutate(P = as.numeric(as.character(P)))%>%
      dplyr::mutate(Set = factor(Set, levels =uniGS))%>%
      dplyr::select(Set, gene1, gene2, Cor, P)

    #Warning if 0 exist in p-values
    if(min(SubGeneCorDF$P) ==0){
      print(paste(length(SubGeneCorDF$P[SubGeneCorDF$P==0]),
                  "p-values = 0. In plots, these will be replaced with the lowest non-zero value",
                  formatC(sort(unique(SubGeneCorDF$P))[2], digits=2, format="e")))
    }


    #Boxplot of all pairwise correlations per module
    coherence_boxplot_cor<-
      SubGeneCorDF%>%
      ggplot2::ggplot(ggplot2::aes(y=Cor, x=Set))+
      ggplot2::geom_hline(yintercept =0, color="black")+
      ggplot2::geom_boxplot(outlier.shape = NA)+
      ggplot2::geom_hline(yintercept = r_cutoff, linetype = "dashed", color = "red")+
      ggplot2::scale_y_continuous(breaks=c(seq(-1,1,0.2), r_cutoff))+
      ggplot2::ylab("Correlation Strength")+
      ggplot2::xlab(paste0(module_set, " Modules"))+
      ggplot2::labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
                    subtitle = "Inter-Gene PearsonCorrelation")+
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust=0.5),
                     plot.subtitle = ggplot2::element_text(hjust=0.5),
                     axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))

    labels_df_P<-data.frame(name=rep("negLogP", length(p_cutoff)),
                            h_pos = -log10(p_cutoff)-0.2,
                            h_label = paste("p=",p_cutoff,sep=""))

    coherence_boxplot_p<-
      SubGeneCorDF%>%
      dplyr::mutate(P.est = ifelse(P==0, sort(unique(SubGeneCorDF$P))[2], P)) %>% #Fill in true 0 with lowest P-value in dataset
      ggplot2::ggplot(ggplot2::aes(y=-log10(P.est), x=Set))+
      ggplot2::geom_hline(yintercept =0, color="black")+
      ggplot2::geom_boxplot(outlier.shape = NA)+
      ggplot2::geom_hline(yintercept = -log10(p_cutoff), color="red", linetype="dashed")+
      ggplot2::geom_text(data=labels_df_P, ggplot2::aes(label = h_label, x=3, y=h_pos), color="red")+
      ggplot2::ylab("-log10(Correlation P Value)")+
      ggplot2::xlab(paste0(module_set, " Modules"))+
      ggplot2::labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
                    subtitle = "Inter-Gene PearsonCorrelation")+
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust=0.5),
                     plot.subtitle = ggplot2::element_text(hjust=0.5),
                     axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))

    # Or with facets
    labels_df<-data.frame(name=c(rep("Cor", length(r_cutoff)), rep("negLogP", length(p_cutoff))),
                          h_line = c(r_cutoff, -log10(p_cutoff)),
                          h_pos = c(r_cutoff-0.02, -log10(p_cutoff)-0.2),
                          h_label = c(paste("r=",r_cutoff,sep=""), paste("p=",p_cutoff,sep="")))

    coherence_boxplot_faceted<-
      SubGeneCorDF%>%
      dplyr::mutate(P.est = ifelse(P==0, sort(unique(SubGeneCorDF$P))[2], P)) %>% #Fill in true 0 with lowest P-value in dataset
      dplyr::mutate(negLogP = -log10(P.est))%>%
      tidyr::pivot_longer(cols = c(Cor, negLogP))%>%
      ggplot2::ggplot(ggplot2::aes(y=value, x=Set))+
      ggplot2::geom_hline(yintercept =0, color="black")+
      ggplot2::geom_boxplot(outlier.shape = NA)+
      ggplot2::geom_hline(data=labels_df, ggplot2::aes(yintercept = h_line), color = "red", linetype="dashed")+
      ggplot2::geom_text(data=labels_df, ggplot2::aes(label = h_label, x=3, y=h_pos), color="red")+
      ggplot2::xlab(paste0(module_set, " Modules"))+
      ggplot2::labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
                    subtitle = "Inter-Gene PearsonCorrelation")+
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust=0.5),
        plot.subtitle = ggplot2::element_text(hjust=0.5),
        axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))+
      ggplot2::facet_wrap(~name, scales = "free_y", strip.position = "left",
                          labeller = ggplot2::labeller(name = c("Cor" = "Correlation Strength", "negLogP" = "-log10(Correlation P Value)")))


    if(return_plot){print(coherence_boxplot_faceted)}

    if(length(unlist(GSabsent))>0){
      return(list(coherence_boxplot_combined = coherence_boxplot_faceted,
                  coherence_boxplot_cor = coherence_boxplot_cor,
                  coherence_boxplot_p = coherence_boxplot_p,
                  subgene_correlation_df = SubGeneCorDF,
                  GSabsent = GSabsent))} else

                    return(list(coherence_boxplot_combined = coherence_boxplot_faceted,
                                coherence_boxplot_cor = coherence_boxplot_cor,
                                coherence_boxplot_p = coherence_boxplot_p,
                                subgene_correlation_df = SubGeneCorDF))

  }
