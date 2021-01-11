#############################################################
# Master Equation Related Functions
#############################################################
#' sample EVF from tree
SampleSubtree <- function(par,depth,anc_state,edges,ncells,neutral=NA){
  # depth equals to 0 for root cell
  # get the children of the current node
  children <- edges[edges[,2]==par,3]
  result<-lapply(c(1:length(children)),function(j){
    # given the parent and child, find the edge, of the form index, start, end, length
    edge<-edges[edges[,2]==par & edges[,3]==children[j],]
    # if the current node is a leaf, no edge going out of that node
    if(sum(edges[,2]==children[j])==0){
      if(is.na(neutral[1])){
        result <- SampleEdge(edge,depth,anc_state,edges,ncells)}else{
          t_sample <- neutral[neutral[,1]==edge[2] & neutral[,2]==edge[3],3]
          result <- SampleEdge(edge,depth,anc_state,edges,ncells,t_sample)
        }
      result <- result[c(1:(length(result[,1]-1))),]
    }else{
      # if current node is not a leaf, still edges here
      if(is.na(neutral[1])){
        result <- SampleEdge(edge,depth,anc_state,edges,ncells)}else{
          t_sample <- neutral[neutral[,1]==edge[2] & neutral[,2]==edge[3],3]
          result <- SampleEdge(edge,depth,anc_state,edges,ncells,t_sample)
        }
      anc_state <- result[length(result[,1]),4]
      result <- result[c(1:(length(result[,1]-1))),]
      depth <- depth + edge[4]
      # recursively call the sample subtree
      result1 <- SampleSubtree(children[j],depth,anc_state,edges,ncells,neutral)
      result <- rbind(result,result1)
    }
    return(result)
  })
  result<-do.call(rbind,result)
  return(result)
}

#' sample EVF from edge
SampleEdge <- function(edge,depth,anc_state,edges,ncells,t_sample=NA){
  if(is.na(t_sample[1])){
    #t_sample <- c(0,sort( runif(round(edge[4]*ncells/sum(edges[,4])),0,edge[4]) ))
    if (ceiling(edge[4]*ncells/sum(edges[,4]))-1 < 0) {stop("the total number of cells is too few.")}
    # ceiling(edge[4]*ncells/sum(edges[,4]))-1) -> number of cell included in the current edge
    # edge[4]/ceiling(edge[4]*ncells/sum(edges[,4]))-1) -> intervel size
    t_sample <- c(0,seq(0, edge[4], edge[4]/(ceiling(edge[4]*ncells/sum(edges[,4]))-1)))
    # add in 0 and edge[4] at the beginning and the end
    t_sample<-c(t_sample,edge[4])
  }else{
    t_sample<-sort(c(0,t_sample-depth))
  }

  # calculate the difference between every two adj element in an array e.g. a = [1,2,3] diff(a) = [1,1]
  t_interval<-diff(t_sample)
  # change as gaussian random variable with variance related to the interval size
  x_change <- sapply(t_interval,function(sig){rnorm(1,0,sqrt(sig))})
  # for each element in x_change, sum over current and all previous elements, input [1,2,3], output [1.3.6]
  x_sample <- cumsum(x_change)
  # add all x_sample with anc_state, t_sample[-1]: t_sample remove the first 0 element
  # different depths correspond to different cells, anc_state was evf center, add in new move
  result<-cbind(depth+t_sample[-1],x_sample+anc_state)
  # result of the form start-edge, end-edge, depth, evf
  result <- cbind(rep(edge[2],length(result[,1])),rep(edge[3],length(result[,1])),result)
  return(result)
}

#' sample from truncated normal distribution
#' @param a the minimum value allowed
#' @param b the maximum value allowed
#' @return samples
#' @export
rnorm_trunc <- function(n, mean, sd, a, b){
  vec1 <- rnorm(n, mean = mean, sd=sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i){
      while (TRUE){
        temp <- rnorm(1, mean = mean, sd=sd)
        if (temp > a & temp < b) {break}}
      return(temp)} )
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}

#' Getting GeneEffects matrices
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from
#' @param evf_res the EVFs generated for cells
#' @param is_in_module a vector of length ngenes. 0 means the gene is not in any gene module. In case of non-zero values, genes with the same value in this vector are in the same module.
#' @return a list of 3 matrices, each of dimension ngenes * nevf
#' @export
GeneEffects <- function(ngenes,nevf,randseed,prob,geffect_mean,geffect_sd, evf_res, is_in_module){
  set.seed(randseed)

  gene_effects <- lapply(c('kon','koff','s'),function(param){
    effect <- lapply(c(1:ngenes),function(i){
      nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
      nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=geffect_mean,sd=geffect_sd)
      return(nonzero)
    })
    return(do.call(rbind,effect))
  })


  mod_strength <- 0.8
  if (sum(is_in_module) > 0){ # some genes need to be assigned as module genes; their gene_effects for s will be updated accordingly
    for (pop in unique(is_in_module[is_in_module>0])){ # for every population, find the
      #which(evf_res[[2]]$pop == pop)
      for (iparam in c(1,3)){
        nonzero_pos <- order(colMeans(evf_res[[1]][[iparam]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[iparam]]),
                             decreasing = T)[1:ceiling(nevf*prob)]
        geffects_centric <- rnorm(length(nonzero_pos), mean = 0.5, sd=0.5)
        for (igene in which(is_in_module==pop)){
          for (ipos in 1:length(nonzero_pos)){
            gene_effects[[iparam]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
          }
          gene_effects[[iparam]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
        }
      }

      nonzero_pos <- order(colMeans(evf_res[[1]][[2]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[2]]),
                           decreasing = T)[1:ceiling(nevf*prob)]
      geffects_centric <- rnorm(length(nonzero_pos), mean = -0.5, sd=0.5)
      for (igene in which(is_in_module==pop)){
        for (ipos in 1:length(nonzero_pos)){
          gene_effects[[2]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
        }
        gene_effects[[2]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
      }
    }
  }
  return(gene_effects)
}

#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample,den_fun){
  probs <- den_fun$y/sum(den_fun$y)
  bw <- den_fun$x[2]-den_fun$x[1]
  bin_id <- sample(size=nsample,x=c(1:length(probs)),prob=probs,replace=T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)),function(j){
    runif(n=counts[j],min=(den_fun$x[sampled_bins[j]]-0.5*bw),max=(den_fun$x[sampled_bins[j]]+0.5*bw))
  })
  samples <- do.call(c,samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows.
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from
#' @param bimod the bimodality constant
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @import plyr
#' @return params a matrix of ngenes * 3
#' @export
Get_params <- function(gene_effects,evf,match_param_den,bimod,scale_s){
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    # X=matrix(data=c(1:10),ncol=2)
    # this line is to check that the row and columns did not flip
    temp <- alply(X, 1, function(Y){Y})
    values <- do.call(c,temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample=max(ranks),den_fun=match_param_den[[i]]))
    temp3 <- matrix(data=sorted[ranks],ncol=length(X[1,]),byrow=T)
    return(temp3)
  })

  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x}))*scale_s

  return(scaled_params)
}


#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function),
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows.
#' @param param_realdata the fitted parameter distribution to sample from
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @export
Get_params2 <- function(gene_effects,evf,bimod,ranges){
  print("calculate paras")
  # params <- lapply(gene_effects,function(X){evf %*% t(X)})
  # updated a little bit
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})

  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    temp <- apply(X,2,function(x){1/(1+exp(-x))})
    temp2 <- temp*(ranges[[i]][2]-ranges[[i]][1])+ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
}

#' Generating EVFs for cells sampled along the trajectory of cell development
#' @param phyla tree for cell developement
#' @param ncells number of cells
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
#' @import phytools
#' @export
ContinuousEVF <- function(phyla,ncells,n_nd_evf,n_de_evf,evf_center=1,vary='s',
                          Sigma,seed){
  set.seed(seed)
  edges <- cbind(phyla$edge,phyla$edge.length)
  edges <- cbind(c(1:length(edges[,1])),edges)
  # number of edges that connect to the current node
  connections <- table(c(edges[,2],edges[,3]))

  # ROOT: no incoming edges
  root <- as.numeric(which(lapply(seq(1,length(connections)), function(x){is.na(match(x, edges[,3]))}) == TRUE))
  #TIPS: no outgoing edges
  tips <- as.numeric(which(lapply(seq(1,length(connections)), function(x){is.na(match(x, edges[,2]))}) == TRUE))
  # INTERNAL: other
  internal <- as.numeric(seq(1, length(connections))[-c(root,tips)])


  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }
  # evf_center mean of evf of non-Diff EVFs.
  neutral <- SampleSubtree(root,0,evf_center,edges,ncells)
  param_names <- c("kon", "koff", "s")
  # sample nd_evf for k_on, k_off and s, the nd_evf value is shared for all cells, with center and variance given by user
  # parami equals to k_on, k_off and s
  evfs <- lapply(c(1:3),function(parami){

    # sample nd_evf
    nd_evf <- lapply(c(1:N_ND_evfs[parami]),function(ievf){
      rnorm(ncells,evf_center,Sigma)
    })
    # ncells by nd_evf
    nd_evf <- do.call(cbind,nd_evf)

    # sample de_evf
    if(N_DE_evfs[parami]!=0){
      #if there is more than 1 de_evfs for the parameter we are looking at
      de_evf <- lapply(c(1:N_DE_evfs[parami]),function(evf_i){
        SampleSubtree(root,0,evf_center,edges,ncells,neutral=neutral)
      })

      de_evf <- lapply(de_evf,function(X){X[,4]})
      de_evf <- do.call(cbind,de_evf)
      de_evf <- de_evf[c(1:ncells),] # the de_evf for more cells are generated. Here we need to truncate the last ones.
      evfs <- cbind(nd_evf,de_evf)
      colnames(evfs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(de_evf[1,])),c(1:length(de_evf[1,])),sep='_'))
    }else{
      evfs <- nd_evf
      colnames(evfs)<-paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_')
    }
    return(evfs)
  })

  names(evfs) <- c("k_on", "k_off", "s")

  meta <- data.frame(pop=apply(neutral[,c(1:2)],1,function(X){paste0(X,collapse='_')}),depth=neutral[,3])
  return(list(evfs,meta[c(1:ncells),]))
  # note previously the number of sampled evfs and meta isn't necessarily ncells?
}


#' Generating first half cycle of EVFs for cells sampled along cycle trajectory
#' @param ncells Number of cells
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param depth depth of the sampling process
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
GenerateForwardEVF <- function(ncells,n_nd_evf,n_de_evf, depth=1,
                               evf_center=1,vary='s',Sigma,seed){
  set.seed(seed)

  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }

  param_names <- c("kon", "koff", "s")

  evfs <- lapply(c(1:3),function(parami){

    # sample nd_evf
    nd_evf <- lapply(c(1:N_ND_evfs[parami]),function(ievf){
      rnorm(ncells,evf_center,Sigma)
    })
    # make matrix of dimension (ncells, N_ND_evfs[parami])
    nd_evf <- do.call(cbind,nd_evf)


    # sample de_evf
    if(N_DE_evfs[parami]!=0){
      #if there is more than 1 de_evfs for the parameter we are looking at
      de_evf <- lapply(c(1:N_DE_evfs[parami]),function(evf_i){
        # random time step
        # t_sample <- c(0,sort( runif(ncells-2,0,depth),depth))
        # uniform time step, of length ncells + 1
        t_sample <- seq(0, depth, depth/ncells)

        # interval of length ncells
        t_interval<-diff(t_sample)
        # change as gaussian random variable with variance related to the interval size
        x_change <- sapply(t_interval,function(sig){rnorm(1,0,sqrt(sig))})
        # for each element in x_change, sum over current and all previous elements, input [1,2,3], output [1.3.6]
        x_sample <- cumsum(x_change)

        # evf_center is the root cell evf, need to remove the first value of t_sample, which correspond to evf_center directly
        # depth, evf, of length ncells
        result <- cbind(t_sample[-1],x_sample+evf_center)
        return(result)
      })

      # extract de_evf value, and make a matrix, of dimension ncells by n_de_evf
      de_evf <- lapply(de_evf,function(X){X[,2]})
      # should be of dimension ncells by n_de_evf
      de_evf <- do.call(cbind,de_evf)
      # total evfs
      evfs <- cbind(nd_evf,de_evf)
      colnames(evfs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(de_evf[1,])),c(1:length(de_evf[1,])),sep='_'))
    }else{
      evfs <- nd_evf
      colnames(evfs)<-paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_')
    }

    return(evfs)
  })
  return(list(evfs, as.data.frame(list(pop = rep("cycle", ncells), depth = rep(0.0, ncells)))))
}

#' Generating second half cycle of EVFs for cells sampled along cycle trajectory
#' @param evf_res first half cycle of EVFs
#' @param n_nextra Number of cells in second cycle
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
GenerateBackwardEVF <- function(evf_res, n_nextra, n_nd_evf,n_de_evf,
                                evf_center=1,vary='s',Sigma,seed){

  set.seed(seed)
  kinetic_evf <- evf_res[[1]]
  backbone_info <- evf_res[[2]]
  param_names <- c("k_on", "k_off", "s")
  names(kinetic_evf) <- param_names

  # number of cells for the backward process
  ncells <- dim(kinetic_evf[["k_on"]])[[1]] - 1

  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }

  # random ordering the backward stepsize
  shuffled_order <- sample(seq(1, ncells))

  evfs <- lapply(c(1:3),function(parami){

    # sample nd_evf
    nd_evf <- lapply(c(1:N_ND_evfs[parami]),function(ievf){
      rnorm(ncells,evf_center,Sigma)
    })
    # make matrix of dimension (ncells, N_ND_evfs[parami])
    nd_evf <- do.call(cbind,nd_evf)


    # sample de_evf
    if(N_DE_evfs[parami]!=0){

      # get the forward evf values the correspond to de
      diff_evf <- kinetic_evf[[parami]][2:(ncells+1),(n_nd_evf+1):(n_nd_evf + n_de_evf)] - kinetic_evf[[parami]][1:ncells,(n_nd_evf+1):(n_nd_evf + n_de_evf)]
      # backward of diff_evf
      diff_evf_back <- -diff_evf[shuffled_order,]
      # start evf of backward process
      start_evf <- kinetic_evf[[parami]][ncells,(n_nd_evf+1):(n_nd_evf + n_de_evf)]
      # accumulate the diff value
      cum_evf_back <- apply(diff_evf_back, 2, cumsum)
      # all row plus the starting evf, calculate the current evf
      evf_back <- sweep(cum_evf_back, 2, start_evf, "+")
      # bind nd evf with de evf to form the complete evf
      evfs <- cbind(nd_evf,evf_back)

      # name the evfs accordingly
      colnames(evfs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(evf_back[1,])),c(1:length(evf_back[1,])),sep='_'))

    }else{
      evfs <- nd_evf
      colnames(evfs)<-paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_')
    }

    if(n_nextra > 0){
      # truncate the backbward cells to make the number of cells sum up to ncells_total
      evfs <- rbind(kinetic_evf[[parami]], evfs[1:(ncells-1),], kinetic_evf[[parami]][1:n_nextra,])
    }else{
      evfs <- rbind(kinetic_evf[[parami]], evfs[1:(ncells-1),])
    }

    return(evfs)
  })

  names(evfs) <- c("k_on", "k_off", "s")

  return(list(evfs, as.data.frame(list(pop = rep("cycle", dim(evfs[["k_on"]])[1]), depth = rep(0.0, dim(evfs[["k_on"]])[1])))))

}

#' Generating EVFs for cells sampled along cycle trajectory
#' @param ncells_total number of cells
#' @param nextra number of cells in second cycle
#' @param n_unstable number of cells before reaching steady state
#' @param n_evf Number of EVFs
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param depth depth of the sampling process
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
#' @export
GenerateCycleEVF <- function(ncells_total, nextra, n_unstable, nevf,
                             n_de_evf, vary = "all", depth = 1, evf_center = 1, Sigma, seed){
  # when reach this number of cells, start going back
  turning_cell <- (ncells_total - nextra - n_unstable)/2+1

  # calculate the evf value for each cell before (include) turning cell, and also include n_unstable
  evf_res <-  GenerateForwardEVF(ncells = turning_cell+n_unstable,n_nd_evf = nevf-n_de_evf,
                                 n_de_evf=n_de_evf, depth=depth, evf_center=evf_center,vary=vary,
                                 Sigma=Sigma,seed=seed[1])

  # remove the evf of unstable cell here
  evf_res[[1]] <- lapply(evf_res[[1]], function(i){
    return(i[(n_unstable+1):(n_unstable+turning_cell),])
  })
  evf_res[[2]] <- evf_res[[2]][(n_unstable+1):(n_unstable+turning_cell),]

  # evf_res is a list of two elements, the first element is a list correspond to the kinetic evf, with three element correspond
  # to the kon, koff and s. The second element correspond to the pop and depth in the backbone.
  evf_res <- GenerateBackwardEVF(evf_res=evf_res, n_nextra=nextra, n_nd_evf=nevf-n_de_evf,
                                 n_de_evf=n_de_evf, evf_center=evf_center,vary=vary,
                                 Sigma=Sigma,seed=seed[1])

  return(evf_res)
}


