#' Generate true transcript counts for linear structure
#' @param kinet_params kinetic parameters, include k_on, k_off, s and beta
#' @param start_state the starting state: on or off of each gene
#' @param start_s spliced count of the root cell in the branch
#' @param start_u unspliced count of the root cell in the branch
#' @param randpoints1 the value which evf mean is generated from
#' @param ncells1 number of cells in the branch
#' @param ngenes number of genes
#' @param beta_vec splicing rate of each gene
#' @param d_vec degradation rate of each gene
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
#' @export
gen_1branch <- function(kinet_params, start_state, start_s, start_u, randpoints1, ncells1, ngenes, beta_vec, d_vec){

  # totaltime equals to total numeber of cells
  totaltime <- ncells1

  # store unspliced count
  counts_u1 <- matrix(0, ngenes, ncells1)
  # store spliced count
  counts_s1 <- matrix(0, ngenes, ncells1)
  # store state matrix
  state_mat <- matrix(0, ngenes, ncells1)

  # store the kinetic values, include k_on, k_off, s and beta

  kinet_params1 <- list(k_on=kinet_params$k_on, k_off=kinet_params$k_off, s=kinet_params$s,
                        beta=beta_vec, d=d_vec)
  xs <- list()
  ys <- list()
  which_cells <- list()

  for (igene in 1:ngenes){
    k_on <- kinet_params1[["k_on"]][igene,1]; k_off <- kinet_params1[["k_off"]][igene,1];
    s <- kinet_params1[["s"]][igene,1]; beta <- beta_vec[igene]; d <- d_vec[igene];

    cycle_length <- 1/k_on + 1/k_off
    min_wtime <- min(1/k_on, 1/k_off)
    npart <- max(ceiling(cycle_length/min_wtime)*2, cycle_length*2)
    stepsize <- cycle_length/npart

    nsteps <- ceiling(totaltime/stepsize)
    x <- numeric(nsteps); x[1] <- start_u[igene]
    y <- numeric(nsteps); y[1] <- start_s[igene]
    # which_cell stores with length nsteps stores the cell correspond to current time step.
    which_cell <- numeric(nsteps); which_cell[1] <- 1
    curr_time <- numeric(nsteps);
    p_table <- matrix(0, 2, 2)
    p_table[1,2] <- 1/(npart*(1/k_off/cycle_length)); p_table[1,1] <- 1-p_table[1,2];
    p_table[2,1] <- 1/(npart*(1/k_on/cycle_length)); p_table[2,2] <- 1-p_table[2,1]

    curr_state <- numeric(nsteps); curr_state[1] <- start_state[igene] # 1 means on state, 2 means off state
    t <- 1
    while (which_cell[t] < ncells1){
      t <- t+1
      curr_time[t] <- curr_time[t-1] + stepsize
      if ( (curr_time[t] - which_cell[t-1]) > 0){
        which_cell[t] <- which_cell[t-1] + 1
        # also evolve the kinetic parameters a little bit
        k_on <- kinet_params1[["k_on"]][igene, which_cell[t]];
        k_off <- kinet_params1[["k_off"]][igene, which_cell[t]];
        s <- kinet_params1[["s"]][igene, which_cell[t]];

        # check npart code here, on p_table, update the step size
        cycle_length <- 1/k_on + 1/k_off
        min_wtime <- min(1/k_on, 1/k_off)
        npart <- max(ceiling(cycle_length/min_wtime)*2, cycle_length*2)
        stepsize <- cycle_length/npart
        p_table <- matrix(0, 2, 2)
        p_table[1,2] <- 1/(npart*(1/k_off/cycle_length)); p_table[1,1] <- 1-p_table[1,2];
        p_table[2,1] <- 1/(npart*(1/k_on/cycle_length)); p_table[2,2] <- 1-p_table[2,1]
      } else {which_cell[t] <- which_cell[t-1]}

      if (runif(1, 0, 1) > p_table[curr_state[t-1],1]){
        curr_state[t] <- 2
      } else {curr_state[t] <- 1}

      if (curr_state[t] == 1) {
        x[t] <- x[t-1] + s*stepsize - beta*x[t-1]*stepsize
        if (x[t] < 0) {x[t] <- 0}
      } else {x[t] <- x[t-1] -beta*x[t-1]*stepsize}
      y[t] <- y[t-1] + beta*x[t-1]*stepsize - d*y[t-1]*stepsize
      if (y[t] < 0) {y[t] <- 0}
    }

    xs[[igene]] <- x; ys[[igene]] <- y; which_cells[[igene]] <- which_cell
    # extract value for each cell
    for (icell in 1:ncells1){
      all_idx <- which(which_cell == icell)
      closest <- which.min(abs((curr_time[all_idx] - (icell-1)) - randpoints1[icell]))
      counts_u1[igene,icell] <- x[all_idx[closest]]
      counts_s1[igene,icell] <- y[all_idx[closest]]
      state_mat[igene,icell] <- curr_state[all_idx[closest]]
    }
  }
  cell_time <- randpoints1+(0:(ncells1-1))

  # calculate the ground truth velocity
  ncells2rm <- 0
  velo_mat <- beta_vec*counts_u1 - d_vec*counts_s1

  dynamics <- list()
  dynamics[["unspliced_steps"]] <- xs
  dynamics[["spliced_steps"]] <- ys
  dynamics[["which_cells"]] <- which_cells

  return(list(counts_u = counts_u1, counts_s=counts_s1, kinet_params=kinet_params1, state_mat=state_mat, cell_time=cell_time, velocity=velo_mat, dynamics=dynamics))
}


#' Update the kinetic parameters of highly expressed genes
#' @param prop_hge proportion of the highly expressed gene
#' @param mean_hge parameter of highly expressed gene
#' @param params kinetic parameters
#' @return updated kinetic parameters
#' @import phytools
#' @export
includeHge <- function(prop_hge=0.015, mean_hge=5, params){
  # number of genes
  ngenes <- dim(params[["k_on"]])[1]
  # number of cells
  ncells_total <- dim(params[["k_on"]])[2]

  if (prop_hge>0){
    # select hge from all possible genes
    chosen_hge <- sample(ngenes, ceiling(ngenes*prop_hge), replace = F)
    # vector of length hge
    multi_factors <- numeric(length(chosen_hge))
    rank_sum <- rank(rowSums(params[[3]][chosen_hge,]))
    multi_factors <- sapply(1:length(chosen_hge), function(igene){
      tosubstract <- -rank_sum[igene]*1/length(chosen_hge)+1
      if (runif(1, 0, 1) < 1){
        multi_factor <- mean_hge - tosubstract
      } else {multi_factor <- mean_hge}
      return(multi_factor)
    })
    new_s <- matrix(0, length(chosen_hge), ncells_total)
    for (i in 1:length(chosen_hge)){
      new_s[i,] <- params[[3]][chosen_hge[i],] * (2^multi_factors[i])
    }
    params[[3]][chosen_hge,] <- new_s # update s

  } else {chosen_hge <- NULL}

  return(list(params = params, chosen_hge = chosen_hge))
}


#' Generate both evf and gene effect and simulate true transcript counts, tree structure
#' @param ncells_total number of cells
#' @param ngenes number of genes
#' @param start_s initial spliced count, NULL if not specified
#' @param start_u initial unspliced count, NULL if not specified
#' @param start_stats initial kinetic state, NULL if not specified
#' @param evf_center the value which evf mean is generated from
#' @param nevf number of evfs
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param randseed random seed
#' @param n_de_evf number of differential evfs between populations
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param Sigma parameter of the std of evf values within the same population
#' @param geffect_mean the mean of the normal distribution where the non-zero gene effect sizes are sampled from
#' @param gene_effect_sd the standard deviation of the normal distribution where the non-zero gene effect sizes are sampled from
#' @param gene_effect_prob the probability that the effect size is not 0
#' @param bimod the amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal
#' @param param_realdata pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param prop_hge the proportion of very highly expressed genes
#' @param mean_hge the parameter to amplify the gene-expression levels of the very highly expressed genes
#' @param n_unstable the number of cells to be removed at the beginning of the simulation, sutract from ncells_total
#' @param  plot plot the kinetic or not
#' @return a list of 6 elements: unspliced true counts, spliced true counts, kinetic parameter, state matrix, cell pseudo-time, velocity
#' @import phytools
#' @export
SimulateVeloTree <- function(ncells_total, ngenes, start_s = NULL, start_u = NULL, start_state = NULL, evf_center=1,nevf=20, phyla,
                             randseed, n_de_evf=12,vary='s',Sigma=0.1,
                             geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                             bimod=0,param_realdata="zeisel.imputed",scale_s=1,
                             prop_hge=0.015, mean_hge=5, n_unstable=25, plot=FALSE){
  set.seed(randseed)
  seed <- sample(c(1:1e5),size=2)



  # calculate the evf value for each cell
  evf_res <- ContinuousEVF(phyla,ncells_total,n_nd_evf=nevf-n_de_evf,n_de_evf=n_de_evf,
                           evf_center=evf_center,vary=vary,
                           Sigma=Sigma,seed=seed[1])

  if(plot == TRUE){
    plotEVF(evf_res, width = 15, height = 15, units = "in", dpi = 1000)
  }
  # assign modules automatically and return module idx for genes.
  is_in_module <- numeric(ngenes)

  # generate the gene effect matrix
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd, evf_res=evf_res, is_in_module=is_in_module)

  # match distribution and generate kinetic values: kon, koff, s, beta, d
  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      data(param_realdata.zeisel.imputed)
    } else {stop("wrong input for parameter param_realdata")}

    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params[,3]=log(base=10,match_params[,3])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod,scale_s=scale_s)
  }else{
    params <- Get_params2(gene_effects,evf_res[[1]],bimod,ranges)
  }
  # kinetic params
  names(params) <- c("k_on", "k_off", "s")

  # update the kinetic value of the highly expressed genes
  hge <- includeHge(prop_hge=prop_hge, mean_hge=mean_hge, params)
  chosen_hge <- hge[["chosen_hge"]]
  params <- hge[["params"]]

  # degradation rate of the size ngenes * 1
  d_vec <- rnorm(n=ngenes, mean=0.4, sd=0.1)
  # splicing rate of the size ngenes * 1
  beta_vec <- rnorm(n=ngenes, mean=0.8, sd=0.1)


  # process tree file
  # number of edges that connect to the current node
  degrees <- table(c(phyla$edge[,1],phyla$edge[,2]))

  # root node is the node with no incoming edges
  root <- as.numeric(which(lapply(seq(1,length(degrees)), function(x){is.na(match(x, phyla$edge[,2]))}) == TRUE))

  # find starting and ending cells for each branches
  end_cell <- vector()
  branches <- vector()
  num_cells <- vector()

  for(pop in unique(evf_res[[2]]$pop)){
    # find all cell indices correspond to the current population
    index <- which(evf_res[[2]]$pop == pop)
    # the index is ordered according to depth originally
    end_cell <- c(end_cell, index[length(index)])
    new_branch <- as.numeric(strsplit(pop, split ="_")[[1]])
    branches <- rbind(branches, new_branch)
    num_cells <- c(num_cells,length(index))
  }
  # create a backbone_info matrix
  backbone <- cbind(start_node=branches[,1], end_node = branches[,2], end_cell = end_cell, start_cell = rep(0, length(end_cell)), num_cells = num_cells)
  # print(backbone)

  for(i in 1:dim(backbone)[1]){
    if(backbone[i,"start_node"] == root){
      backbone[i,"start_cell"] <- 0
    }else{
      index <- which(backbone[,"end_node"] == backbone[i, "start_node"])
      backbone[i,"start_cell"] <- backbone[index,"end_cell"]
    }
  }

  # print(backbone)
  # generate counts
  # unspliced counts matrix
  counts_u <- vector()
  # spliced counts matrix
  counts_s <- vector()
  # cell developmental time
  cell_time <- vector()
  # count step changes
  steps_u <- vector()
  steps_s <- vector()
  which_cells <- vector()
  # promoter states
  state_mat <- vector()
  # RNA velocity
  velocity <- vector()

  # generate the promoter states of the root cell with dimension ngenes * 1
  root_state <- sample(c(1,2), size = ngenes, replace = TRUE)

  cell2rm <- vector()

  # loop through all branches
  for(i in 1:dim(backbone)[1]){
    kinetic_params <- list()
    # indices of the cells that belongs to current backbone
    index <- which(evf_res[[2]]$pop == paste0(backbone[i, "start_node"], "_", backbone[i, "end_node"]))
    # subset k_on, k_off and s

    for(j in names(params)){
      # dimension ngenes by #index(subset of cells)
      kinetic_params[[j]] <- params[[j]][,index]
    }
    # the points we use to get snapshot for each cell
    randpoints <- runif(n=backbone[i, "num_cells"], min=0, max=1)

    if(backbone[i, "start_node"] == root){
      if(is.null(start_s) | is.null(start_u) | is.null(start_state)){
        start_s <- rep(0, ngenes)
        start_u <- rep(0,ngenes)
        start_state <- root_state
      }
      # if use dim, then the starting would be NA
      cells2rm <- c(cell2rm, seq(as.integer(length(counts_s)/ngenes) + 1,as.integer(length(counts_s)/ngenes) + n_unstable))
      start_cell_time <- 0
    }else{
      # start cell must have already been generated, and note that the index must be larger than 1
      start_s <- counts_s[,backbone[i,"start_cell"]]
      start_u <- counts_u[,backbone[i,"start_cell"]]
      start_state <- state_mat[,backbone[i,"start_cell"]]
      start_cell_time <- cell_time[backbone[i,"start_cell"]]
    }

    result <- gen_1branch(kinet_params = kinetic_params, start_state = start_state, start_u = start_u, start_s = start_s,
                          randpoints1 = randpoints, ncells1 = backbone[i, "num_cells"], ngenes = ngenes, beta_vec = beta_vec, d_vec = d_vec)

    counts_u <- cbind(counts_u, result$counts_u)
    counts_s <- cbind(counts_s, result$counts_s)
    state_mat <- cbind(state_mat, result$state_mat)

    steps_s <- cbind(steps_s, result$dynamics[["spliced_steps"]])
    steps_u <- cbind(steps_u, result$dynamics[["unspliced_steps"]])
    which_cells <- cbind(which_cells, result$dynamics[["which_cells"]])

    # true time of each cell, start from zero for each branch, need to update when merged together
    cell_time <- c(cell_time, result$cell_time + start_cell_time)
    velocity <- cbind(velocity, result$velocity)
  }

  counts_u <- counts_u[,-cells2rm]
  counts_s <- counts_s[,-cells2rm]
  state_mat <- state_mat[,-cells2rm]
  cell_time <- cell_time[-cells2rm]
  velocity <- velocity[,-cells2rm]
  final_kinetics <- list(k_on = params[[1]][,-cells2rm], k_off = params[[2]][,-cells2rm],
                         s = params[[3]][,-cells2rm], beta = beta_vec, degrade = d_vec)

  dynamics <- list(unspliced_steps = steps_u, spliced_steps = steps_s, which_cells = which_cells)

  return(list(counts_u = counts_u, counts_s=counts_s, kinet_params=final_kinetics, state_mat=state_mat,
              cell_time=cell_time, velocity=velocity, dynamics=dynamics,
              backbone = evf_res[[2]]$pop[-cells2rm], chosen_hge = chosen_hge))
}


#' Generate both evf and gene effect and simulate true transcript counts, cycle structure
#' @param ncells_total number of cells
#' @param ngenes number of genes
#' @param start_s initial spliced count, NULL if not specified
#' @param start_u initial unspliced count, NULL if not specified
#' @param start_stats initial kinetic state, NULL if not specified
#' @param evf_center the value which evf mean is generated from
#' @param nevf number of evfs
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param randseed random seed
#' @param n_de_evf number of differential evfs between populations
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param Sigma parameter of the std of evf values within the same population
#' @param geffect_mean the mean of the normal distribution where the non-zero gene effect sizes are sampled from
#' @param gene_effect_sd the standard deviation of the normal distribution where the non-zero gene effect sizes are sampled from
#' @param gene_effect_prob the probability that the effect size is not 0
#' @param bimod the amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal
#' @param param_realdata pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param prop_hge the proportion of very highly expressed genes
#' @param mean_hge the parameter to amplify the gene-expression levels of the very highly expressed genes
#' @param n_unstable the number of cells to be removed at the beginning of the simulation, sutract from ncells_total
#' @param  plot plot the kinetic or not
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
#' @export
SimulateVeloCycle <- function(ncells_total, ngenes, start_s = NULL, start_u = NULL, start_state = NULL, evf_center=1,nevf=20,
                              randseed, n_de_evf=12,vary='all',Sigma=0.1,
                              geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                              bimod=0,param_realdata="zeisel.imputed",scale_s=1,
                              prop_hge=0.015, mean_hge=5, nextra=0, n_unstable=25,
                              plot = FALSE){
  if(nextra < 0 | nextra >= (ncells_total-n_unstable)/2 | n_unstable < 0 | n_unstable >= ncells_total){
    stop("nextra should be selected from [0, (ncells_total-n_unstable)/2), and n_unstable from [0, ncells_total)")
  }

  set.seed(randseed)
  seed <- sample(c(1:1e5),size=2)


  evf_res <- GenerateCycleEVF(ncells_total, nextra, n_unstable, nevf,
                              n_de_evf, vary = vary, depth = 10, evf_center = evf_center, Sigma = Sigma, seed = seed[1])

  if(plot == TRUE){
    plotEVF(evf_res, width = 15, height = 15, units = "in", dpi = 1000)
  }

  # assign modules automatically and return module idx for genes.
  is_in_module <- numeric(ngenes)

  # generate the gene effect matrix
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd, evf_res=evf_res, is_in_module=is_in_module)

  # match distribution and generate kinetic values: kon, koff, s, beta, d
  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      data(param_realdata.zeisel.imputed)
    } else {stop("wrong input for parameter param_realdata")}

    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params[,3]=log(base=10,match_params[,3])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod,scale_s=scale_s)
  }else{
    params <- Get_params2(gene_effects,evf_res[[1]],bimod,ranges)
  }
  # kinetic params
  names(params) <- c("k_on", "k_off", "s")

  # update the kinetic value of the highly expressed genes
  hge <- includeHge(prop_hge=prop_hge, mean_hge=mean_hge, params)
  chosen_hge <- hge[["chosen_hge"]]
  params <- hge[["params"]]

  # degradation rate of the size ngenes * 1
  d_vec <- rnorm(n=ngenes, mean=0.8, sd=0.1)
  # splicing rate of the size ngenes * 1
  beta_vec <- rnorm(n=ngenes, mean=0.8, sd=0.1)

  # if the start state is not provided
  if(is.null(start_s) | is.null(start_u) | is.null(start_state)){
    root_state <- sample(c(1,2), size = ngenes, replace = TRUE)
    start_s <- rep(0, ngenes)
    start_u <- rep(0,ngenes)
    start_state <- root_state
  }

  start_cell_time <- 0
  # remove unstable cells
  ncells_total <- ncells_total -n_unstable
  randpoints <- runif(n=ncells_total, min=0, max=1)

  result <- gen_1branch(kinet_params = params, start_state = start_state, start_s = start_s, start_u = start_u,
                        randpoints1 = randpoints, ncells1 = ncells_total, ngenes = ngenes, beta_vec = beta_vec,
                        d_vec = d_vec)

  final_kinetics <- list(k_on = result$kinet_params$k_on, k_off = result$kinet_params$k_off,
                         s = result$kinet_params$s, beta = result$kinet_params$beta, degrade = result$kinet_params$d)

  dynamics <- list(unspliced_steps = result$dynamics[["unspliced_steps"]], spliced_steps = result$dynamics[["spliced_steps"]], which_cells = result$dynamics[["which_cells"]])

  return(list(counts_u = result$counts_u, counts_s=result$counts_u, kinet_params=final_kinetics, state_mat=result$state_mat,
              cell_time=result$cell_time, velocity=result$velocity, dynamics=dynamics,
              backbone = evf_res[[2]]$pop, chosen_hge = chosen_hge))
}

#' Generate both evf and gene effect and simulate true transcript counts, cycle-tree structure
#' @param ncells_total number of cells
#' @param ncells_cycle number of cells in cell-cycle stage
#' @param ngenes number of genes
#' @param evf_center the value which evf mean is generated from
#' @param nevf number of evfs
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param randseed random seed
#' @param n_de_evf number of differential evfs between populations
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param Sigma parameter of the std of evf values within the same population
#' @param geffect_mean the mean of the normal distribution where the non-zero gene effect sizes are sampled from
#' @param gene_effect_sd the standard deviation of the normal distribution where the non-zero gene effect sizes are sampled from
#' @param gene_effect_prob the probability that the effect size is not 0
#' @param bimod the amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal
#' @param param_realdata pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param prop_hge the proportion of very highly expressed genes
#' @param mean_hge the parameter to amplify the gene-expression levels of the very highly expressed genes
#' @param n_unstable the number of cells to be removed at the beginning of the simulation, sutract from ncells_total
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
#' @export
SimulateCycleTree <- function(ncells_total, ncells_cycle, ngenes, evf_center=1, nevf=20,
                              phyla, randseed, n_de_evf=12, vary="s", Sigma=0.1,
                              geffect_mean=0, gene_effects_sd=1, gene_effect_prob=0.3,
                              bimod=0, param_realdata="zeisel.imputed",scale_s=1,
                              prop_hge=0.015, mean_hge=5, nextra=0, n_unstable=25, plot = TRUE,
                              start_s = NULL, start_u = NULL, start_state = NULL, is_in_module){

  if(nextra < 0 | nextra >= (ncells_cycle-n_unstable)/2 | n_unstable < 0 | n_unstable >= ncells_cycle){
    stop("nextra should be selected from [0, (ncells_total-n_unstable)/2), and n_unstable from [0, ncells_total)")
  }

  set.seed(randseed)
  seed <- sample(c(1:1e5),size=2)

  evf_res_cycle <- GenerateCycleEVF(ncells_cycle, nextra, n_unstable, nevf,
                              n_de_evf, vary = vary, depth = 10, evf_center = evf_center, Sigma = Sigma, seed = seed[1])

  evf_res_tree <- ContinuousEVF(phyla,ncells_total-ncells_cycle,n_nd_evf=nevf-n_de_evf,n_de_evf=n_de_evf,
                           evf_center=evf_center,vary=vary,
                           Sigma=Sigma,seed=seed[1])

  ncells_cycle = dim(evf_res_cycle[[1]]$k_on)[1]
  ncells_tree = dim(evf_res_tree[[1]]$k_on)[1]

  evf_res <- list(list(k_on = rbind(evf_res_cycle[[1]]$k_on, evf_res_tree[[1]]$k_on), k_off = rbind(evf_res_cycle[[1]]$k_off, evf_res_tree[[1]]$k_off),
                  s = rbind(evf_res_cycle[[1]]$s, evf_res_tree[[1]]$s)), rbind(evf_res_cycle[[2]], evf_res_tree[[2]]))

  if(plot == TRUE){
    plotEVF(evf_res, width = 15, height = 15, units = "in", dpi = 1000)
  }

  # assign modules automatically and return module idx for genes.
  is_in_module <- numeric(ngenes)

  # generate the gene effect matrix
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd, evf_res=evf_res, is_in_module=is_in_module)

  # match distribution and generate kinetic values: kon, koff, s, beta, d
  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      data(param_realdata.zeisel.imputed)
    } else {stop("wrong input for parameter param_realdata")}

    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params[,3]=log(base=10,match_params[,3])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod,scale_s=scale_s)
  }else{
    params <- Get_params2(gene_effects,evf_res[[1]],bimod,ranges)
  }
  # kinetic params
  names(params) <- c("k_on", "k_off", "s")

  # update the kinetic value of the highly expressed genes
  hge <- includeHge(prop_hge=prop_hge, mean_hge=mean_hge, params)
  chosen_hge <- hge[["chosen_hge"]]
  params <- hge[["params"]]


  # degradation rate of the size ngenes * 1
  d_vec <- rnorm(n=ngenes, mean=0.8, sd=0.1)
  # splicing rate of the size ngenes * 1
  beta_vec <- rnorm(n=ngenes, mean=1.2, sd=0.1)

  params_cycle = list(k_on = params$k_on[,1:ncells_cycle], k_off = params$k_off[,1:ncells_cycle], s = params$s[,1:ncells_cycle])
  params_tree = list(k_on = params$k_on[,(ncells_cycle+1):(ncells_cycle+ncells_tree)], k_off = params$k_off[,(ncells_cycle+1):(ncells_cycle+ncells_tree)],
                     s = params$s[,(ncells_cycle+1):(ncells_cycle+ncells_tree)])

  # if the start state is not provided
  if(is.null(start_s) | is.null(start_u) | is.null(start_state)){
    root_state <- sample(c(1,2), size = ngenes, replace = TRUE)
    start_s <- rep(0, ngenes)
    start_u <- rep(0,ngenes)
    start_state <- root_state
  }

  start_cell_time <- 0
  # remove unstable cells
  randpoints <- runif(n=ncells_cycle, min=0, max=1)

  result <- gen_1branch(kinet_params = params_cycle, start_state = start_state, start_s = start_s, start_u = start_u,
                        randpoints1 = randpoints, ncells1 = ncells_cycle, ngenes = ngenes, beta_vec = beta_vec,
                        d_vec = d_vec)

  final_kinetics_cycle <- list(k_on = result$kinet_params$k_on, k_off = result$kinet_params$k_off,
                               s = result$kinet_params$s, beta = result$kinet_params$beta,
                               degrade = result$kinet_params$d)

  dynamics_cycle <- list(unspliced_steps = result$dynamics[["unspliced_steps"]], spliced_steps = result$dynamics[["spliced_steps"]], which_cells = result$dynamics[["which_cells"]])

  result_cycle <-list(counts_u = result$counts_u, counts_s=result$counts_u, kinet_params=final_kinetics_cycle, state_mat=result$state_mat,
                      cell_time=result$cell_time, velocity=result$velocity, dynamics = dynamics_cycle,
                      backbone = evf_res_cycle[[2]]$pop, chosen_hge = chosen_hge)

  # Tree part
  # number of edges that connect to the current node
  degrees <- table(c(phyla$edge[,1],phyla$edge[,2]))

  # root node is the node with no incoming edges
  root <- as.numeric(which(lapply(seq(1,length(degrees)), function(x){is.na(match(x, phyla$edge[,2]))}) == TRUE))

  # find starting and ending cells for each branches
  end_cell <- vector()
  branches <- vector()
  num_cells <- vector()

  for(pop in unique(evf_res_tree[[2]]$pop)){
    # find all cell indices correspond to the current population
    index <- which(evf_res_tree[[2]]$pop == pop)
    # the index is ordered according to depth originally
    end_cell <- c(end_cell, index[length(index)])
    branches <- rbind(branches, as.numeric(strsplit(pop, split ="_")[[1]]))
    num_cells <- c(num_cells,length(index))
  }
  # create a backbone_info matrix
  backbone <- cbind(start_node=branches[,1], end_node = branches[,2], end_cell = end_cell, start_cell = rep(0, length(end_cell)), num_cells = num_cells)

  for(i in 1:dim(backbone)[1]){
    if(backbone[i,"start_node"] == root){
      backbone[i,"start_cell"] <- 0
    }else{
      index <- which(backbone[,"end_node"] == backbone[i, "start_node"])
      backbone[i,"start_cell"] <- backbone[index,"end_cell"]
    }
  }
  # order backbone according to the start cells, in case the start cells are always generated
  backbone <- backbone[order(backbone[,"start_cell"]),,drop = FALSE]


  # generate counts
  # unspliced counts matrix
  counts_u <- vector()
  # spliced counts matrix
  counts_s <- vector()
  # cell developmental time
  cell_time <- vector()
  # count step changes
  steps_u <- vector()
  steps_s <- vector()
  which_cells <- vector()
  # promoter states
  state_mat <- vector()
  # RNA velocity
  velocity <- vector()

  # generate the promoter states of the root cell with dimension ngenes * 1
  root_state <- sample(c(1,2), size = ngenes, replace = TRUE)

  # loop through all branches
  for(i in 1:dim(backbone)[1]){
    kinetic_params <- list()
    # indices of the cells that belongs to current backbone
    index <- which(evf_res_tree[[2]]$pop == paste0(backbone[i, "start_node"], "_", backbone[i, "end_node"]))
    # subset k_on, k_off and s

    for(j in names(params)){
      # dimension ngenes by #index(subset of cells)
      kinetic_params[[j]] <- params_tree[[j]][,index]
    }
    # the points we use to get snapshot for each cell
    randpoints <- runif(n=backbone[i, "num_cells"], min=0, max=1)


    start_s = result_cycle$counts_s[,ncells_cycle]
    start_u = result_cycle$counts_u[,ncells_cycle]
    start_state = result_cycle$state_mat[,ncells_cycle]

    if(backbone[i, "start_node"] == root){
      if(is.null(start_s) | is.null(start_u) | is.null(start_state)){
        start_s <- rep(0, ngenes)
        start_u <- rep(0,ngenes)
        start_state <- root_state
      }
      start_cell_time <- 0
    }else{
      # start cell must have already been generated, and note that the index must be larger than 1
      start_s <- counts_s[,backbone[i,"start_cell"]]
      start_u <- counts_u[,backbone[i,"start_cell"]]
      start_state <- state_mat[,backbone[i,"start_cell"]]
      start_cell_time <- cell_time[backbone[i,"start_cell"]]
    }

    result <- gen_1branch(kinet_params = kinetic_params, start_state = start_state, start_u = start_u, start_s = start_s,
                          randpoints1 = randpoints, ncells1 = backbone[i, "num_cells"], ngenes = ngenes, beta_vec = beta_vec, d_vec = d_vec)

    counts_u <- cbind(counts_u, result$counts_u)
    counts_s <- cbind(counts_s, result$counts_s)
    state_mat <- cbind(state_mat, result$state_mat)

    steps_s <- cbind(steps_s, result$dynamics[["spliced_steps"]])
    steps_u <- cbind(steps_u, result$dynamics[["unspliced_steps"]])
    which_cells <- cbind(which_cells, result$dynamics[["which_cells"]])

    # true time of each cell, start from zero for each branch, need to update when merged together
    cell_time <- c(cell_time, result$cell_time + start_cell_time)
    velocity <- cbind(velocity, result$velocity)
  }

  final_kinetics_tree <- list(k_on = params_tree[[1]], k_off = params_tree[[2]],
                              s = params_tree[[3]], beta = beta_vec, degrade = d_vec)

  dynamics_tree <- list(unspliced_steps = steps_u, spliced_steps = steps_s, which_cells = which_cells)

  result_tree <- list(counts_u = counts_u, counts_s=counts_s, kinet_params=final_kinetics_tree, state_mat=state_mat,
                      cell_time=cell_time, velocity=velocity, dynamics = dynamics_tree,
                      backbone = evf_res_tree[[2]]$pop, chosen_hge = chosen_hge)


  # merge results
  final_kinetics <- list(k_on = cbind(result_cycle$kinet_params$k_on, result_tree$kinet_params$k_on),
                         k_off = cbind(result_cycle$kinet_params$k_off, result_tree$kinet_params$k_off),
                         s = cbind(result_cycle$kinet_params$s, result_tree$kinet_params$s),
                         beta = result_cycle$kinet_params$beta, degrade = result_cycle$kinet_params$d)

  final_dynamcis <- list(unspliced_steps = cbind(dynamics_cycle$unspliced_steps, dynamics_tree$unspliced_steps),
                         spliced_steps = cbind(dynamics_cycle$spliced_steps, dynamics_tree$spliced_steps),
                         which_cells = cbind(dynamics_cycle$which_cells, dynamics_tree$which_cells))

  result <- list(counts_u = cbind(result_cycle$counts_u, result_tree$counts_u),
                 counts_s=cbind(result_cycle$counts_s, result_tree$counts_s),
                 kinet_params=final_kinetics,
                 state_mat=cbind(result_cycle$state_mat, result_tree$state_mat),
                 cell_time=c(result_cycle$cell_time, result_cycle$cell_time[length(result_cycle$cell_time)] + result_tree$cell_time),
                 velocity=cbind(result_cycle$velocity, result_tree$velocity),
                 dynamics=final_dynamcis,
                 backbone=c(result_cycle$backbone, result_tree$backbone, chosen_hge = chosen_hge))

  return(result)
}

#' Add technical noise to the data
#' @param result result produced by simulation
#' @param capture.rate capture rate of RNA molecules, between 0 and 1
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @export
technicalNoise <- function(result, capture.rate = 0.2){
  counts_u <- round(result$counts_u)
  counts_s <- round(result$counts_s)

  counts_u2 <- apply(counts_u, 1:2, function(x){
    if(x!=0){
      x = rbinom(n=1, size=x, p=capture.rate)
    }
    return(x)
  })

  counts_s2 <- apply(counts_s, 1:2, function(x){
    if(x!=0){
      x = rbinom(n=1, size=x, p=capture.rate)
    }
    return(x)
  })
  result$counts_u = counts_u2
  result$counts_s = counts_s2
  return(result)
}
