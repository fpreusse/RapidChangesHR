####
### Tree BH procedure
####
## Input:
# pvals: Vector of N p-values (some of which may be NA) for most granular level of tree
# groups: N x L matrix with groups[n, l] = k means that hypothesis n belongs 
#         to the kth group in level l. Final column should be 1:N since most granular grouping is individual hypotheses
# q: Vector of L error rates to be targeted for each level in tree
# test: Vector of L-1 combination methods to obtain level l-1 p-values from level l  p-values. 
#       Alternatives are "simes", "fisher", "arbitrary". Default is "simes" at each level
# procedure:Vector of L procedures used to control the FDR. Possible: "BH" or "BY" or NULL (this allows for weak dependency at each level)
## Output:
#   sel:      N x L binary matrix where sel[n, l] = 1 for first element of 
#             level l group indicates that group was selected
## Remark:
# This function is an adaptation of the get_TreeBH_selection function by Christine B. Peterson, available at https://github.com/cbpeterson/TreeBH
# Exhaustive information can be found in Bogomolov et al (2021) Hypotheses on a tree: new error rates and testing strategies. https://doi.org/10.1093/biomet/asaa086

get_TreeBH_updated <- function(pvals, groups, q, test = "simes", procedure= NULL) {  
  # Total number of hypotheses at most granular level
  N <- length(pvals)
  
  # Total number of levels
  L <- ncol(groups)
  
  # Input checks
  if (length(q) != L) {
    stop("Must specify target q for each level of hierarchy")
  }
  if (min(q) < 0 || max(q) > 1) {
    stop("Target q must be in the range 0 - 1")
  }
  if (min(pvals, na.rm = TRUE) < 0 || max(pvals, na.rm = TRUE) > 1) {
    stop("P-values must be in the range 0 - 1")
  }
  if (length(pvals) != nrow(groups)) {
    stop("Dimension mismatch between pvals and groups")
  }
  if(!is.null(procedure)& length(procedure)==1){
    procedure <- rep(procedure, L)
  }
  if(is.null(procedure)){
    procedure <- rep("BH", L)
  }
  if(length(procedure)!= L){
    stop("If Procedure change between levels: procedure has to be defined for each level")
  }
  for (cur_level in 2:L) {
    cur_groups <- unique(sort(groups[ , cur_level]))
    for (cur_group in cur_groups) {
      cur_group_inds <- which(groups[, cur_level] == cur_group)
      parent_groups <- groups[cur_group_inds, cur_level - 1]
      if (length(unique(parent_groups)) != 1) {
        stop("Groups must be nested within hierarchy")
      }
    }
  }
  
  if (!length(unique(groups[ , L])) == nrow(groups)) {
    stop("Assumption is that lowest level in tree corresponds to individual hypotheses")
  }

  sel <- matrix(0, nrow = N, ncol = L)

  if (length(test) == 1) {
    test <- rep(test, L - 1)
  }
  q_adj <- matrix(NA, nrow = N, ncol = L)
  q_adj[ , 1] <- q[1]
  group_pvals <- matrix(NA, nrow = N, ncol = L)
  group_pvals[ , L] <- pvals
  aggregate_levels <- (L-1):1
  for (cur_level in aggregate_levels) {
    cur_groups <- unique(sort(groups[ , cur_level]))
    for (cur_group in cur_groups) {
      cur_group_inds <- which(groups[, cur_level] == cur_group)
      cur_pvals <- group_pvals[cur_group_inds, cur_level + 1]
      if (test[cur_level] == "simes") {
        group_pvals[cur_group_inds[1], cur_level] <- get_simes_p(cur_pvals)
      } else if (test[cur_level] == "fisher") {
        group_pvals[cur_group_inds[1], cur_level] <- get_fisher_p(cur_pvals)
      } else if(test[cur_level] =="arbitrary"){
        group_pvals[cur_group_inds[1], cur_level] <- get_heller_p(cur_pvals)
      } else {
        stop("Options for parameter test are 'simes', 'fisher' and 'arbitrary'")
      }
    }
  }
  for (cur_level in 1:L) {
    # Stop if there were no selections in previous level
    if (cur_level > 1 && sum(sel[ , cur_level - 1]) == 0) {
      break
    }
    if (cur_level == 1) {
      sel_prev <- 1 
    } else {
      sel_prev <- which(sel[ , cur_level - 1] == 1)
    }
    for (parent_sel in sel_prev) {
      # Identify selected children of current parent
      if (cur_level == 1) {
        parent_group_ind <- 1
        child_inds <- which(!is.na(group_pvals[ , cur_level]))
      } else {
        parent_group_num <- groups[parent_sel, cur_level - 1]
        parent_group_ind <- min(which(groups[ , cur_level - 1] == parent_group_num))
        child_inds <- which(groups[ , cur_level - 1] == parent_group_num &
                              !is.na(group_pvals[ , cur_level]))
      }
      if (length(child_inds) > 1) {
        if(is.null(procedure)){
          sel_ind_within_group <- which(qvalue(group_pvals[child_inds, cur_level],
                                               lambda = 0)$qvalue <=
                                          q_adj[parent_group_ind, cur_level])
        }else{
          if(procedure[cur_level]== "BH"){
            fam_size <- length(group_pvals[child_inds, cur_level])
            t_vals <- q_adj[parent_group_ind, cur_level]*(1:fam_size)/fam_size
            rej_in_group <- sort(group_pvals[child_inds, cur_level])[sort(group_pvals[child_inds, cur_level])<= t_vals]
            sel_ind_within_group <- which(group_pvals[child_inds, cur_level] %in% rej_in_group)
          }
          if(procedure[cur_level]=="BY"){
            fam_size <- length(group_pvals[child_inds, cur_level])
            t_vals <- (1:fam_size)/fam_size*q_adj[parent_group_ind, cur_level]/sum(1/(1:fam_size))
            rej_in_group <- sort(group_pvals[child_inds, cur_level])[sort(group_pvals[child_inds, cur_level])<= t_vals]
            sel_ind_within_group <- which(group_pvals[child_inds, cur_level] %in% rej_in_group)
          }
        }
      } else {
        sel_ind_within_group <- which(group_pvals[child_inds, cur_level] <=
                                        q_adj[parent_group_ind, cur_level])
      }
      sel[child_inds[sel_ind_within_group], cur_level] <- 1
      if (cur_level < L) {
        R_parent_sel <- length(sel_ind_within_group)
        n_parent_sel <- length(child_inds)
        
        q_adj[child_inds[sel_ind_within_group], cur_level + 1] <-
          q[cur_level + 1] *
          q_adj[parent_sel, cur_level] / q[cur_level] *
          R_parent_sel / n_parent_sel
      }
    }
  }
  sel
}

### Functions to combine p-values:
# using Simes approach:
get_simes_p <- function(pvals) {
  # Don't include NAs
  pvals <- pvals[!is.na(pvals)]
  
  if (length(pvals) == 0) {
    # Return NA if input had length 0 or was all NAs
    NA
  } else {
    # Else compute Simes p
    pvals <- sort(as.numeric(pvals))
    length(pvals) * min(pvals / seq(1:length(pvals)))
  }
}

# using Fisher approach:
get_fisher_p <- function(pvals) {
  # Don't include NAs
  pvals <- pvals[!is.na(pvals)]
  
  if (length(pvals) == 0) {
    # Return NA if input had length 0 or was all NAs
    NA
  } else {
    # Else compute Fisher p-value
    chisq_df <- 2 * length(pvals)
    pchisq(-2 * sum(log(pvals)), chisq_df, lower.tail = FALSE)
  }
}

# For arbitrary dependency (by Heller et al, also given in Dickhaus, 2014)
get_heller_p <- function(pvals) {
  # Don't include NAs
  pvals <- pvals[!is.na(pvals)]
  
  if (length(pvals) == 0) {
    # Return NA if input had length 0 or was all NAs
    NA
  } else {
    # Else compute p-value under arbitrary dependence
    pvals <- sort(as.numeric(pvals))
    length(pvals) * pvals[1]
  }
}

####
### Inheritance Procedure
####
## Input:
# family_size: vector which includes the number of hypotheses that belong to one family at each level of the hierarchical structure
#              That is, if 5 elementary hypotheses belong to level L-1, at the elementary level we have 5.
#               Note: besides the first and last entry, the family_size should be equal to 2.
# test:        Vector of L-1 combination methods to obtain level l-1 p-values from
#              level l  p-values. Alternatives are "simes", "fisher", "arbitrary". Default is
#              "simes" at each level (i.e., approach by Heller et al for PRDS)
#               If the dependence is the same at all levels: ok to just use one entry in vector
# alpha:       (targeted) level of FWER control
# p_vals:       Vector of the p-values of the elementary hypotheses. Need to be ordered as they would be in a tree (i.e., family-wise grouping)
## Notes:
# In this application, all groups have the same size. Therefore, 
# we can use how many hypotheses are inlcuded per level and start "counting" from the beginning.
# Entries in family_size and test correspond to the level of the hierarchy, i.e., last entry is the elementary level.
# This code is only valid if the hypotheses at the levels between the first (excluding the global null) and the last only have two siblings.

inheritance_selection <- function(family_size, test=c("arbitrary"), alpha, p_vals){
  # check if family_size and elementary_p work together
  # for this, p_vals needs to be a multiple of the size of the families at each level
  n_levels <- length(family_size)
  n_elementary <- length(p_vals)
  for(i in 1:n_levels){
    if(n_elementary%%family_size[i]!=0){
      stop("number of elementary hypotheses and size of families does not match.")
    }
  }
  if(length(test)==1){
    test <- rep(test, n_levels)
  }
  # Compute the p-values and critical values for each hypothesis in the tree
  # we start "bottom up" in a recursive manner
  # last one is the global null
 
  p_vals_tree <-critical_value <- vector("list",n_levels+1)
  # this is a list which contains the p-values at each level. The last entry of the list are the elementary hypotheses
  
  # in our case: each branch in the tree has the same number of leaf nodes (i.e., number of elementary hypotheses)
  cv <- alpha/cumprod(c(1,family_size))
  
  n_families <- length(p_vals)/family_size[n_levels]
  p_vals_tree[[n_levels+1]] <- data.frame(p_values=p_vals, 
                                          id_parent= rep(c(1:n_families), each=family_size[n_levels]))
  critical_value[[n_levels+1]] <- rep(cv[n_levels+1], length(p_vals))
  for(i in (n_levels):1){
    # how many hypotheses are there at the current level of the tree?
    n_families <- nrow(p_vals_tree[[i+1]])/family_size[i]
    if(n_families%%1!=0)stop("Size of families dos not match")
    # initialize the critical values for each hypothesis in the tree
    critical_value[[i]] <- rep(cv[i], n_families)
    # compute the p-values for each hypothesis in the tree
    p_families <- numeric(n_families)
    if(test[i]=="simes"){
      for(j in 1:n_families){
        p_families[j] <- get_simes_p(pvals= p_vals_tree[[i+1]][((j-1)*family_size[i]+1):(j*family_size[i]),1])
      }
    }else if(test[i]=="fisher"){
     for(j in 1:n_families){
       p_families[j] <- get_fisher_p(pvals= p_vals_tree[[i+1]][((j-1)*family_size[i]+1):(j*family_size[i]),1])
     }
    }else{
      for(j in 1:n_families){
        p_families[j] <- get_heller_p(pvals= p_vals_tree[[i+1]][((j-1)*family_size[i]+1):(j*family_size[i]),1])
      }
    }
    if(i>=2){
      n_parent <- length(p_families)/family_size[i-1]
      p_vals_tree[[i]] <- data.frame(p_values= p_families,
                                     id_parent= rep(c(1:n_parent), each=family_size[i-1]))
    }else{
      p_vals_tree[[i]] <- data.frame(p_values= p_families,
                                     id_parent= rep(0, length(p_families)))
    }
  }
  remove(cv)
  # Check in top-down approach which hypotheses are rejected. The critical values are the same for all hypotheses at one level
  # is global hypothesis rejected?
  # elementary hypotheses are not included here, is done seperatly to account for Shaffer's factor:
  if(p_vals_tree[[1]][,1]<=critical_value[[1]][1]){
    rejected_level <- rejection_new(critical_values= critical_value, p_vals_tree=p_vals_tree, considered_hypotheses=NULL)
    if(rejected_level$new_rejections){
      rejected_level <- rejected_level$rejected_level
    }else{
      return(list(p_values = p_vals_tree[[1]],
                  critical_value = critical_value[[1]][1],
                  rejected = c(1)))
    }
  }else{
    # in this case, not even the global hypothesis is rejected, so we can terminate the function. Return 0 to indicate that nothing has been rejected.
    return(list(p_values = p_vals_tree[[1]],
                critical_value = critical_value[[1]][1],
                rejected = c(0)))
  }
  rejected_per_iteration <- critical_vals <- vector("list",0)
  update <- T
  rej_elementary <- numeric(0)
  while(update){
    # only consider those elementary hypotheses whose parent hypotheses have been rejected in the new round
    cons_fam_element <- rejected_level[[n_levels]] 
    within_family_heritance <- heritance_func(prev_rejected=rej_elementary, 
                                              p_values_tree= p_vals_tree[[n_levels+1]]$p_values, 
                                              family_size= family_size, 
                                              n_levels= n_levels, 
                                              critical_values= critical_value[[n_levels]],
                                              considered_families=cons_fam_element)
    rejected_level[[n_levels+1]] <- within_family_heritance$elementary_rejections
    # update the rejected elementary hypotheses: 
    rej_elementary <- unique(c(rej_elementary, rejected_level[[n_levels+1]]))
    # save the rejections
    rejected_per_iteration <- c(rejected_per_iteration, list(rejected_level))
    # save the utilized critical values:
    critical_value[[n_levels+1]]<- within_family_heritance$updated_critical_values
    critical_vals <- c(critical_vals, list(critical_value))
    # If there is no wealth to be inherited by others in a next step: break the while loop
    if(sum(within_family_heritance$free_wealth, na.rm=T)==0){
      break}
    
    # otherwise, update the critical values using the update_critical_vals function
    # check if anything changes here using update_needed
    # careful: rejected_level is relevant only for the considered hypotheses relevant, not for every single one.
    cvs <- update_critical_vals(prev_rejected=rej_elementary, 
                                p_values_tree= p_vals_tree, 
                                family_size= family_size, 
                                n_levels= n_levels, 
                                critical_vals= critical_value)
    critical_value <- cvs$updated_cv
    # if no update is needed because there are no changes:
    if(!cvs$update_needed)break
    # else: update the rejected_per_iteration list using rejection_new
    # check if anything has changed using new_rejections (if this is false, no new hypothesis has been rejected)
    rejection_updated <- rejection_new(critical_values= critical_value, 
                                       p_vals_tree= p_vals_tree, 
                                       considered_hypotheses=cvs$new_hypotheses)
    rejected_level <- rejection_updated$rejected_level
    update <- rejection_updated$new_rejections
  }
  
  # now summarize all rejections in one single list:
  n_iterations <- length(rejected_per_iteration)
  rejection_full <- vector("list", length=n_levels) # we utilize n_levels instead of n_levels +1 because it is not necessary to indicate that global hypothesis has been rejected
  for(i in 1:n_levels){
    rej_level <- numeric(0)
    for(j in 1:n_iterations){
      rej_level <- unique(c(rej_level,rejected_per_iteration[[j]][[i+1]]))
    }
    rejection_full[[i]] <- sort(rej_level)
  }
  
  return(list(rejections= rejection_full,
              rejections_per_iteration=rejected_per_iteration,
              critival_values_per_iteration=critical_vals ))
  
}

#### Helper functions needed for inheritance_selection:
### Function for inheritance at the elementary level (within a family):
## Input
# prev_rejected: vector indicating which elementary hypotheses are already rejected
# p_values_tree: p_values of the elementary hypotheses
# family_size: vector indicating the size of the families per level
# n_levels: number of levels in the hierarchy, corresponds to the length of the family_size
# critical_values: vector of wealth per family at the elementary level. If it is the same for all: only use one value
# considered_families: vector indicating which families are considered, should coincide with considered parent hypotheses. Corresponds to the (newly) rejected parents in each round
## Output
# vector indicating which elementary hypotheses are rejected after passing on wealth within the family
# logical vector indicating which family is rejected fully and can pass on wealth to other family
# save_critical_vals: vector containing the critical values used to reject elementary hypotheses
heritance_func <- function(prev_rejected, p_values_tree, family_size, n_levels, critical_values, considered_families=NULL){
  # first, check at the elementary level which hypotheses have been already rejected and if any wealth leaves the family:
  if(length(critical_values)==1){
    n_families <- length(p_values_tree)/family_size[n_levels]
    critical_values <- rep(critical_values, n_families)
  }else{
    n_families <- length(critical_values)
  }
  updated_rejections <- numeric(0)
  pass_on <- logical(n_families)
  if(is.null(considered_families)){
    fam_cons <- c(1:n_families)
  }else{
    fam_cons <- considered_families
  }
  save_critical_vals <- numeric(0)
  ind <- 1
  for(i in fam_cons){
    id_fam <- c(((i-1)*family_size[n_levels]+1):(i*family_size[n_levels]))
    # check how many elementary hypotheses have been rejected in this family. 
    n_rej_start <- n_rej <- sum(id_fam%in% prev_rejected)
    
    # if no elementary hypothesis has been rejected in a family: update the critical values using Shaffe's factor:
    if(n_rej==0){
      new_critical_value <- critical_values[ind]/(family_size[n_levels]-1)
      update_rej <- id_fam[which(p_values_tree[id_fam]<=new_critical_value)]
      n_rej <- length(update_rej)
      save_critical_vals <- c(save_critical_vals, new_critical_value)
    }
    
    # If more than one but not all elementary hypotheses in a family have been rejected: update the critical values at the elementary_level
    if(n_rej >=1 & n_rej < family_size[n_levels]){
      new_round <- T
    }else{new_round <- F}
    # now start within the family at the elementary level to make new contributions
    while(new_round){
      new_size <- family_size[n_levels]-n_rej
      # update the critical value, i.e. distribute the wealth of the parent anew
      # we test all hypotheses because the previously rejected hypotheses will also be included with the updated critical value.
      new_critical_value <- critical_values[ind]/new_size
      save_critical_vals <- c(save_critical_vals, new_critical_value)
      update_rej <- id_fam[which(p_values_tree[id_fam]<=new_critical_value)]
      n_rej_new <- length(update_rej)
      if(n_rej_new >=2 & n_rej_new < family_size[n_levels]){
        new_round <- T
      }else{new_round <- F}
      if(n_rej==n_rej_new){new_round <- F}
      n_rej <- n_rej_new
    }
    # 1) update the rejections
    if(n_rej >n_rej_start){
      updated_rejections <- c(updated_rejections, update_rej)
    }
    # 2) indicate if wealth is "free" to be passed on for each considered family
    if(n_rej == family_size[n_levels]){
      pass_on[ind] <- T
    }else{pass_on[ind] <- F}
    ind <- ind+1
  }
  return(list(elementary_rejections= updated_rejections,
              free_wealth= pass_on,
              updated_critical_values= save_critical_vals))
}

### Function for inheritance across families
## Input
# prev_rejected: vector including all rejected elementary hypotheses
# p_values_tree: the p-values of all hypotheses in the hierarchy, as well as the corresponding parent hypotheses
# critical_values: vector of critical values per level (or per hypothesis?)
## Output:
# updated_cv: updated critical values, corresponds to the critical values that have been used to test the hypotheses in the second round
# new_hypotheses: which hypotheses have inherited wealth from their siblings?
# update_needed: True or False, states if anything can be changed
update_critical_vals <- function(prev_rejected, p_values_tree, family_size, n_levels, critical_vals){
  # dead_families indicate at the level n_levels-1 which families are rejected and can pass on their wealth to other families
  # check: 
  # 1) which parent hypotheses do the rejected families have (so what are potential siblings)
  # initialize a list which remembers at each level in the hierarchy which hypothesis have surviving descendants (on the elementary level: which hypotheses are not rejected)
  alive_families <- vector("list", n_levels+1)
  alive_families[[n_levels+1]] <- c(1:nrow(p_values_tree[[n_levels+1]]))
  # if rejected/ dead family: use zero in dead_families
  alive_families[[n_levels+1]][prev_rejected] <- 0
  
  # This gives the number of "surviving" children per parent, so how many hypotheses can inherit the wealth IF the parent is rejected
  # if no survive: use zero
  for(i in n_levels:1){
    # to get "alive" families at level i, we check the number of dead children in alive_families[[i+1]]
    n_families <- nrow(p_values_tree[[i+1]])/family_size[i]
    alive <- numeric(n_families)
    for(j in 1:n_families){
      id_fam <- c(((j-1)*family_size[i]+1):(j*family_size[i]))
      alive[j] <- sum(alive_families[[i+1]][id_fam]!=0)
    }
    alive_families[[i]] <- alive
  }
  # at each level of the tree: know for each hypothesis the number of children that it can inherit to
  # If a hypothesis has zero heirs, it itself as value 0 and will not be considered when updating the rejections
  # If the critical vector of an hypothesis does not change, it will not be included in further consideration
  # also: initialize list which includes indices of the hypotheses that have an updated critical value and need to be tested again:
  consider_next <- vector("list", n_levels+1)
  consider_next[[1]] <- 1
  new_critical_values <- critical_vals
  update_needed <- F
  for(i in 2:(n_levels+1)){
    # check: are there any hypotheses with surviving heirs at the current level of the hierarchy? 
    consider_next[[i-1]] <-id_new <- which(alive_families[[i-1]]>0)
    # this is the number of hypotheses with surviving heirs at the current level of the hierarchy
    # Split the corresponding wealth of the parent hypothesis (critical_value[[i-1]]) by the number of "alive" children (id_fam) to get the critical value at this level of the hierarchy
    # the number of alive children is given in alive_families
    if(length(id_new)!=0){
      update_needed <- T
      new_cv_helper <- numeric(0)
      for(j in 1:length(id_new)){
        #1) Update the critical values, so give wealth of parent only to the surviving siblings
        # We only have one sibling per level in our simulation study
        interim <- rep(new_critical_values[[i-1]][j]/alive_families[[i-1]][id_new[j]], alive_families[[i-1]][id_new[j]])
        new_cv_helper <- c(new_cv_helper, interim)
      }
      new_critical_values[[i]] <- new_cv_helper
    }
  }
  consider_next[[(n_levels+1)]] <- which(alive_families[[n_levels+1]]!=0)
  return(list(updated_cv= new_critical_values,
         new_hypotheses= consider_next,
         update_needed = update_needed))
}

### Function to test new hypotheses if anything is rejected:
## Input:
# critical_values: needs to be a list o length n_levels+1, including the critical values per hypothesis at each level of the tree
# p_vals_tree: the p_values for each hypothesis in the tree
# considered_hypotheses: list of length n_levels+1, indicating which hypotheses we are testing. If NULL: all hypotheses are considered
## Output
# rejected_level: the indices of the rejected hypotheses at this level.
# new_rejections: true or false if any new hypothesis has been rejected and could pass on wealth
rejection_new <- function(critical_values, p_vals_tree, considered_hypotheses=NULL){
  n_levels <- length(p_vals_tree)
  rejected_level <- vector("list", n_levels)
  rejected_level[[1]] <- 1 #the global hypothesis is always rejected, otherwise we would not even run this function
  new_rejections <- F
  
  if(is.null(considered_hypotheses)){
    for(i in 2:(n_levels-1)){
      # which hypotheses have a rejected parent?
      id_test <- which(p_vals_tree[[i]][,2]%in%rejected_level[[i-1]])
      # But also: I need the index of the hypothesis that is rejected (and not the index of the hypotheses with rejected parent)
      rejected_level[[i]] <- id_test[which(p_vals_tree[[i]][id_test,1]<=critical_values[[i]][1])]
      if(length(rejected_level[[i]]!=0))new_rejections <- T
    }
  }else{
    for(i in 2:(n_levels-1)){
      # which hypotheses are considered at each level:
      id_cons <- considered_hypotheses[[i]]
      # only test those hypotheses that have been considered AND have a rejected parent
      id_test <- id_cons[which(p_vals_tree[[i]][id_cons,2]%in%rejected_level[[i-1]])]
      # Return the index of the hypothesis that is rejected (and not the index of the hypotheses with rejected parent)
      rejected_level[[i]] <- id_test[which(p_vals_tree[[i]][id_test,1]<=critical_values[[i]][1])]
      if(length(rejected_level[[i]]!=0))new_rejections <- T
    }
  }
  return(list(rejected_level= rejected_level,
         new_rejections= new_rejections))
}

