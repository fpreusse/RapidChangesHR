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