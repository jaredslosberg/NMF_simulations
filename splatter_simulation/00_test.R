if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("Oshlack/splatter", dependencies = TRUE,
#                      build_vignettes = TRUE)

require(splatter)
require(scater)




params <- newSplatParams()
getParam(params, "nGenes")

set.seed(155)
sim_single <- splatSimulate(params, method = "single")

sim_single <- scater::logNormCounts(sim_single)
sim_single <- runPCA(sim_single)
plotPCA(sim_single)

sim_path <- splatSimulate(params, method = "path")
sim_path <- scater::logNormCounts(sim_path)
sim_path <- runPCA(sim_path)
plotPCA(sim_path, colour_by = "Step")

#Documentation for more complex (path) simulation
#http://oshlacklab.com/splatter/articles/splat_params.html#path-parameters-1

params_branching <- newSplatParams(batchCells = 1000, nGenes = 1000)
sim_path_branching <- splatSimulate(params_branching,
                                    method = "path",
                                    group.prob = rep(0.25,4),
                                    de.prob = 0.5,
                                    path.from = c(0,1,2,2))
sim_path_branching <- scater::logNormCounts(sim_path_branching)
sim_path_branching <- runPCA(sim_path_branching)
plotReducedDim(sim_path_branching, dimred = "PCA", colour_by = "Step")

sim_path_branching <- runUMAP(sim_path_branching, ncomponents = 50)
plotReducedDim(sim_path_branching, dimred = "UMAP", colour_by = "Group")


#1,2,6,8,23
gw <- cbind(gene_weights[,"NMF_pattern14"],gene_weights[,"NMF_pattern15"])
