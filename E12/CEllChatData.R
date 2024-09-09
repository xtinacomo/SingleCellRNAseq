
setwd("/Users/christinacomo/Desktop/E12.5 SOR data")
getwd()

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(uwot)
library(patchwork)
library(CellChat)

.error_if_no_Seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat installation required for working with Seurat objects")
  }
}

fibroblasts.E12.5 <- readRDS("E12.5Fibro.rds")
DimPlot(fibroblasts.E12.5)

# You can directly pass a Seurat Object to createCellChat:
cellchat <- createCellChat(object = fibroblasts.E12.5, group.by = "ident")

# Set the ligand-receptor interaction database

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)


# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

saveRDS(cellchat, file = "E12.5_fibroblasts_CellChat")
# Extract the inferred cellular communication network as a data frame  (i.e. easily access the inferred cell-cell communications of interest)

# 3 different examples of how you can use this:

# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set `slot.name = "netP"` to access the the inferred communications at the level of signaling pathways
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = "CellChat_R03",
          col.names = TRUE, row.names = TRUE, append = FALSE)


#LOOK AT FIGURES FROM ORIGINAL PAPER
# gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net2 <- subsetCommunication(cellchat, signaling = c("WNT"))


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#significantly different pathways
#"WNT"    "ncWNT"  "CXCL"   "MIF"    "ANGPTL" "MK"     "PTN"
cellchat@netP$pathways

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
### Bubble plot
#{r, fig.width=4,fig.height = 3, fig.wide = TRUE, fig.align = "center"}
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("WNT"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("WNT", "ncWNT"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
```


