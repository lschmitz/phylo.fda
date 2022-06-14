# Phylogenetic flexible disciminant analysis

# Lars Schmitz, 2014

# 1 Preliminaries
# 2 Analysis

################################################################################################

# 1 Preliminaries

  # Loading libraries

    library(ape)
    library(class)
    library(geiger)
    library(lattice)
    library(mda)
    library(nnet)
    source("phylo.fda.v0.2.R")

  # Reading in the tree.

    treA <- read.tree("example.tree.tre")
    if(!is.binary.tree(treA)) treA <- multi2di(treA, random = TRUE) # Randomly resolving polytomies if there are any

  # Reading in the data and assigning rownames.

    ddA <- read.csv("measurements.csv")
    rownames(ddA) <- ddA$taxon

  # Ordering data to match tip labels.
  # Note that for this tutorial data and phylogeny match.
  # Normally that's not the case. You can use 'geiger's treedata() function for that, as explained earlier.

    ddA <- ddA[match(treA$tip.label,rownames(ddA)),]


################################################################################################


# 2 Analysis

  # Defining groups and taxa.
    gA <- ddA$groups # contains data on ecology/behavior
    taxaA <- ddA$taxon # species names
    rownames(ddA) <- taxaA # assigning species names to rows       

  # Tree and data preparation: all taxa (test and training, e.g., fossil and living) or training only (living).
    XA <- ddA[,2:4] # we are selecting 3 variables from the dataset
    XA <- log10(XA) #log10 transformation
    XA <- signif(XA, 4)
    testtaxa <- rownames(ddA[gA=="unknown",]) # specifying taxa with unknown group, e.g. fossils
    testtaxan <- row(ddA)[gA=="unknown",1]
    trainingtaxa <- rownames(ddA[-testtaxan,]) # creating a dataframe that only contains taxa with known group affiliation
    X <- XA[-testtaxan,]
    dd <- ddA[-testtaxan,]
    g <- gA[-testtaxan]
    tre <- drop.tip(treA, testtaxa  )
    is.ultrametric(tre) # Tree should be ultrametric after removing fossils.
    
    #A <- ddA[ddA$groups=="unknown",]
    #A <- A[,2:4]
    #A <- log10(A)
    #A <- signif(A, 4)

  # Identifying optimal lambda: where is the strongest correlation between form and ecology among living taxa?

    filename_stem <- "ecology" # A plot will appear in separate window; a pdf is saved in your working directory.
    ol1 <- optLambda(X,g,tre,idc=filename_stem )
    ol1$optlambda # displaying the optimal lambda value

  # Performing the discriminant analysis

    pri  <- c(0.1427,0.5864,0.2709) # The prior probabilities will vary with your datasets. Replace by equal priors if unknown.
    optl <- 0.08 # Replace with the optimal lambda value from above.
    pfda <- phylo.fda.pred(XA,gA,taxaA,treA,testtaxan,val=optl,priin=pri)

  # Let's Compile cross-classification matrix how well pFDA performed.

    pfda$confusion # Misclassified poroportion listed as "error".

  # Now, let's extract predictions for the fossils.

    test.class <- as.character(predict(pfda, pfda$DATAtest, type="class"))

  # Discriminant scores scores for the fossils are retreived this way.

    test.variates <- predict(pfda, pfda$DATAtest, type="variates")

  # Posterior probabilities of group affiliations are available, too.

    test.prob <- predict(pfda, pfda$DATAtest, type="posterior")

  # And now let's put everything together (here shown for the fossils):

    test.results <- cbind(test.class, test.prob, test.variates)
    colnames(test.results) <- c("predicted class", "P(c)", "P(d)", "P(n)", "DA1", "DA2")
    rownames(test.results) <- testtaxa

    test.results # Let's take a look what we got.

  # The same can be done for the training data (living species), which is useful for evaluating how well pFDA performed.

    training.class <- as.character(predict(pfda, pfda$DATA, type="class"))
    training.variates <- predict(pfda, pfda$DATA, type="variates")
    training.prob <- predict(pfda, pfda$DATA, type="posterior")
    training.results <- cbind(as.character(g), training.class, training.prob, training.variates)
    colnames(training.results) <- c("true class", "predicted class", "P(c)", "P(d)", "P(n)", "DA1", "DA2")
    rownames(training.results) <- trainingtaxa

  # The DA scores (colored by groups) can be plotted using the following lines.
  # First some data organization of the lving species

    training <- as.data.frame(cbind(training.variates, as.character(dd$groups)))
    colnames(training) <- c("DA1", "DA2", "groups")
    rownames(training) <- dd$taxon

  # Followed by organization of fossils

    unknown <- as.character(rep("unknown", times=length(testtaxan)))
    test <- as.data.frame(cbind(test.variates, unknown))
    colnames(test) <- c("DA1", "DA2", "groups")
    rownames(test) <- testtaxa

  # All put together in one object, before broken into groups

    scatter <- rbind(training, test)
    scatter[, 1:2] <- lapply(scatter[,1:2], as.character)
    scatter[, 1:2] <- lapply(scatter[,1:2], as.numeric)

    c <- scatter[scatter$groups=="cathemeral",] 
    d <- scatter[scatter$groups=="diurnal",] 
    n <- scatter[scatter$groups=="nocturnal",]
    f <- scatter[scatter$groups=="unknown",] 

  # Finally we can plot the discriminant scores. 
  # Note that this is no longer a classical morphospace.
  # In addition to shape and size, phylogenetic covariance influences position of species.

    plot(c$DA1, c$DA2, 
         asp=1, 
         xlim=range(scatter[,1]),
         ylim=range(scatter[,2]), 
         pch=21, col="black", bg="gray", 
         cex=1.25, cex.lab=1.5, 
         xlab="DA 1", ylab="DA 2");
    points(d$DA1, d$DA2, pch=21, cex=1.25, col="black", bg="gold");
    points(n$DA1, n$DA2, pch=19, cex=1.25, col="black");
    points(f$DA1, f$DA2, pch=24, cex=1.25, col="black", bg="green");
    box(lwd=2); axis(1, lwd=2, lwd.ticks=2); axis(2, lwd=2, lwd.ticks=2)

# It is also highly recommended to explore how choice of optimal lambda influences the classification of fossils.
# Let me know if you have any questions.

# Key references:
# Hastie, T., Tibshirani, R., Buja, A. (1994) Flexible disriminant analysis by optimal scoring.
#   J Americ Stat Assoc, 89, 1255-1270.
# Motani, R., Schmitz, L. (2011) Phylogenetic versus functional signals in the evolution of form-function relationships in terrestrial vision. 
#   Evolution, 65, 2245-2257.
# Schmitz, L., Motani, R. (2011) Nocturnality in dinosaurs inferred from scleral ring and orbit morphology. 
#   Science, 332, 705-708.

################################################################################################
