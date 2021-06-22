# //
# //  Function_singleton_calculation_Mezza.r
# //  
# //
# //  Created by Massimo Mezzavilla on 24/06/2020.
# //

#ifndef Function_singleton_calculation_Mezza_h
#define Function_singleton_calculation_Mezza_h


#endif /* Function_singleton_calculation_Mezza_h */



ZTRF<-function (formula, data, family = gaussian){
          # formula=X~Y
          # family = gaussian
          if (is.character(family))
              family <- get(family, mode = "function", envir = parent.frame())
          if (is.function(family))
              family <- family()
          if (is.null(family$family)) {
              print(family)
              stop("'family' not recognized")
          }
          if (is(try(formula, silent = TRUE), "try-error")) {
              formula <- data[[as(match.call()[["formula"]], "character")]]
          }
          if (is(formula, "formula")) {
              mf <- model.frame(formula, data, na.action = na.pass,
                  drop.unused.levels = TRUE)
              mids <- complete.cases(mf)
              mf <- mf[mids, ]
              y <- model.response(mf)
              desmat <- model.matrix(formula, mf)
              lmf <- glm.fit(desmat, y, family = family)
              resid <- lmf$resid
          }
          else if (is(formula, "numeric") || is(formula, "integer") ||
              is(formula, "double")) {
              y <- formula
              mids <- (!is.na(y))
              y <- y[mids]
              resid <- y
              if (length(unique(resid)) == 1)
                  stop("trait is monomorphic")
              if (length(unique(resid)) == 2)
                  stop("trait is binary")
          }
          else {
              stop("formula argument must be a formula or one of (numeric, integer, double)")
          }
          y <- (resid - mean(resid))/sd(resid)
          tmeas <- as.logical(mids)
          out <- rep(NA, length(mids))
          out[tmeas] <- y
          out
      }


NTRF<-function (formula, data, family = gaussian){
    var <- ZTRF(formula, data, family)
    out <- rank(var) - 0.5
    out[is.na(var)] <- NA
    mP <- 0.5/max(out, na.rm = T)
    out <- out/(max(out, na.rm = T) + 0.5)
    out <- qnorm(out)
    out
}

##estimate DSC and SSC score
#Set arguments to use from command line
# 
#read commandline args
args <- commandArgs(trailing=TRUE)

#set population files names and paths
singleton_table <- args[[1]]
out_file <- args[[2]]

# singleton_table <- "mysingleton_data.txt"

a=read.table(singleton_table, h=T)

#the header file is preformatted as:
# GENE_ID GENE_NAME       SING_GENE       GENE_LENGTH     GENE_SING_DENSITY       TRANSCRIPT_ID   SING_TRANSCRIPT TRANSCRIPT_LENGTH       CDS_SING_DENSITY
X=(a$SING_GENE-a$SING_TRANSCRIPT)
Y=(a$GENE_LENGTH-a$TRANSCRIPT_LENGTH)
SC_ncds=NTRF(X~Y)
SC_cds=NTRF(a$SING_TRANSCRIPT~ a$TRANSCRIPT_LENGTH)
Atot3=data.frame(a[,1:9], SC_cds, SC_ncds)
DSC_score=rowDiffs(as.matrix(data.frame(Atot3$SC_ncds, Atot3$SC_cds)))
SSC_score=rowSums(data.frame(Atot3$SC_cds, Atot3$SC_ncds))
Atot4=data.frame(Atot3,DSC_score,SSC_score)

#write output file
out_folder <- dirname(out_file)
dir.create(out_folder, recursive=T)
write.table(Atot4,file=paste(out_file,sep=""),sep="\t",col.names=T,quote=F,row.names=F)