############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:
####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.
####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.
####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.
#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

#################################################################################################
###### This is for a demonstration of how SignatureAnalzyer works in 35 Billiary PCAWG WGS
#################################################################################################

CURRENT <- paste(getwd(),"/",sep="")
INPUT <- paste(CURRENT,"INPUT_SignatureAnalyzer/",sep="") ### Directory for INPUT data
OUTPUT <- paste(CURRENT,"OUTPUT_DEMO/",sep="") #### Directory for OUTPUT data
TEMPORARY <- paste(CURRENT,"TEMPORARY_SignatureAnalzyer/",sep="")
system(paste("mkdir",OUTPUT,sep=" "))
system(paste("mkdir",TEMPORARY,sep=" "))

source("SignatureAnalyzer.PCAWG.function.R")

###### loading input mutation counts matrix along 1536 SNV + 78 DNP + 83 INDEL features across 2780 PCAWG samples 
###### SNV features in 96 contexts were defined as four-letters such as "CTAA", in which first two letters refers to the reference and altered bases, and
###### next two letters represent a base at -1 and +1 position at mutated pyrimidines (C or T), e.g, C>T at ACA tri-nucleotide sequence motif.
###### SNV features in 1536 contexts were defined as six-letters such as "CTAATT", in which first two letters refers to the reference and altered bases, and
###### next four letters represent a base at -2, -1, +1, and +2 position at mutated pyrimidines (C or T), e.g, C>T at AACTT penta-nucleotide sequence motif.
###### SNV features in 96 contexts were defined as four-letters such as "CTAT", in which first two letters refers to the reference and altered bases, and
###### next two letters represent a base at -1 and +1 position at mutated pyrimidines (C or T), e.g, C>T at ACT tri-nucleotide sequence motif.

load(file=paste(INPUT,"lego1536.PAN.SNV.091217.RData",sep="")) ### SNV lego-matrix in 1536 contexts
load(file=paste(INPUT,"lego96.PAN.SNV.091217.RData",sep="")) ### SNV lego-matrix in 96 contexts
load(file=paste(INPUT,"lego1536.DNP.091217.RData",sep="")) ### DNP lego-matrix
load(file=paste(INPUT,"lego1536.INDEL.091217.RData",sep="")) ### INDEL lego-matrix

lego1536.PAN <- rbind(lego1536.SNV,lego1536.DNP,lego1536.INDEL) ### COMPOSITE lego-matrix in 1697 features (1536 SNV + 78 DNP + 83 INDEL)

###### loading information on hyper- or ultra-mutated samples; POLE (n=9), MSI (n=39), all SKIN (n=107), and single TMZ (CNS_GBM__SP25494)
load(file=paste(INPUT,"sample.POLE.RData",sep=""))
load(file=paste(INPUT,"sample.MSI1.RData",sep=""))
load(file=paste(INPUT,"sample.remove.RData",sep=""))
sample.TMZ <- c("CNS_GBM__SP25494")

###### We subset 35 Biliary_AdenoCA samples from the "lego96.SNV" (96 trinucleotide sequence contexts) and perform a signature extraction.
ttype <- sapply(colnames(lego96.SNV),function(x) strsplit(x,"__")[[1]][1])
lego96.demo <- lego96.SNV[,ttype=="Biliary_AdenoCA"]

Kcol <- 25         ### maximum number of signatures
tol <- 1.e-07	   ### tolerance for convergence; we lower this threshold to 1.e-05 for the PCAWG analysis.
method <- "L1W.L2H.LEGO96.Biliary"

###### signature extraction
###### about 26 minutes for 10 indepenent BayesNMF runs on the processon 2.6 GHz Intel Cori7 in MacBook (OS X El Captian)
for (i in 1:10) {
        res <- BayesNMF.L1W.L2H(as.matrix(lego96.demo),200000,10,5,tol,Kcol,Kcol,1)
	W <- res[[1]]
	K <- sum(colSums(W)>1)
        save(res,file=paste(OUTPUT,paste(method,i,"RData",sep="."),sep=""))
}

###### plotting signatures
n.run <- 10
summary.run <- array(0,dim=c(n.run,3))
for (i in 1:n.run) {
	load(file=paste(OUTPUT,paste(method,i,"RData",sep="."),sep=""))  ### loading RData file
	W <- res[[1]]  ### signature loading
	H <- res[[2]]  ### activity loading
	index.W <- colSums(W)>1 ### only keep columns in W with non-zero contributions
	W <- W[,index.W]
	H <- H[index.W,]
	colsum <- colSums(W)
	rowsum <- rowSums(H)
	### By scaling the signature loading matrix has all mutation burdens - each signture (column in W) now represents a number of mutations
	### assigned to each signature.
	for (j in 1:ncol(W)) {
		W[,j] <- W[,j]*rowsum[j]
		H[j,] <- H[j,]*colsum[j]
	}
	K <- ncol(W) ### number of extracted signatures
	colnames(W) <- paste("W",seq(1:K),sep="")
	p <- plot.signature.SNV(as.matrix(W),paste("DEMO",i,sep="."))  ### plotting signatures; res[[4]] = -log(posterior)
	pdf(file=paste(OUTPUT,paste("signature",method,K,round(res[[4]]),i,"pdf",sep="."),sep=""),width=(12),height=0.75*K)
		plot(p)
	dev.off()
	summary.run[i,1] <- i
	summary.run[i,2] <- K
	summary.run[i,3] <- -res[[4]]
}
colnames(summary.run) <- c("Run","K","posterior")
summary.run <- data.frame(summary.run)
summary.run <- summary.run[order(summary.run$posterior,decreasing=T),] ### summary.run is ordered by a posterior

###### from "summary.run" we note that seven runs converged to K=12 and three runs converged to K=13. 
###### we have chosen the solution (Run 7) having a maximum posterior at K=12; "L1W.L2H.LEGO96.Biliary.7.RData"
