suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
src <- c("Integration_Manuscript/Figure5_SharedNetworkArchitecture/Figure5_MasterData.tsv",
         "Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv")
paths <- c(file.path(root, src), "../01_baseline/results/sample_manifest.tsv")
m <- fread(paths[1], na.strings="NA")
a <- fread(paths[2], na.strings="NA")
b <- fread(paths[3], na.strings=c("NA",""))[sample_type=="Primary Tumor"]
x <- merge(m,a,by.x="PatientID",by.y="patient_id")
x <- merge(x,b[,.(patient_barcode,primary_site,remove_base_and_mixed,p16,ISH)],
           by.x="PatientID",by.y="patient_barcode")
x[smoking_bin.x=="",smoking_bin.x:=NA_character_]
x[, smoke:=factor(smoking_bin.x,levels=c("Never","Ever"))]
x[, sex:=factor(Sex)]; x[,stage:=factor(Stage)]; x[,alcohol:=factor(alcohol_bin)]
state <- function(z) fifelse(z %in% c("Negative","Positive"),z,NA_character_)
x[,both_disagree:=!is.na(state(p16)) & !is.na(state(ISH)) & state(p16)!=state(ISH)]
x[,hpv_status:=factor(fcase(p16 %in% "Positive"|ISH %in% "Positive","Positive",
                            p16 %in% "Negative"|ISH %in% "Negative","Negative",
                            default=NA_character_),levels=c("Negative","Positive"))]
x[,site:=factor(fcase(grepl("Base of tongue",primary_site),"base_tongue",
                      grepl("pharynx",primary_site),"mixed",
                      grepl("tongue",primary_site),"other_tongue",
                      grepl("Floor of mouth",primary_site),"floor",
                      default="other_oral"))]
sig <- c("Exposure2","Proliferation2","Stromal5")
base <- c("Age","sex","stage","site")
fit_group <- function(d,label,adjust_alcohol=FALSE) {
  vars <- c(base,if(adjust_alcohol) "alcohol")
  d <- d[complete.cases(d[,c("smoke",sig,vars),with=FALSE])]
  stopifnot(nrow(d)>ncol(model.matrix(reformulate(c("smoke",vars),"Exposure2"),d)),
            all(c("Ever","Never") %in% as.character(d$smoke)))
  rbindlist(lapply(sig,function(y){
    f <- lm(reformulate(c("smoke",vars),y),data=d)
    X <- model.matrix(f); h <- hatvalues(f)
    stopifnot(qr(X)$rank==ncol(X),max(h)<1-1e-8)
    w <- residuals(f)/(1-h); bread <- solve(crossprod(X))
    V <- bread %*% crossprod(X*w) %*% bread
    beta <- unname(coef(f)["smokeEver"]); se <- sqrt(V["smokeEver","smokeEver"])
    df <- df.residual(f)
    data.table(cohort=label,signature=y,N=nrow(d),ever=sum(d$smoke=="Ever"),
      never=sum(d$smoke=="Never"),difference_Ever_minus_Never=beta,
      lower95=beta-qt(.975,df)*se,upper95=beta+qt(.975,df)*se,
      p=2*pt(-abs(beta/se),df),adjustment=paste(vars,collapse=" + "))
  }))
}
oral <- x[remove_base_and_mixed %in% TRUE]
assay_negative <- oral[hpv_status=="Negative" & !both_disagree]
strict_both_negative <- oral[p16=="Negative" & ISH=="Negative"]
cat("oral_any_assay_negative_rows=",nrow(assay_negative),"; strict_both_assays_negative_rows=",nrow(strict_both_negative),"\n",sep="")
out <- rbindlist(list(
 fit_group(assay_negative,"Oral, no positive HPV assay; at least one negative assay"),
 fit_group(strict_both_negative,"Oral, p16 and ISH both negative")
),fill=TRUE)
out[,FDR_3_scores:=p.adjust(p,"BH"),by=cohort]
fwrite(out,"hpv_audit/oral_hpv_negative_score_sensitivity.tsv",sep="\t")
fwrite(data.table(cohort=c("Oral, no positive HPV assay; at least one negative assay","Oral, p16 and ISH both negative"),
 N_initial=c(nrow(assay_negative),nrow(strict_both_negative)),
 ever_initial=c(sum(assay_negative$smoke=="Ever",na.rm=TRUE),sum(strict_both_negative$smoke=="Ever",na.rm=TRUE)),
 never_initial=c(sum(assay_negative$smoke=="Never",na.rm=TRUE),sum(strict_both_negative$smoke=="Never",na.rm=TRUE))),
 "hpv_audit/oral_hpv_negative_cohort_counts.tsv",sep="\t")
fwrite(data.table(input=c(paths,"../01_baseline/results/TCGA_HNSC_HPV_cBioPortal.tsv"),
 md5=unname(tools::md5sum(c(paths,"../01_baseline/results/TCGA_HNSC_HPV_cBioPortal.tsv")))),
 "hpv_audit/oral_hpv_negative_input_hashes.tsv",sep="\t")
capture.output(sessionInfo(),file="hpv_audit/oral_hpv_negative_session.txt")

