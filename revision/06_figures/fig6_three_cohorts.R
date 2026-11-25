set.seed(20261121)
old <- parse("revised_figures.R",encoding="UTF-8")
for (expr in old) {
  if (is.call(expr) && identical(expr[[1]],as.name("<-")) &&
      identical(expr[[2]],as.name("fig6"))) break
  eval(expr,envir=.GlobalEnv)
}
src <- "../09_GSE85446"
cc <- fread(file.path(src,"Final8_discrimination_cc.tsv"))
au <- fread(file.path(src,"Final8_discrimination_au.tsv"))
hr <- fread(file.path(src,"Final8_clinical_adjustment.tsv"))
stopifnot(cc$valid==2000L,identical(au$months,c(12L,36L,60L)),nrow(hr)==1L)
label <- "GSE85446, frozen Final8"
c_data[,cohort:=as.character(cohort)]
c_data <- rbindlist(list(c_data,data.table(cohort=label,C=cc$C,lower=cc$lower95,
          upper=cc$upper95,group="External Final8")),fill=TRUE)
c_data[,cohort:=factor(cohort,levels=rev(cohort))]
p_d$data <- c_data
auc_data[,cohort:=as.character(cohort)]
auc_data <- rbindlist(list(auc_data,data.table(cohort=label,years=levels(auc_data$years),
          AUC=au$AUC,lower=au$lower95,upper=au$upper95)))
auc_data[,cohort:=factor(cohort,levels=unique(cohort))]
auc_data[,years:=factor(years,levels=c("1 year","3 years","5 years"))]
pd <- position_dodge(width=.37)
p_e <- ggplot(auc_data,aes(years,AUC,color=cohort,group=cohort,shape=cohort)) +
  geom_hline(yintercept=.5, linetype=2, color="grey65",linewidth=.45) +
  geom_line(linewidth=.55,position=pd) +
  geom_errorbar(aes(ymin=lower,ymax=upper),position=pd,width=.08,linewidth=.48) +
  geom_point(size=2,position=pd) +
  scale_colour_manual(values=c(amber,red,blue,"#327864"),
                      labels=c("TCGA pipeline","GSE41613","GSE42743","GSE85446")) +
  scale_shape_manual(values=c(16,17,15,18),
                     labels=c("TCGA pipeline","GSE41613","GSE42743","GSE85446")) +
  scale_y_continuous(limits=c(.35,1),breaks=seq(.4,1,.2)) +
  labs(x=NULL,y="Time-dependent AUC (95% CI)") +base_theme +
  guides(color=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))
adj_data[,cohort:=as.character(cohort)]
adj_data <- rbindlist(list(adj_data,data.table(cohort="GSE85446, adjusted",
        HR=hr$HR_per_SD,lower=hr$low95,upper=hr$high95)),fill=TRUE)
adj_data[,cohort:=factor(cohort,levels=rev(cohort))]
p_f$data <- adj_data
p_f <- p_f + scale_colour_manual(values=c("GSE41613, adjusted"=red,
        "GSE42743, adjusted"=blue,"GSE85446, adjusted"="#327864"))
stopifnot(nrow(c_data)==6L,nrow(auc_data)==12L,nrow(adj_data)==3L)
fig <- (p_a|p_b)/(p_c|p_d)/(p_e|p_f) + plot_annotation(tag_levels="A") &
       theme(plot.tag=element_text(face="bold",size=14,color="black"))
ggsave("Fig6_three_cohorts.pdf",fig,width=12,height=9,device=cairo_pdf,bg="white")
set.seed(20261121)
ggsave("Fig6_three_cohorts.png",fig,width=12,height=9,dpi=300,device=grDevices::png,type="cairo",bg="white")
