# complex_networks_plots.R
# plot the pretty pictures used in the paper.


# Show the various types of protein targets and their freqs
bar_plot_drugtargets <- function(){
  tempnames <- table(drug_targets$TargetClass)
  mp <- barplot(table(drug_targets$TargetClass),col="grey50",space=1,
                ylab="Number of drug targets (proteins)",ylim=c(0,6000),main="",xlab = "",xaxt = "n")
                
  mytable <- length(table(drug_targets$TargetClass))
  end_point <- 0.5 + (mytable) + (mytable)-1
  
  text(seq(1.5,end_point,by=2), par("usr")[3]-0.25, # -0.25
       srt = 60, adj= 1, xpd = TRUE,
       labels = paste(names(tempnames)),cex=1.0)
}


