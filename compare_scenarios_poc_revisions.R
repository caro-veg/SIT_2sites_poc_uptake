
library(reshape2)
library(ggplot2)
library(grid)
library(cowplot)
library(XLConnect)

path <- "D:\\NumericBehaviourODE\\AdditionaPOCT_scenarios\\"
filename <- "Revisions.Results.xlsx"

wb <- loadWorkbook(paste0(path, filename))
df <- readWorksheet(wb, sheet="POCT_ParCD86N_0.01", header=TRUE)

df$Mutation.Rate <- as.numeric(df$Mutation.Rate)
df$Uptake <- as.factor(df$Uptake)
df$Init.Prev <- factor(df$Init.Prev, levels=c(0, 0.669, 3, 13, 38.6))
df$Res.5 <- as.numeric(df$Res.5)

dl <- split(df, df$Mutation.Rate)

plist <- list()

for(i in 1:length(dl))
{
	p <- ggplot(data=dl[[i]], aes(x=Init.Prev, y=Res.5, fill=Uptake))
	p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=16), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16))  
	p <- p + geom_bar(stat="identity", position=position_dodge()) + coord_flip() 
	p <- p + xlab("Prevalence mutation 1 (%)") + ylab("No. simulations with >5% resistance")
	p <- p + guides(fill=guide_legend(title="Uptake"))
	p <- p + ylim(0, 100)

	plist[[i]] <- p
}

pic <- plot_grid(plotlist=plist, nrow=2, labels=c("a", "b", "c", "d", "e", "f", "g"), label_size=18)
ggsave(filename=paste0(path, "FigSX_poct_1pcGyrAmut.pdf"), dpi=600, plot=pic, width=16, height=9) 


