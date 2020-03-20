
library(reshape2)
library(ggplot2)


path <- "D:\\NumericBehaviourODE\\AdditionaPOCT_scenarios\\"
filename <- "SIT_2Sites-6.02e-008-0.0054-0.0833-1.778-1-1-0.1-2.45e-008-2.45e-008-1.23e-007-1.23e-007-1478000-13288-8492-220-0-1825-100.csv"

uptake <- seq(0, 0.9, 0.1)
mutRates <- c("2.45e-008", "4.9e-008", "1.23e-007", "2.45e-007", "2.45e-006", "2.45e-005", "7.95e-005")

filenames <- list()
params <- list()
pos <- 0

for(i in 1:7)
{
	for(j in 1:10)
	{
		pos = pos + 1
		file <- paste0("SIT_2Sites-6.02e-008-0.0054-0.0833-1.778-1-1-", toString(uptake[j]), "-2.45e-008-2.45e-008-", mutRates[i], "-", mutRates[i], "-1478000-21780-0-220-0-1825-100.csv")
		filenames[[pos]] <- file
		params[[pos]] <- list(file, as.numeric(mutRates[i]), uptake[j], 0)
	}
}


for(i in 1:7)
{
	for(j in 1:10)
	{
		pos = pos + 1
		file <- paste0("SIT_2Sites-6.02e-008-0.0054-0.0833-1.778-1-1-", toString(uptake[j]), "-2.45e-008-2.45e-008-", mutRates[i], "-", mutRates[i], "-1478000-21633-147-220-0-1825-100.csv")
		filenames[[pos]] <- file
		params[[pos]] <- list(file, as.numeric(mutRates[i]), uptake[j], 0.669)
	}
}


for(i in 1:7)
{
	for(j in 1:10)
	{
		pos = pos + 1
		file <- paste0("SIT_2Sites-6.02e-008-0.0054-0.0833-1.778-1-1-", toString(uptake[j]), "-2.45e-008-2.45e-008-", mutRates[i], "-", mutRates[i], "-1478000-21120-660-220-0-1825-100.csv")
		filenames[[pos]] <- file
		params[[pos]] <- list(file, as.numeric(mutRates[i]), uptake[j], 3)
	}
}


for(i in 1:7)
{
	for(j in 1:10)
	{
		pos = pos + 1
		file <- paste0("SIT_2Sites-6.02e-008-0.0054-0.0833-1.778-1-1-", toString(uptake[j]), "-2.45e-008-2.45e-008-", mutRates[i], "-", mutRates[i], "-1478000-18920-2860-220-0-1825-100.csv")
		filenames[[pos]] <- file
		params[[pos]] <- list(file, as.numeric(mutRates[i]), uptake[j], 13)
	}
}


for(i in 1:7)
{
	for(j in 1:10)
	{
		pos = pos + 1
		file <- paste0("SIT_2Sites-6.02e-008-0.0054-0.0833-1.778-1-1-", toString(uptake[j]), "-2.45e-008-2.45e-008-", mutRates[i], "-", mutRates[i], "-1478000-13288-8492-220-0-1825-100.csv")
		filenames[[pos]] <- file
		params[[pos]] <- list(file, as.numeric(mutRates[i]), uptake[j], 38.6)
	}
}




############################################################################
# GET NUMBER OF RUNS IN WHICH RESISTANCE PREVALENCE REACHED 5%
############################################################################


outList <- list()
for(i in 1:length(filenames))
{
	df <- read.csv(paste0(path, filenames[[i]]), header=TRUE)
	df$repeat. <- as.factor(df$repeat.)
	split.df <- split(df, f=df$repeat.)

	
	r11 <- rep(0, length(unique(df$repeat.)))
	r10 <- rep(0, length(unique(df$repeat.)))
	r01 <- rep(0, length(unique(df$repeat.)))

	for(j in 1:length(split.df))
	{
		I_total <- apply(split.df[[j]][nrow(split.df[[j]]), c(4:15)], 1, sum)[[1]]
	
		T_11 <- split.df[[j]][nrow(split.df[[j]]), "T11"]
		T_10 <- split.df[[j]][nrow(split.df[[j]]), "T10"]
		T_01 <- split.df[[j]][nrow(split.df[[j]]), "T01"]

		T_alt11 <- split.df[[j]][nrow(split.df[[j]]), "Ta11"]
		T_alt10 <- split.df[[j]][nrow(split.df[[j]]), "Ta10"]
		T_alt01 <- split.df[[j]][nrow(split.df[[j]]), "Ta01"]

		I_11 <- split.df[[j]][nrow(split.df[[j]]), "I11"]
		I_10 <- split.df[[j]][nrow(split.df[[j]]), "I10"]
		I_01 <- split.df[[j]][nrow(split.df[[j]]), "I01"]

		IT_11 <- I_11 + T_11 + T_alt11
		IT_10 <- I_10 + T_10 + T_alt10
		IT_01 <- I_01 + T_01 + T_alt01

	
		if(I_total > 0)
		{
			r11[j] <- IT_11 / I_total * 100
			r10[j] <- IT_10 / I_total * 100
			r01[j] <- IT_01 / I_total * 100
		}	
		else
		{
			r11[j] <- 0
			r10[j] <- 0
			r01[j] <- 0
		}
	}
	
	outList[[i]] <- c(params[[i]], length(r11[r11>1]), length(r11[r11>5]), length(r11[r11>10]))
}

out <- do.call(rbind, outList)
out <- as.data.frame(out)
names(out) <- c("Scenario", "Mutation.Rate", "Uptake", "Init.Prev", "Res.1", "Res.5", "Res.10")

wb.out <- loadWorkbook(paste0(path, "Revisions.Results.xlsx"), create=TRUE)
createSheet(wb.out, name="POCT_ParCD86N_0.01")
writeWorksheet(wb.out, out, sheet="POCT_ParCD86N_0.01")
saveWorkbook(wb.out)




