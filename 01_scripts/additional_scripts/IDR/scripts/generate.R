if (!require("reshape2"))
	install.packages("reshape2", repos='https://ftp.gwdg.de/pub/misc/cran/')
library(reshape2)

data <- read.table("temp/all_tools", header= TRUE)

idr_input <- dcast(data, Location ~ Tool, value.var = "QM", fun.aggregate = max)

idr_input[idr_input=="-Inf"]<- 0 #(-.Machine$integer.max)

write.table(idr_input, "temp/idr_Input")

