args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)



visualize_jxn_and_exon_locations <- function(df, pos,output_file,title) {
	miny <- pos - 30
	maxy <- pos + 30
	num_rows <- dim(df)[1]
	for (row_num in 1:num_rows) {
		if (df[row_num,1] < miny) {
			df[row_num,1] = miny
		}
		if (df[row_num,1] > maxy) {
			df[row_num,1] = maxy
		}
		if (df[row_num,2] < miny) {
			df[row_num,2] = miny
		}

		if (df[row_num,2] > maxy) {
			df[row_num, 2] = maxy
		}

	}


	df2 <- data.frame(position=c(df$start, df$end), y = c(df$y, df$y), version=as.factor(c(as.character(df$version), as.character(df$version))), group =as.factor(c(as.character(df$group), as.character(df$group))))
	vecy <- paste0(as.character(df$version),df$y)
	ploty <- ggplot(df2, aes(x=position, y=y, colour = factor(version,levels=c("jxn","exon")), group=group)) + geom_line() +
		theme_bw() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=10),axis.text=element_text(size=9),panel.grid.minor = element_blank(),panel.grid.major=element_blank(),, axis.line = element_line(colour = "black"), legend.text = element_text(size=9), legend.title = element_text(size=10)) +
		geom_vline(xintercept=pos, colour="black", size=1) + 
		labs(x = "Position (BP)", y = "",colour="",title=title)

	ggsave(ploty, file=output_file,width = 16,height=10.5,units="cm")


}


variant_position_enrichment_debug_dir <- args[1]
pvalue_threshold <- args[2]
distance <- args[3]


input_file <- paste0(variant_position_enrichment_debug_dir, "outlier_distance_to_ss_differences_between_methods_distance_", distance, "_pvalue_thresh_0", pvalue_threshold, ".txt")

data <- read.table(input_file, header=TRUE)
num_rows <- dim(data)[1]



for (row_num in 1:num_rows) {
	line_start <- c()
	line_end <- c()
	version <- c()	
	y_positions <- c()
	groups <- c()
	variant_position <- data[row_num, 1]
	jxn_vec <- as.character(data[row_num, 4])
	exon_vec <- as.character(data[row_num, 5])
	exon_arr <- strsplit(exon_vec, ",")[[1]]
	jxn_arr <- strsplit(jxn_vec, ",")[[1]]
	y_pos <- 1
	for (index in 1:length(exon_arr)) {
		exon <- exon_arr[index]
		if (exon != "none") {
			exon_info <- strsplit(exon,':')[[1]]
			start <- exon_info[1]
			end <- exon_info[2]
			line_start <- c(line_start, as.numeric(start))
			line_end <- c(line_end, as.numeric(end))
			version <- c(version, "exon")
			y_positions <- c(y_positions, y_pos)
			groups <- c(groups, paste0("group", y_pos))
			y_pos <- y_pos + 1
		}
	}
	for (index in 1:length(jxn_arr)) {
		exon <- jxn_arr[index]
		exon_info <- strsplit(exon,':')[[1]]
		start <- exon_info[1]
		end <- exon_info[2]
		line_start <- c(line_start, as.numeric(start))
		line_end <- c(line_end, as.numeric(end))
		version <- c(version, "jxn")
		y_positions <- c(y_positions, y_pos)
		groups <- c(groups, paste0("group", y_pos))
		y_pos <- y_pos + 1
	}

	output_file <-paste0(variant_position_enrichment_debug_dir, "debug_outlier_distance_to_ss_differences_between_methods_distance_", distance, "_pvalue_thresh_0", pvalue_threshold, "_", row_num,".pdf")
	df <- data.frame(start=line_start, end=line_end, version=factor(version), y=y_positions, group=factor(groups))
	title=paste0("observed=",data[row_num,2], " / exon=", data[row_num,3])
	visualize_jxn_and_exon_locations(df, variant_position, output_file,title)
}


