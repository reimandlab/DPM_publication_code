# perform DPM and Brown's p-value merging on simulated datasets

setwd("~/")
library(ActivePathways)
library(ggplot2)
library(patchwork)
library(ggrastr)

set.seed(54321)

count_res = function(alpha, res, method, direction) {
	
	n_res = sum(res <= alpha)
	n_FDR = sum(p.adjust(res, method = "fdr") <= alpha)
	
	data.frame(alpha, n_res, n_FDR, method, direction, stringsAsFactors = FALSE)
}

master_results_table = NULL


for (pval_distr in c("uniform_uniform", "signf_signf", "uniform_signf")) {
	
	# the second set of Z-scores differs from the first set by some added noise
	a_bit_of_noise = rnorm(mean = 0, sd = 0.2, n = 10000)

	if (pval_distr == "uniform_uniform") {
		# sample z-scores for first set of P-values
		z1 = rnorm(mean = 0, sd = 1, n = 10000)
		z2 = z1 + a_bit_of_noise
		z3 = rnorm(mean = 0, sd = 1, n = 10000)
		# Z-scores to p-values
		p1 = pnorm(z1)
		p2 = pnorm(z2)
		p3 = pnorm(z3)
	}
	if (pval_distr == "signf_signf") {
		# sample z-scores for first set of P-values
		z1 = rnorm(mean = 1, sd = 1, n = 10000)
		z2 = z1 + a_bit_of_noise
		z3 = rnorm(mean = 1, sd = 1, n = 10000)
		# Z-scores to p-values
		p1 = pnorm(z1, lower.tail = FALSE)
		p2 = pnorm(z2, lower.tail = FALSE)
		p3 = pnorm(z3, lower.tail = FALSE)
	}
	if (pval_distr == "uniform_signf") {
		# sample z-scores: one uniform, one highly-significant
		z1 = rnorm(mean = 0, sd = 1, n = 10000)
		z3 = rnorm(mean = 1, sd = 1, n = 10000)
		# Z-scores to p-values
		p1 = pnorm(z1, lower.tail = FALSE)
		p3 = pnorm(z3, lower.tail = FALSE)
	}	
	
	for (analysis_type in c("dependent", "independent")) {
		
		# this option does not make sense
		if (analysis_type == "dependent" && pval_distr == "uniform_signf") {
			next()
		}
	
		# p-values for merging in a matrix
		p_mat = cbind(p1 = p1, p2 = p2)
		if (analysis_type == "independent") {
			p_mat = cbind(p1 = p1, p2 = p3)		
		}
		
		corval = cor(p_mat[,1], p_mat[,2])#, method = "spearman")
		cat(analysis_type, "", pval_distr, ": ", corval, "\n")
		
		# directions for merging in a matrix
		d_mat = matrix(1, 10000, 2)
		colnames(d_mat) = c("p1", "p2")
		
		# create another set of directions where half of directions are mismatching
		d_mat_mixed = d_mat
		d_mat_mixed[,1] = rbinom(10000, 1, 0.5) * -2 + d_mat_mixed[,1]
		
		# constraints_vector: (1,1) means that expected directions are matching, (-1,1) means that mismatches are expected
		CV_match = c("p1" = 1, "p2" = 1)
		CV_mismatch = c("p1" = -1, "p2" = 1)
		
		# merge P-values such that all input genes have matching directions
		Pmerge_DPM_agree = merge_p_values(p_mat, scores_direction = d_mat, 
				constraints_vector = CV_match, method = "DPM")
		Pmerge_Strube_agree = merge_p_values(p_mat, scores_direction = d_mat, 
				constraints_vector = CV_match, method = "Strube_directional")
		# merge P-values such that all input genes have mismatching directions
		Pmerge_DPM_disagree = merge_p_values(p_mat, scores_direction = d_mat, 
				constraints_vector = CV_mismatch, method = "DPM")
		Pmerge_Strube_disagree = merge_p_values(p_mat, scores_direction = d_mat, 
				constraints_vector = CV_mismatch, method = "Strube_directional")
		# merge P-values such that 50% input genes have matching directions and 50% have mismatching directions
		Pmerge_DPM_mixed = merge_p_values(p_mat, scores_direction = d_mat_mixed, 
				constraints_vector = CV_match, method = "DPM")
		Pmerge_Strube_mixed = merge_p_values(p_mat, scores_direction = d_mat_mixed, 
				constraints_vector = CV_match, method = "Strube_directional")

		alphas = c(0.01, 0.05, 0.1, 0.2)
		
		res_DPM_mixed = do.call(rbind, lapply(alphas, count_res, Pmerge_DPM_mixed, "DPM", "mixed_50%"))
		res_Strube_mixed = do.call(rbind, lapply(alphas, count_res, Pmerge_Strube_mixed, "Strube_directional", "mixed_50%"))
		res_DPM_agree = do.call(rbind, lapply(alphas, count_res, Pmerge_DPM_agree, "DPM", "agree"))
		res_Strube_agree = do.call(rbind, lapply(alphas, count_res, Pmerge_Strube_agree, "Strube_directional", "agree"))
		res_DPM_disagree = do.call(rbind, lapply(alphas, count_res, Pmerge_DPM_disagree, "DPM", "disagree"))
		res_Strube_disagree = do.call(rbind, lapply(alphas, count_res, Pmerge_Strube_disagree, "Strube_directional", "disagree"))
		
		res_combined = rbind(
				res_DPM_mixed, res_Strube_mixed,
				res_DPM_agree, res_Strube_agree,
				res_DPM_disagree, res_Strube_disagree)
		res_combined$direction = factor(res_combined$direction, levels = c("agree", "mixed_50%", "disagree"))
		
		ggplt_alpha = ggplot(res_combined, aes(factor(alpha), n_res, fill = method)) + 
				geom_bar(stat = "identity", position = "dodge") + 
				facet_wrap(~direction) + 
				ggtitle(paste0("Directional merging of P-values"), paste(analysis_type, ", ", pval_distr, ", rho=", signif(corval, 2))) + 
				scale_fill_manual(values = c(DPM = "orange", Strube_directional = "darkblue")) + 
				theme_bw()
				
		dfr_p = data.frame(p_mat, 
				Pmerge_DPM_disagree, Pmerge_Strube_disagree, 
				Pmerge_DPM_agree, Pmerge_Strube_agree, 
				stringsAsFactors = FALSE)
				
		dfr_p$mergedP_signf_disagree = "none"
		dfr_p[dfr_p$Pmerge_DPM_disagree < 0.05, "mergedP_signf_disagree"] = "DPM"
		dfr_p[dfr_p$Pmerge_Strube_disagree < 0.05, "mergedP_signf_disagree"] = "Strube"
		dfr_p[dfr_p$Pmerge_Strube_disagree < 0.05 & dfr_p$Pmerge_DPM_disagree < 0.05, "mergedP_signf_disagree"] = "BOTH"
		dfr_p$mergedP_signf_disagree = factor(dfr_p$mergedP_signf_disagree, levels = c("DPM", "Strube", "none", "BOTH"))

		dfr_p[dfr_p$Pmerge_DPM_agree < 0.05, "mergedP_signf_agree"] = "DPM"
		dfr_p[dfr_p$Pmerge_Strube_agree < 0.05, "mergedP_signf_agree"] = "Strube"
		dfr_p[dfr_p$Pmerge_Strube_agree < 0.05 & dfr_p$Pmerge_DPM_agree < 0.05, "mergedP_signf_agree"] = "BOTH"
		dfr_p$mergedP_signf_agree = factor(dfr_p$mergedP_signf_agree, levels = c("DPM", "Strube", "none", "BOTH"))
		
		plot_subtitle = paste("n_p1=", sum(dfr_p$p1 < 0.05), "; n_p2=", sum(dfr_p$p2 < 0.05), "(p<0.05)")
		p1_label = paste0("P1, -log10, ", gsub("(.+)_(.+)", "\\1", pval_distr))
		p2_label = paste0("P2, -log10, ", gsub("(.+)_(.+)", "\\2 ", pval_distr))
		
		ggplt_distr_disagree = ggplot(dfr_p, aes(-log10(p1), -log10(p2), color = mergedP_signf_disagree)) + 
				geom_point_rast(alpha = 0.333, size = 1) + 
				ggtitle(paste0("P scatter: disagree"), plot_subtitle) + 
				scale_color_manual("merged P<0.05,", drop = FALSE,
						values = c(DPM = "orange", Strube = "darkblue", BOTH = "darkred", none = "lightgrey")) + 
				scale_x_continuous(p1_label, limits = c(0,6)) +
				scale_y_continuous(p2_label, limits = c(0,6)) +
				theme_bw()

		ggplt_distr_agree = ggplot(dfr_p, aes(-log10(p1), -log10(p2), color = mergedP_signf_agree)) + 
				geom_point_rast(alpha = 0.333, size = 1) + 
				ggtitle(paste0("P scatter: agree"), plot_subtitle) + 
				scale_color_manual("merged P<0.05", drop = FALSE,
						values = c(DPM = "orange", Strube = "darkblue", BOTH = "darkred", none = "lightgrey")) + 
				scale_x_continuous(p1_label, limits = c(0,6)) +
				scale_y_continuous(p2_label, limits = c(0,6)) +
				theme_bw()
				
		ggplt_combined = ggplt_alpha + (ggplt_distr_agree / ggplt_distr_disagree) +  plot_layout(widths = c(2, 1))
		
		fname = paste0("003R_Directional_Pmerge_", analysis_type, "__", pval_distr, "__DPM_Strube.pdf")
		ggsave(ggplt_combined, file = fname, width = 8, height = 5)
		system("sleep 3")
		system(paste("open ", fname))
		
		res_combined1 = cbind(pval_distr, analysis_type, res_combined)
		master_results_table = rbind(master_results_table, res_combined1)
	}
}

save(master_results_table, file = "master_results_table.rsav")
write.csv(master_results_table, file = "master_results_table.csv")
system(paste("open master_results_table.csv"))




#dependent  uniform_uniform :  0.9783065
#independent  uniform_uniform :  -0.0009814493
#dependent  signf_signf :  0.9771874
#independent  signf_signf :  0.0001891833
#independent  uniform_signf :  -0.004741001
