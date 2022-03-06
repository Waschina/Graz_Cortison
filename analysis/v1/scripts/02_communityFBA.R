library(MicrobiomeGS2)
library(ggplot2)
library(parallel)
library(ggpubr)
library(rstatix)
n.cores <- detectCores()
cl <- makeCluster(min(n.cores-1, 3))
cl_big <- makeCluster(n.cores - 1)

min_rel_abun <- 0.001 # minimum relative abundance of species to be included in commFBA

mic <- new("Microbiome",
           uniq.table.file = "data/clean/asv_tab.tsv",
           model.mapping.file = "analysis/v1/files/asvs_to_HRGM.m8",
           sample.description.file = "data/clean/sample_info.csv",
           uniq.table.format = "R_table"
)
mic <- filter_mapping(mic, method.resolve.multiple = "first")
mic <- create_model_table(mic)
mic <- filter_samples(mic, min.seqs = 2000, max.unclassified = 0.3)


# featch relevant models
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/",
                                 IDs = rownames(mic@model.table[-1,]))

models <- parLapply(cl_big, models,
                    function(mod) {
                      require(sybil)
                      
                      der <- deadEndMetabolites(mod)
                      mod <- rmReact(mod, react = der$der)
                      
                      return(mod)
                    })
stopCluster(cl_big)

clusterExport(cl, c("models","mic","min_rel_abun")) # may take a while...

cFBA_out <- parLapply(cl, 1:ncol(mic@model.table),
                      fun = function(i) {
                        require(MicrobiomeGS2)
                        
                        n <- ncol(mic@model.table)
                        
                        rel_models <- rownames(mic@model.table[which(mic@model.table[,i] > 0),])
                        rel_models <- rel_models[rel_models != "_unclassified"]
                        rel_abun <- mic@model.table[rel_models, i]
                        rel_abun <- rel_abun/sum(rel_abun)
                        rel_abun <- rel_abun[rel_abun >= min_rel_abun]
                        
                        tmpres <- communityFBA_FB(models[names(rel_abun)],
                                                         model.prop = rel_abun,
                                                  cpl_c = 100)
                        tmpres$met.interchange[, gr := tmpres$community.growth]
                        
                        system(paste0("echo \"", i,"/",n," (",
                                      colnames(mic@model.table)[i],", ", 
                                      length(rel_abun)," species, mu = ",
                                      round(tmpres$community.growth, digits = 4),
                                      ") done ", tmpres$solj$stat,
                                      "\""))
                        return(tmpres)
                      })

stopCluster(cl)
names(cFBA_out) <- colnames(mic@model.table)

outflow <- lapply(cFBA_out, function(x) x$met.interchange)
outflow <- rbindlist(outflow, idcol = "sample")
outflow[, o.flux.n := o.flux / gr]
#outflow[, o.flux.n := o.flux]
outflow <- dcast(outflow, sample ~ rxn, value.var = "o.flux.n", fill = 0)
outflow <- melt(outflow, id.vars = "sample", variable.name = "compound", value.name = "flux")
outflow <- outflow[compound != "EX_cpd00001_e0"]
outflow[compound == "EX_cpd00221_e0", compound := "EX_cpd00159_e0"] # merge both lactate forms
outflow <- outflow[, .(flux = sum(flux)), by = c("sample","compound")]
outflow[abs(flux)< 1e-6 & flux != 0, flux := 0]

fwrite(outflow, "analysis/v1/files/outflow.csv")
#saveRDS(cFBA_out, file = "analysis/v1/files/cFBA_out.RDS")


# realte with response info
dt_of <- merge(outflow, mic@sample.description)
setkey(dt_of, "patient-id")
cpdOI <- "EX_cpd00211_e0"


my_comparisons <- list( c("Non-response", "Response"))
p_resp <- ggplot(dt_of[compound == cpdOI], aes(response_lichtiger,
                                               flux,
                                               fill = response_lichtiger,
                                               col = response_lichtiger)) +
  geom_boxplot(width = 0.25, position= position_dodge(0.75),
               alpha = 0.25, outlier.shape = NA) +
  geom_violin(alpha = 0.25,  width = 0.75, position= position_dodge(0.75)) +
  geom_point(pch = 21, size = 0.45, alpha = 0.25,
             position = position_jitterdodge(0.1)) +
  scale_fill_manual(values = c("#3c204b","#375e3a")) +
  scale_color_manual(values = c("#3c204b","#375e3a")) +
  facet_grid(. ~ visit) +
  theme_bw() +
  labs(y = "Predicted butyrate production\n[mmol/gDW]",
       x = "Therapy response") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons,  tip.length = 0)


wilcox.test(dt_of[compound == cpdOI & visit == "V2" & response_lichtiger == "Response", flux],
            dt_of[compound == cpdOI & visit == "V2" & response_lichtiger == "Non-response", flux])


# V1 vs V2 withon responders vs. non-responders
pat_pairs <- dt_of[compound == cpdOI,.N,by = `patient-id`][N == 2, `patient-id`]

stat.test <- dt_of[pat_pairs][compound == cpdOI] %>%
  group_by(response_lichtiger) %>%
  wilcox_test(flux ~ visit, paired = TRUE) %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "visit")
stat.test$y.position <- 0.95 * stat.test$y.position


p_visits <- ggplot(dt_of[pat_pairs][compound == cpdOI],
                   aes(visit,
                       flux,
                       fill = response_lichtiger,
                       col = response_lichtiger,
                       group = `patient-id`)) +
  geom_line(alpha = 0.35) + 
  facet_grid((. ~ response_lichtiger)) +
  geom_point(pch = 21, size = 0.45, alpha = 0.25) +
  scale_fill_manual(values = c("#3c204b","#375e3a")) +
  scale_color_manual(values = c("#3c204b","#375e3a")) +
  theme_bw() +
  labs(y = "Predicted butyrate production\n[mmol/gDW]",
       x = "Visit") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  stat_pvalue_manual(stat.test, label = "{p}", tip.length = 0)
# paired wilcox test for test V1 vs V2
setkey(dt_of, "patient-id")

wilcox.test(dt_of[pat_pairs][compound == cpdOI & visit == "V1" & response_lichtiger == "Response", flux],
            dt_of[pat_pairs][compound == cpdOI & visit == "V2" & response_lichtiger == "Response", flux],
            paired = TRUE)

wilcox.test(dt_of[pat_pairs][compound == cpdOI & visit == "V1" & response_lichtiger == "Non-response", flux],
            dt_of[pat_pairs][compound == cpdOI & visit == "V2" & response_lichtiger == "Non-response", flux],
            paired = TRUE)

p_comb <- egg::ggarrange(p_resp, p_visits, ncol = 2)

ggsave("analysis/v1/plots/cFBA_butyrate.pdf", plot = p_comb, width = 8, height = 4.7)
