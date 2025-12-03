Pri_AFR_mctre<-ape::read.tree("/Users/zsong/Downloads/Primates_Dryad_no_scripts/timetree/combined.trees") # Wisniewski Anna L., Lloyd Graeme T. and Slater Graham J. 2022 Extant species fail to estimate ancestral geographical ranges at older nodes in primate phylogenyProc. R. Soc. B.28920212535 http://doi.org/10.1098/rspb.2021.2535
Pri_AFR_mctre1<-Pri_AFR_mctre[[1]]

Pri_Sur_AFR<-read.csv("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/Surv_AFR_Primate.csv")


which(Pri_Sur_AFR$tip_label %ni% Pri_AFR_mctre[[1]]$tip.label)

random_numbers<- sample(1:1000, 100)
Pri_AFR_mctre<-Pri_AFR_mctre[random_numbers]
library(ape)
for(i in 1: 100){
  Pri_AFR_mctre[[i]]<-drop.tip(Pri_AFR_mctre[[i]], Pri_AFR_mctre[[i]]$tip.label[! Pri_AFR_mctre[[i]]$tip.label %in% Pri_Sur_AFR$tip_label ], root.edge = F, rooted = is.rooted(Pri_AFR_mctre[[i]]))
  
}

#write.tree(Pri_AFR_mctre,"/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/Pri_AFR_mctre.tre")
Pri_ASR_MCCtre<- read.nexus("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/Pri_ASR_MCCtre.tre") 
Pri_Sur_AFR<-read.csv("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/Surv_AFR_Primate.csv")

#Caculate Residual Brain size

rownames(Pri_Sur_AFR)<-Pri_Sur_AFR$tip_label

Pri_Sur_AFR$Residual_BS<-NA
a<-lm(log10(Br.ADF..g.)~log10(WT.ADF..kg.*1000),Pri_Sur_AFR )
summary(a)
Pri_Sur_AFR$Residual_BS<-a$residuals


Pri_Sur_AFR$Residual_Adult_LifeSpan<-NA
a<-lm(log10( Pri_Sur_AFR$Max_lifespan - Pri_Sur_AFR$AFR..yr.) ~log10(WT.ADF..kg.*1000),Pri_Sur_AFR )
summary(a)
Pri_Sur_AFR$Residual_Adult_LifeSpan<-a$residuals

#data transferred
library(car)
library(phylolm)
#
Pri_Sur_AFR$Body<-scale( log10(Pri_Sur_AFR$WT.ADF..kg.), center = T)
Pri_Sur_AFR$Brain<-scale( log10(Pri_Sur_AFR$Br.ADF..g.), center = T)
Pri_Sur_AFR$AFR<-scale( log10(Pri_Sur_AFR$AFR..yr.), center = T)
Pri_Sur_AFR$Residual_BS1<-scale( Pri_Sur_AFR$Residual_BS, center = T)
Pri_Sur_AFR$Surv_AFR<-scale( logit(Pri_Sur_AFR$Surv..AFR), center = T)
Pri_Sur_AFR$Surv_1st<-scale(logit( Pri_Sur_AFR$Surv.1st.yr), center = T)

Pri_Sur_AFR$Max_AdultLifeSpan<-scale(log10( Pri_Sur_AFR$Max_lifespan - Pri_Sur_AFR$AFR..yr.), center = T)



#model selection for Surv_AFR

Pri_Sur_AFR_ResiBS_model<- phylolm(  Surv_AFR ~  Residual_BS1 + Body+ AFR , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01)
vif.phyloglm( Pri_Sur_AFR_ResiBS_model)

summary(Pri_Sur_AFR_ResiBS_model)
library(MuMIn)
dd<-dredge(Pri_Sur_AFR_ResiBS_model)
subset(dd, delta < 5)
write.csv(as.data.frame(subset(dd)), "/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/modelResiBS2.csv")
summary(model.avg(dd, subset = delta < 2))

Coeffs(phylolm(Surv_AFR  ~  Residual_BS1 + Body , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))


##Brain log
Pri_Sur_AFR_BS_model<- phylolm(Surv_AFR  ~  Brain + Body+ AFR , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre, model="lambda", lower.bound = 0.01)
vif.phyloglm(Pri_Sur_AFR_BS_model)

summary(Pri_Sur_AFR_BS_model)
library(MuMIn)
dd1<-dredge(Pri_Sur_AFR_BS_model)
subset(dd1, delta < 5)
write.csv(as.data.frame(subset(dd1)), "/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/modelBS2.csv")
summary(model.avg(dd1, subset = delta < 2))

Coeffs(phylolm(Surv_AFR  ~  Brain + Body , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))

library(ggplot2)
pdf(file = paste("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/", "SurAFR ~ Brainlog.pdf", sep = ""),width = 8, height = 6)

ggplot(Pri_Sur_AFR, 
       aes(x = log(Br.ADF..g.) , y = Surv..AFR  )) +
  geom_point(aes(size = log(WT.ADF..kg.)), alpha = 0.8)+
  geom_smooth(method = "lm", fill = NA, color="black" ) +
  
  theme(legend.position="right",
        text = element_text(size=20), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black")) +
  labs(size = "Body size (log)",
       x = "Brain size (log)",
       y = "Survival rate until AFR")
dev.off()

pdf(file = paste("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/", "SurAFR ~ BrainResi.pdf", sep = "") ,width = 8, height = 6)

ggplot(Pri_Sur_AFR, 
       aes(x = Residual_BS , y = Surv..AFR  )) +
  geom_point(aes(size = log(WT.ADF..kg.)), alpha = 0.8)+
  geom_smooth(method = "lm", fill = NA, color="black" ) +
  
  theme(legend.position="right",
        text = element_text(size=20), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black")) +
  labs(size = "Body size (log)",
       x = "Brain size (residual)",
       y = "Survival rate until AFR")
dev.off()


#Surival 1 year
Pri_Sur_AFR_ResiBS_model_1srt<- phylolm(Surv_1st  ~  Residual_BS1 + Body, data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01)
Coeffs(Pri_Sur_AFR_ResiBS_model_1srt)
summary(Pri_Sur_AFR_ResiBS_model_1srt)
vif.phyloglm(Pri_Sur_AFR_ResiBS_model_1srt)
library(MuMIn)
dd_1st<-dredge(Pri_Sur_AFR_ResiBS_model_1srt)
subset(dd_1st, delta < 5)
write.csv(as.data.frame(subset(dd_1st)), "/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/modelResiBS_1st2.csv")


Pri_Sur_AFR_BS_model_1st<- phylolm(Surv_1st  ~  Brain + Body , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01)
Coeffs(Pri_Sur_AFR_BS_model_1st)
summary(Pri_Sur_AFR_BS_model_1st)
vif.phyloglm(Pri_Sur_AFR_BS_model_1st)
library(MuMIn)
dd_1st_2<-dredge(Pri_Sur_AFR_BS_model_1st)
subset(dd_1st_2, delta < 5)
write.csv(as.data.frame(subset(dd_1st_2)), "/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/modelLogBS_1st2.csv")




pdf(file = paste("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/", "Sur1st ~ BrainResi1.pdf", sep = ""),width = 8, height = 6)

ggplot(Pri_Sur_AFR, 
       aes(x = Residual_BS , y = Surv.1st.yr  )) +
  geom_point(aes(size = log(WT.ADF..kg.)), alpha = 0.8)+
  geom_smooth(method = "lm", fill = NA, color="black" ) +
  
  theme(legend.position="right",
        text = element_text(size=20), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black")) +
  labs(size = "Body size (log)",
       x = "Brain size (residual)",
       y = "Survival rate until age 1")
dev.off()




pdf(file = paste("/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/ModelSel/", "Sur1st ~ Brainlog.pdf", sep = ""),width = 8, height = 6)

ggplot(Pri_Sur_AFR, 
       aes(x = log(Br.ADF..g.)  , y = Surv.1st.yr  )) +
  geom_point(aes(size = log(WT.ADF..kg.)), alpha = 0.8)+
  geom_smooth(method = "lm", fill = NA, color="black" ) +
  
  theme(legend.position="right",
        text = element_text(size=20), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black")) +
  labs(size = "Body size (log)",
       x = "Brain size (log)",
       y = "Survival rate until age 1")
dev.off()


############################################################
############################################################
############################################################

##

summary(phylolm(Surv_AFR  ~  Brain+Body , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))
summary(phylolm(Surv_AFR  ~  Residual_BS1+Body , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))
summary(phylolm(Surv_AFR  ~  Residual_BS1+AFR , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))
summary(phylolm(Surv_AFR  ~AFR + Body + Residual_BS1 , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))
summary(phylolm(Surv_AFR  ~AFR + Body + Brain , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))
summary(phylolm(Surv_AFR  ~AFR , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))
summary(phylolm(Surv_AFR  ~Residual_BS1 , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))



summary(phylolm(Surv_1st ~ Residual_BS1 , data=Pri_Sur_AFR, phy=Pri_ASR_MCCtre,model="lambda", lower.bound = 0.01))


nrepo<-3

inv.phylo <- inverseA((Pri_ASR_MCCtre) ,nodes="TIPS",scale=TRUE)

prior2 <- list(
  R = list(V = diag(nrepo), nu = nrepo-1+0.002),  # 残差协方差，因变量 2 个
  G = list(G1 = list(V = diag(nrepo), nu = nrepo-1+0.002))  # 组间方差
)

MRBPMM_Pri_Sur_AFR <- MCMCglmm(
  cbind(Body,Residual_BS1, AFR) ~ trait:Surv_AFR-1 ,  # 多因变量
  random = ~ us(trait):tip_label,  # 允许因变量间协方差
  rcov = ~ us(trait):units,  # 误差协方差结构
  ginverse = list(tip_label = inv.phylo$Ainv),  # 逆协方差矩阵
  prior = prior2,
  family = c("gaussian","gaussian" ,"gaussian" ),
  data = Pri_Sur_AFR,
  nitt = 50000, burnin = 10000, thin = 40  # MCMC 采样参数
  
)

summary(MRBPMM_Pri_Sur_AFR)



MRBPMM_Pri_Sur_AFR_Log <- MCMCglmm(
  cbind(Body,Brain, AFR) ~ trait:Surv_AFR-1 ,  # 多因变量
  random = ~ us(trait):tip_label,  # 允许因变量间协方差
  rcov = ~ us(trait):units,  # 误差协方差结构
  ginverse = list(tip_label = inv.phylo$Ainv),  # 逆协方差矩阵
  prior = prior2,
  family = c("gaussian","gaussian" ,"gaussian" ),
  data = Pri_Sur_AFR,
  nitt = 50000, burnin = 10000, thin = 40  # MCMC 采样参数
  
)

summary(MRBPMM_Pri_Sur_AFR_Log)

nrepo<-2


prior2 <- list(
  R = list(V = diag(nrepo), nu = nrepo-1+0.002),  # 残差协方差，因变量 2 个
  G = list(G1 = list(V = diag(nrepo), nu = nrepo-1+0.002))  # 组间方差
)
MRBPMM_Pri_Sur_FR_Log <- MCMCglmm(
  cbind(Body,Brain) ~ trait:Surv_1st-1 ,  # 多因变量
  random = ~ us(trait):tip_label,  # 允许因变量间协方差
  rcov = ~ us(trait):units,  # 误差协方差结构
  ginverse = list(tip_label = inv.phylo$Ainv),  # 逆协方差矩阵
  prior = prior2,
  family = c("gaussian", "gaussian" ),
  data = Pri_Sur_AFR,
  nitt = 50000, burnin = 10000, thin = 40  # MCMC 采样参数
  
)

summary(MRBPMM_Pri_Sur_FR_Log)

############################################################
############################################################
############################################################


##Max_lifespan
Pri_Sur_AFR$Max_lifespan<-NA
Pri_Sur_AFR$Max_lifespan_Ref<-NA
for(i in 1: nrow(Pri_Sur_AFR)){
  Pri_Sur_AFR[i,]$Max_lifespan<-  PrimateLHAllDataCompl[which(PrimateLHAllDataCompl$tip_label == Pri_Sur_AFR[i,]$tip_label ),]$Max..Lifespan..y.
  Pri_Sur_AFR[i,]$Max_lifespan_Ref<-"Kopie von PrimateLHAllDataCompl_J_02.12.csv"
}

Pri_Sur_AFR[which(Pri_Sur_AFR$tip_label=="Gorilla_beringei"),]$Max_lifespan<-50.45
Pri_Sur_AFR[which(Pri_Sur_AFR$tip_label=="Gorilla_beringei"),]$Max_lifespan_Ref<-"Myhrvold2015"


summary(phylolm(Surv_AFR ~ Residual_BS1 +Max_AdultLifeSpan, Pri_Sur_AFR , Pri_ASR_MCCtre, model = "lambda"))
summary(phylolm(Surv_1st ~ Residual_BS1 +Max_AdultLifeSpan, Pri_Sur_AFR , Pri_ASR_MCCtre, model = "lambda"))


##path analysis
library(phylopath)
library(igraph)



# Path analysis
m_PriSuv <- define_model_set(
  ##1##
  # Body -> Brain , Brain -> AFR;
  # Body !-> Sur
  m1  = c( Brain ~ Body   , 
           AFR ~ Brain+  Body    ,
           Surv_AFR ~  Brain+ AFR 
  ),
  
  m2  = c( Brain ~ Body   , 
           AFR ~ Brain    ,
           Surv_AFR ~  Brain+ AFR 
  ),
  # AFR !-> Sur
  m3  = c( Brain ~ Body   , 
           AFR ~ Brain+  Body    ,
           Surv_AFR ~  Brain+ Body 
  ),
  m4  = c( Brain ~ Body   , 
             AFR ~ Brain     ,
             Surv_AFR ~  Brain+ Body 
  ),
  
  # Brain !-> Sur
  m5  = c( Brain ~ Body   , 
             AFR ~ Brain+  Body    ,
             Surv_AFR ~  AFR+ Body 
  ),
  m6  = c( Brain ~ Body   , 
             AFR ~ Brain    ,
             Surv_AFR ~  AFR+ Body 
  ),
  

    ##2##
  # Body -> Brain , AFR -> Brain
  # Body !-> Sur
  m7 = c( Brain ~ Body +AFR  , 
            AFR ~   Body    ,
             Surv_AFR ~  Brain+ AFR 
  ),
  m8  = c( Brain ~ Body +AFR  , 
          
           Surv_AFR ~  Brain+ AFR 
  ),
  # AFR !-> Sur
  m9  = c( Brain ~ Body +AFR  , 
           AFR ~   Body    ,
           Surv_AFR ~  Brain+ Body 
  ),
  
  m10  = c( Brain ~ Body +AFR  , 
             Surv_AFR ~  Brain+ Body 
  ),
  
  # Brain !-> Sur
  m11  = c( Brain ~ Body +AFR  , 
             AFR ~   Body    ,
             Surv_AFR ~  AFR+ Body 
  ),
  m12  = c( Brain ~ Body +AFR  , 
            
             Surv_AFR ~  AFR+ Body 
  ),
  ##4##
  # Body -> Brain , Brain -> AFR
  # Body !-> Sur
  m13  = c( Brain ~ Body +AFR  , 
             Body ~ AFR,
             Surv_AFR ~  Brain+ AFR 
  ),
  
  # AFR !-> Sur
  m14 = c( Brain ~ Body +AFR  , 
           Body ~ AFR,
           Surv_AFR ~  Brain  +Body
  ),

  # Brain !-> Sur
  m15  = c( Brain ~ Body +AFR  , 
             Body ~ AFR,
             Surv_AFR ~   AFR +Body
  )
  
)



position_map <- data.frame(
  term = c("Brain", "Body","AFR","Surv_AFR" ),
  pos_y = c(2, 1, 3,2),
  pos_x = c(1, 2, 2, 3)
)


# 打开 PDF 设备
library(ggplot2)
library(dplyr)
library(patchwork)  # 用于图像拼接

fixed_shrink <- 0.5

# 提取所有模型共有的边
extract_edges <- function(model) {
  edge_mat <- which(model == 1, arr.ind = TRUE)
  from <- rownames(model)[edge_mat[, 1]]
  to <- colnames(model)[edge_mat[, 2]]
  paste(from, to, sep = "->")
}
all_edges_list <- lapply(m_PriSuv, extract_edges)
shared_edges <- Reduce(intersect, all_edges_list)

# 存储所有 ggplot 图
plot_list <- list()

for (i in 1:length(m_PriSuv)) {
  adj_mat <- m_PriSuv[[i]]
  
  edge_df <- which(adj_mat == 1, arr.ind = TRUE)
  edge_df <- data.frame(
    from = rownames(adj_mat)[edge_df[, 1]],
    to   = colnames(adj_mat)[edge_df[, 2]]
  )
  
  edge_df$estimate <- 2
  edge_df <- edge_df %>%
    mutate(
      path = paste(from, to, sep = "->"),
      color = ifelse(path %in% shared_edges, "black", "#4F8095")
    )
  
  plot_data1 <- edge_df %>%
    left_join(position_map, by = c("from" = "term")) %>%
    dplyr::rename(from_x = pos_x, from_y = pos_y) %>%
    left_join(position_map, by = c("to" = "term")) %>%
    dplyr::rename(to_x = pos_x, to_y = pos_y) %>%
    mutate(
      dx = to_x - from_x,
      dy = to_y - from_y,
      norm = sqrt(dx^2 + dy^2),
      ux = dx / norm,
      uy = dy / norm,
      from_x_adj = from_x + (fixed_shrink / 2) * ux,
      from_y_adj = from_y + (fixed_shrink / 2) * uy,
      to_x_adj   = to_x   - (fixed_shrink / 2) * ux,
      to_y_adj   = to_y   - (fixed_shrink / 2) * uy,
      linewidth = abs(estimate) * 6
    )
  
  p <- ggplot(plot_data1) +
    geom_segment(aes(x = from_x_adj, y = from_y_adj, 
                     xend = to_x_adj, yend = to_y_adj, 
                     color = color),
                 lineend = "round") +
    geom_segment(aes(x = to_x_adj - 0.5 * dx, 
                     y = to_y_adj - 0.5 * dy, 
                     xend = to_x_adj, 
                     yend = to_y_adj,
                     color = color),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 inherit.aes = FALSE,
                 lineend = "round") +
    geom_text(data = position_map, 
              aes(x = pos_x, y = pos_y, label = term),
              fontface = "bold", size = 3,
              check_overlap = FALSE) +
    ggtitle(names(m_PriSuv)[i] ) +
    scale_color_identity() +
    
    theme_void() +
    coord_equal( xlim = c(min(position_map$pos_x) - 0.5, max(position_map$pos_x) + 0.5),
                 ylim = c(min(position_map$pos_y) - 0.5, max(position_map$pos_y) + 0.5),
                 expand = FALSE) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 20, b = 10, l = 20) 
    )
  
  plot_list[[i]] <- p
}

# 拼接为 4 列 × 3 行
group_title<-c("Body -> Brain;  Brain + Body -> AFR",
               "Body -> Brain;          Body -> AFR",
               "Body -> Brain;         Brain -> AFR",
               "Body -> Brain;         Brain -> AFR",
               )

combined_plot <- wrap_plots(plot_list, ncol = 6)+  

  plot_annotation(
    title = group_title,
    theme = list(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  )+
  
  theme(
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20)  # 给每张子图更多上边距
  )

combined_plot
# 保存为一页 PDF
PathLocation <- "/Users/zsong/Documents/MPI 2021-2024/Primate Survival Rate/MS/AJP/R1/Pri_Sur_Path/"
ggsave(paste0(PathLocation, "all_models_plot.pdf") , combined_plot, width = 20, height = 10,
       limitsize = FALSE)

#

#Primate Path
library(car)
Pri_Sur_AFR_Path <-Pri_Sur_AFR[,c("tip_label","WT.ADF..kg.", "Br.ADF..g.", "AFR..yr.", 
                                  "Residual_BS","Surv..AFR", "Surv.1st.yr", "Residual_BS") ]
Pri_Sur_AFR_Path$Body<-scale( log10(Pri_Sur_AFR_Path$WT.ADF..kg.), center = T)
Pri_Sur_AFR_Path$Brain<-scale( log10(Pri_Sur_AFR_Path$Br.ADF..g.), center = T)
Pri_Sur_AFR_Path$AFR<-scale( log10(Pri_Sur_AFR_Path$AFR..yr.), center = T)
#Pri_Sur_AFR_Path$Brain<-scale( Pri_Sur_AFR_Path$Residual_BS, center = T)
Pri_Sur_AFR_Path$Surv_AFR<-scale( logit(Pri_Sur_AFR_Path$Surv..AFR), center = T)
#Pri_Sur_AFR_Path$Surv_AFR<-scale(logit( Pri_Sur_AFR_Path$Surv.1st.yr), center = T)

summary(Pri_Sur_AFR_Path)

nrow(Pri_Sur_AFR) # 18
p_Pri<-NA
p_Pri <- phylo_path(m_PriSuv, Pri_Sur_AFR_Path, Pri_ASR_MCCtre); p_Pri
s <- summary(p_Pri); s
#write.csv(s, paste0(PathLocation, "Mammal_Path_model.csv"))
plot(s)

b <- best(p_Pri)
plot(b)
b

avg_mammal <- average(p_Pri, cut_off = 2, avg_method = "conditional")
avg_mammal1<- average(p_Pri,avg_method = "full" )
#plot(avg_mammal, algorithm = 'mds', curvature = 0.1)
coef_plot(avg_mammal1, error_bar = "ci") +
  ggplot2::coord_flip()+
  theme(legend.position = "right",
        text = element_text(size = 15),
        legend.title = element_text(colour = "steelblue", face = "bold.italic", size = 14),
        legend.text = element_text(face = "italic", colour = "steelblue4", size = 10),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

avg_mammal <- average(p_Pri, cut_off = 2, avg_method = "conditional")

coef_plot(avg_mammal, error_bar = "ci") +
  ggplot2::coord_flip()+
  theme(legend.position = "right",
        text = element_text(size = 15),
        legend.title = element_text(colour = "steelblue", face = "bold.italic", size = 14),
        legend.text = element_text(face = "italic", colour = "steelblue4", size = 10),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
