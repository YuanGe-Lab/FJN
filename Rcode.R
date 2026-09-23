setwd("work space")
#1.Calculation of the MF index for all bacterial ASVs and its relationship with plant diversity based on ST group
library(readxl)
library(dplyr)
library(broom)
bt<-read_excel("dataset.xlsx", sheet = "bt")
bt[, 2:433] <-bt[, 2:433] + 0.1
bt<-as.data.frame(bt)
rownames(bt)<-bt[,1] 
bt<-bt[,-1] 
tbt<-as.data.frame(t(bt))
result_bt <- data.frame(matrix(ncol = 29223, nrow = 0))
#2.Calculation of the MF index for all bacterial ASVs
for(i in 1:216){
  row1 <- i
  row2 <- i + 216
  result <- log10(tbt[row1,]/tbt[row2,])
  result_bt <- rbind(result_bt,result)
}
result_bt<-as.data.frame(t(result_bt))
write.csv(result_bt,file="result_bt.csv")
tbt_psf<-read_excel("dataset.xlsx", sheet = "bt_psf")
tbt_psf<-as.data.frame(tbt_psf)
rownames(tbt_psf)<-tbt_psf[,1] 
tbt_psf<-tbt_psf[,-1] 
bt_psf<-as.data.frame(t(tbt_psf))
#3.Calculation of the MF for all bacterial ASVs and its relationship with plant diversity
result_bt_psf<-c()
for(i in 2: 29223){
  cor2<-cor.test(~as.numeric(ft_psf[,i])+Richness,method ="spearman",data=ft_psf)#spearman pearson
  result_ft_psf<-rbind(result_ft_psf,
                       c(colnames(ft_psf)[i],cor2$p.value,cor2$estimate))
}
write.csv(result_ft_psf, 'result_ft_psf_fit.csv')
#4.Calculation of the relationship between biomass, PSFs of plant diversity
bio<-read_excel("dataset.xlsx", sheet = "bio")
qqPlot(lm(biomass~group, data =bio), simulate = TRUE, main = 'QQ Plot', labels = FALSE)
bartlett.test(biomass~group, data =bio)
tbio<-filter(bio,group == "U"|group == "T")
pbio<-filter(bio,group == "U"|group == "P")
#treat fit
tbio$group<-factor(bio$group, levels =c("U", "T"))
p2.1<-ggplot(bio, aes(richness, biomass,linetype=group)) + 
  geom_point(aes(color = group),size = 2,position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(aes(color = group,fill = group),alpha=0.2,size = 1,method = "lm",formula = 'y ~ x', level=0.95,se = T)+
  stat_cor(aes(color = group),method="pearson")+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  scale_color_manual(values = c('#359023', '#ed7d31')) +
  scale_fill_manual(values = c('#359023', '#ed7d31'))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_x_continuous(limits=c(0.5,8.5), breaks=seq(1,8,1))+
  scale_y_continuous(limits=c(-300,2500), breaks=seq(0,2000,1000))+
  theme(legend.position = "none")
#dominants fit
pbio$group<-factor(pbio$group, levels =c("U", "P"))
p2.2<-ggplot(pbio, aes(richness, biomass,linetype=group)) + 
  geom_point(aes(color = group),size = 2,position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(aes(color = group,fill = group),alpha=0.2,size = 1,method = "lm",formula = 'y ~ x', level=0.95,se = T)+
  stat_cor(aes(color = group),method="pearson")+
  scale_linetype_manual(values = c('dashed', 'solid'))+
  scale_color_manual(values = c('#359023', '#0072b2')) +
  scale_fill_manual(values = c('#359023', '#0072b2'))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_x_continuous(limits=c(0.5,8.5), breaks=seq(1,8,1))+
  scale_y_continuous(limits=c(-300,2500), breaks=seq(0,2000,1000))+
  theme(legend.position = "none")
s2<-read_excel("dataset.xlsx", sheet = "psf2")
qqPlot(lm(psf~richness, data = s2), simulate = TRUE, main = 'QQ Plot', labels = FALSE)
ts2<-filter(s2,group == "UT")
ps2<-filter(s2,group == "UP")
#fit
p2.3<-ggplot(ts2, aes(richness, psf))+ 
  geom_point(size = 2,position = position_jitter(width = 0.1, height = 0.1),color = "#ed7d31") +
  geom_hline(yintercept = 0, linetype=2)+
  stat_smooth(method = "lm",formula = y ~ x, se = T,col="#ed7d31",fill="#ed7d31",alpha=0.15)+
  stat_cor(method="pearson")+
  #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label..,sep = '~~~~')), formula = y ~ x, parse = T)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_x_continuous(limits=c(0.5,8.5), breaks=seq(1,8,1))+
  scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2,2,1))+
  theme(legend.position = "top")
p2.4<-ggplot(ps2, aes(richness, psf))+ 
  geom_point(size = 2,position = position_jitter(width = 0.1, height = 0.1),color = "#0072b2") +
  geom_hline(yintercept = 0, linetype=2)+
  stat_smooth(method = "lm",formula = y ~ x, se = T,col="#0072b2",fill="#0072b2",alpha=0.15)+
  stat_cor(method="pearson")+
  #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label..,sep = '~~~~')), formula = y ~ x, parse = T)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_x_continuous(limits=c(0.5,8.5), breaks=seq(1,8,1))+
  scale_y_continuous(limits=c(-2.5,2.5), breaks=seq(-2,2,1))+
  theme(legend.position = "top")
library(cowplot)
p2 <- plot_grid(p2.1,p2.3,p2.2,p2.4,ncol = 2, nrow=2)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggdraw(p2)
ggsave('Fig2.pdf', p2, width = 8, height = 7)
#5.effect：complementarity and selection effect##############
eff<-read_excel("data.xlsx", sheet = "effect")
qqPlot(lm(CE~group, data = eff), simulate = TRUE, main = 'QQ Plot', labels = FALSE)
p2.5<-
  ggplot(data = eff, aes(x = richness, y = CE, linetype=group)) +
  geom_boxplot(aes(color = group,group=group_richness),size=0.8)+
  geom_point(size = 1, shape=1,aes(color = group,group=group_richness),position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(aes(color = group,fill = group),alpha=0.1,span = 0.2, size = 1,method = "lm",formula = 'y ~ x', level=0.95,se = T)+
  stat_cor(aes(color = group),method="pearson")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_poly_eq(aes(label = paste(..p.value.label..)), formula = y ~ x, parse = T,label.x.npc = 0.03, label.y.npc = 0.95) + 
  scale_linetype_manual(values = c('solid', 'solid', 'solid')) +
  scale_color_manual(values = c("#359023","#ed7d31","#0072b2"))+
  scale_fill_manual(values = c("#359023","#ed7d31","#0072b2"))+
  scale_y_continuous(limits=c(-1500,6000), breaks=seq(0,6000,2000))+
  theme(legend.position = "none")+
  facet_wrap(~group)
qqPlot(lm(CE~group, data = eff), simulate = TRUE, main = 'QQ Plot', labels = FALSE)
p2.6<-
  ggplot(data = eff, aes(x = richness, y = SE, linetype=group)) +
  geom_boxplot(aes(color = group,group=group_richness),size=0.8)+
  geom_point(size = 1, shape=1,aes(color = group,group=group_richness),position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(aes(color = group,fill = group),alpha=0.1,span = 0.2, size = 1,method = "lm",formula = 'y ~ x', level=0.95,se = T)+
  stat_cor(aes(color = group),method="pearson")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_poly_eq(aes(label = paste(..p.value.label..)), formula = y ~ x, parse = T,label.x.npc = 0.03, label.y.npc = 0.95) + 
  scale_linetype_manual(values = c('solid', 'solid', 'solid')) +
  scale_color_manual(values = c("#359023","#ed7d31","#0072b2"))+
  scale_fill_manual(values = c("#359023","#ed7d31","#0072b2"))+
  scale_y_continuous(limits=c(-6000,3000), breaks=seq(-6000,0,2000))+
  theme(legend.position = "none")+
  facet_wrap(~group)
library(cowplot)
p2_0 <- plot_grid(p2.5,p2.6,ncol = 2, nrow=1)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggdraw(p2_0)
ggsave('Fig2_0.pdf', p2_0, width = 8, height = 4)
