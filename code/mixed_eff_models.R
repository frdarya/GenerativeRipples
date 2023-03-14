library(lme4)
library(effects)
library(car)
library(tidyverse)
library(MASS)
library(ggpubr)
library(emmeans)
library(pals)
library(optimx)
library(interactions)

df <- read.csv('Vaz_hpf200_25ms_allrip.csv')
df$channel <- factor(df$Channel)
df$subI <- factor(df$subI)
df <- df %>% mutate(pre_post=ifelse(rip_peak<=500,'pre',ifelse(rip_peak > 500, 'post', 'no_rip')))
ms <- seq(from = -1, to = 1, by = 0.002)
df <- df %>% mutate(rip_peak_ms = ms[rip_peak])

df <- df %>% mutate(Speed=1/(rt/1000))
df$entropy.cut = cut(df$entropy, breaks=c(0.6, 1, 1.2, 1.4, 1.6, 1.8, 2))
df$rip_peak_ms.cut <- cut(df$rip_peak_ms, breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
df_prop <-with(df, table(entropy.cut,rip_peak_ms.cut, useNA='no'))


###### rip duration ######
near_stim <- subset(df, c(rip_peak_ms.cut== '(0,0.2]' | rip_peak_ms.cut== '(-0.2,0]'))
m2 <- glmer(rip_dur ~ entropy*rip_peak_ms.cut + surprise*rip_peak_ms.cut+meanEnt+ (1|subI),
             data = near_stim,family = Gamma(link=log),
             control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(m2)
Anova(m2, 3)
sim_slopes(m2b,entropy,rip_peak_ms.cut)
dur_peak_ent_near <- ggplot(data.frame(effect('entropy*rip_peak_ms.cut',m2b)),
                       aes(entropy,fit,group=rip_peak_ms.cut, fill= rip_peak_ms.cut, 
                           color=  rip_peak_ms.cut, ymin=lower,ymax=upper))+
  geom_line(size = 2, alpha = 20)+geom_ribbon(alpha=0.5)+ylab('Ripple duration (s)')+xlab('Entropy')+
  theme_classic(base_size = 30)+labs(fill='Ripple peak time')+ guides(color="none")+
  scale_y_continuous(limits=c(0.025,0.045)) + scale_x_continuous(limits=c(1,2))
dur_peak_surp_near <- ggplot(data.frame(effect('rip_peak_ms.cut*surprise',m2b)),
                        aes(surprise,fit,group=rip_peak_ms.cut,fill=rip_peak_ms.cut,
                            color = rip_peak_ms.cut,ymin=lower,ymax=upper))+
  geom_line(size = 2, alpha=20)+geom_ribbon(alpha=0.5)+xlab('Surprise')+ylab('Ripple duration (s)')+scale_y_continuous(limits=c(0.025,0.045))+
  theme_classic(base_size = 30,base_family='sans')+labs(group='Ripple peak time')+ guides(color="none")
ggarrange(dur_peak_ent_near, dur_peak_surp_near, ncol = 2, nrow = 1, common.legend = TRUE)


# supp fig 2h
m2a <- glmer(rip_dur ~ entropy*rip_peak_ms.cut + surprise*rip_peak_ms.cut+meanEnt+ (1|subI), data = df, family = Gamma(link=log))
summary(m2a)
Anova(m2a, 3)
dur_peak_ent <- ggplot(data.frame(effect('entropy*rip_peak_ms.cut',m2a)),
                       aes(entropy,fit,group=rip_peak_ms.cut, fill= rip_peak_ms.cut, 
                           color=  rip_peak_ms.cut, ymin=lower,ymax=upper))+
  geom_line(size = 2, alpha = 20)+geom_ribbon(alpha=0.05)+theme_classic()+ylab('Ripple duration')+xlab('Entropy')
dur_peak_surp <- ggplot(data.frame(effect('rip_peak_ms.cut*surprise',m2a)),
                        aes(surprise,fit,group=rip_peak_ms.cut,fill=rip_peak_ms.cut,
                            color = rip_peak_ms.cut,ymin=lower,ymax=upper))+
  geom_line(size = 2, alpha=20)+geom_ribbon(alpha=0.05)+theme_classic()+ylab('Ripple duration')+xlab('Surprise')
ggarrange(dur_peak_ent, dur_peak_surp, ncol = 2, nrow = 1, common.legend = TRUE)



###### 2D hist ########
df2 <- read.csv('norm_count_ent_time.csv')
df2 <- na.omit(df2)

#df$entropy <- scale(df$entropy, center = FALSE, scale = TRUE)
df2 <- df2 %>% mutate(pre_post=ifelse(time<0,'pre',ifelse(time >= 0, 'post', 'no_rip')))
df2$pre_post <- factor(df$pre_post)
df2$abs_time <- abs(df$time)
df2 <- df2 %>% mutate(norm_count=ifelse(norm_count==0,0.01,norm_count))
df2$time_ms <- df$time/1000
df2$subI <- factor(df$subI)
df2$entropy <- factor(df$entropy)
df2$time <- factor(df$time)

df2_rest <- subset(df2, norm_count < 20)
df2_pre <- subset(df2, pre_post=='pre')
df2_post <- subset(df2, pre_post=='post')

m_all<- glmer(norm_count~entropy*time+(1|subI), data=df2, family = Gamma(link=log),
              control=glmerControl(tolPwrss=1e-3, optimizer = 'optimx', optCtrl=list(method='nlminb')))
Anova(m_all,3)
emm1 <- emmeans(m_all, specs = pairwise ~ entropy:time)
df_emm <- as.data.frame(emm1)

###### predict Speed with ripples ######
all_trls <-read.csv('Vaz_hpf200_25ms_alltrls.csv')
all_trls$channel <- factor(all_trls$Channel)
all_trls$subI <- factor(all_trls$subI)
all_trls <- all_trls %>% mutate(rip_status=ifelse( is.na(rip_peak),'no_rip',ifelse(rip_peak<=500,'pre', 'post')))
all_trls <- all_trls %>% mutate(Speed=1/(rt/1000))
ms <- seq(from = -1, to = 1, by = 0.002)
all_trls <- all_trls %>% mutate(rip_peak_ms = ms[rip_peak])
all_trls <- all_trls %>% mutate(rip_start_ms = ms[rip_start])
all_trls <- all_trls %>% mutate(trl_blk = trl_num%%40)
all_trls$trl_blk <- ifelse(all_trls$trl_blk==0,40,all_trls$trl_blk)

all_trls$rip_peak_ms.cut <- cut(all_trls$rip_peak_ms, breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))
all_trls$entropy.cut = cut(all_trls$entropy, breaks=c(0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))


pre <- subset(all_trls, rip_status != 'post')
just_pre <- subset(pre, rip_status == 'pre')
nsubs = levels(unique(all_trls$subI))
unq_trls <- tibble()
for (sub in nsubs) {
  print(sub)
  sub_dat <- tibble()
  d <- subset(pre, as.numeric(pre$subI) == as.numeric(sub))
  pre_idx <- subset(d,d$rip_status == 'pre') 
  pre_trls <- pre_idx[!duplicated(pre_idx$trl_num), ]
  no <- subset(d,d$rip_status == 'no_rip')
  sub_dat <- rbind(no,pre_trls)
  unq_trls<- rbind(unq_trls,sub_dat)
}
sub <- 8
sub_dat <- tibble()
d <- subset(pre, as.numeric(pre$subI) == as.numeric(sub))
pre_idx <- subset(d,d$rip_status == 'pre') 
pre_trls <- pre_idx[!duplicated(pre_idx$trl_num), ]
no <- subset(d,d$rip_status == 'no_rip')
sub_dat <- rbind(no,pre_trls)
unq_trls<- rbind(unq_trls,sub_dat)


red_unq <- subset(unq_trls, c(entropy >= 1.6& entropy <= 1.9))
log_rt_entropy <- lmer(log(rt)~ rip_status*entropy +rip_status*surprise + meanEnt+rip_rate+(1|subI), data=red_unq)
Anova(log_rt_entropy,3)

log_rt_ent<- ggplot(data.frame(effect('rip_status*entropy',log_rt_entropy)),aes(x = entropy, y = fit, group= rip_status, colour = rip_status,
                                                                fill = rip_status,ymin = lower, ymax = upper))+ 
  geom_line(size=2) +
  geom_ribbon(alpha = 0.2)+
  scale_color_manual(labels = c("No ripple", "Pre-stim ripple"), values = c("#1F77B4", "#F8766D")) +
  scale_fill_manual(values = c("#1F77B4", "#F8766D")) +
  theme_classic(base_size = 20,base_family='sans')+labs(color='Ripple status')+ guides(fill="none")+
  xlab('Entropy')+ylab('log(RT)')
log_rt_surp<- ggplot(data.frame(effect('rip_status*surprise',log_rt_entropy)),aes(x = surprise, y = fit, group= rip_status, colour = rip_status,
                                                                                fill = rip_status,ymin = lower, ymax = upper))+ 
  geom_line(size=2) +
  geom_ribbon(alpha = 0.2)+
  theme_classic(base_size=18)+xlab('Surprise')+ylab('log(RT)')
ggarrange(log_rt_ent, log_rt_surp, ncol=2, common.legend = TRUE)
sim_slopes(log_rt_entropy,entropy,rip_status)

# sup
log_rt_all <- lmer(log(rt)~ rip_status*entropy +rip_status*surprise + rip_rate+meanEnt+(1|subI), data=unq_trls)
Anova(log_rt_all,3)
log_rt_all_ent<- ggplot(data.frame(effect('rip_status*entropy',log_rt_all)),aes(x = entropy, y = fit, group= rip_status, colour = rip_status,
                                                                                fill = rip_status,ymin = lower, ymax = upper))+ 
  geom_line(size=2) +
  geom_ribbon(alpha = 0.2)+
  scale_color_manual(labels = c("No ripple", "Pre-stim ripple"), values = c("#1F77B4", "#F8766D")) +
  scale_fill_manual(values = c("#1F77B4", "#F8766D")) +
  theme_classic(base_size = 18,base_family='sans')+labs(color='Ripple status')+ guides(fill="none")+
  xlab('Entropy')+ylab('log(RT)')
log_rt_all_surp<- ggplot(data.frame(effect('rip_status*surprise',log_rt_all)),aes(x = surprise, y = fit, group= rip_status, colour = rip_status,
                                                                                  fill = rip_status,ymin = lower, ymax = upper))+ 
  geom_line(size=2) +
  geom_ribbon(alpha = 0.2)+    scale_color_manual(labels = c("No ripple", "Pre-stim ripple"), values = c("#1F77B4", "#F8766D")) +
  scale_fill_manual(values = c("#1F77B4", "#F8766D")) +
  theme_classic(base_size=18)+xlab('Surprise')+ylab('log(RT)')
ggarrange(log_rt_all_ent, log_rt_all_surp, ncol=2, common.legend = TRUE)



