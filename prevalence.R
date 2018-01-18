library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(rstanarm)

# read in the prevalence data and tidy
total_samples <- read.csv("data/tabula-Heffernan_table01.csv", na.strings=c("NA", "NA3"))[,-1]
names(total_samples)[-1] <- c("Campylobacter", "Salmonella", "E. coli", "Enterococci")
total_samples$Animal <- c(rep("Very young calves", 5), rep("Pigs", 4), rep("Poultry", 4), rep("Fresh produce", 2))
total_samples <- total_samples[-c(1,4,6,10,14,15),]
total_samples$Sample_Type = rep(c("Sample", "Isolate", "Selected"), 3)
total_samples <- total_samples[,-1]
total_samples <- total_samples %>% gather(Genus, Count, -Animal, -Sample_Type) %>%
  spread(Sample_Type, Count) %>%
  mutate(Prevalence = Isolate/Sample * 100)

# read in the prevalence per plant and tidy
plant08 <- read.csv("data/tabula-Heffernan_table08.csv")[,c(1:10)]; plant08$Genus="Campylobacter"
plant09 <- read.csv("data/tabula-Heffernan_table09.csv")[,c(1:10)]; plant09$Genus="E.coli"
plant10 <- read.csv("data/tabula-Heffernan_table10.csv")[,c(1:10)]; plant10$Genus="Enterococci"

plant_samples <- rbind(plant08, plant09, plant10) %>%
  gather(Quarter, Value, `Quarter.1`:Total) %>%
  extract(Quarter, into="Quarter", regex='([0-9]+)') %>%
  replace_na(list(Quarter="Total")) %>%
  extract(Value, into=c("Samples", "Isolates", "Selected"), regex=c("([0-9]+) / ([0-9]+) ([0-9]+)"), convert=TRUE) %>%
  mutate(Prevalence=Isolates/Samples * 100)

# check the totals are right (sanity...)
plant_samples %>% group_by(Plant, Genus) %>% summarize(Total = sum(Samples[Quarter != "Total"], na.rm=TRUE), Total2 = unique(Samples[Quarter == "Total"])) %>% filter(Total != Total2)
plant_samples %>% group_by(Plant, Genus) %>% summarize(Total = sum(Isolates[Quarter != "Total"], na.rm=TRUE), Total2 = unique(Isolates[Quarter == "Total"])) %>% filter(Total != Total2)
plant_samples %>% group_by(Plant, Genus) %>% summarize(Total = sum(Selected[Quarter != "Total"], na.rm=TRUE), Total2 = unique(Selected[Quarter == "Total"])) %>% filter(Total != Total2)

# filter out the totals - we don't need them
plant_samples <- plant_samples %>% filter(Quarter != "Total") %>% mutate(Quarter = as.numeric(Quarter))

# work out overall prevalence per plant per genus
prev <- plant_samples %>% group_by(Plant, Genus) %>%
  summarize(Samples = sum(Samples, na.rm=TRUE), Isolates = sum(Isolates, na.rm=TRUE), Prevalence = Isolates/Samples*100) %>%
  left_join(plant_samples %>% select(Plant, Target, Animal, Island) %>% unique)

ggplot(prev, aes(x=Plant, y=Prevalence)) + geom_col() + facet_grid(Genus~Animal, scales = 'free_x')

# fit a GLMER to the data
#mod <- glmer(cbind(Isolates, Samples-Isolates) ~ Animal + (1|Plant), family='binomial', data=subset(prev, Genus == "E.coli"))
#boots <- lme4::bootMer(mod, function(x) { predict(x, re.form=NULL) }, nsim=100, use.u=FALSE, type="parametric")

# try stan!
#ecoli <- prev %>% filter(Genus == "E.coli")
#mod <- glmer(cbind(Isolates, Samples-Isolates) ~ Animal + (1|Plant), family='binomial', data=ecoli)

dat <- prev %>% filter(Genus == "E.coli")
mstan <- stan_glmer(cbind(Isolates, Samples-Isolates) ~ Animal + (1|Plant), family='binomial',
                    data=dat)
stan_out <- data.frame(Plant=1:43, t(apply(posterior_predict(mstan), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  left_join(dat) %>% mutate(LI = `X2.5.`/Samples,
                              LC = `X25.`/Samples,
                              M = `X50.`/Samples,
                              UC = `X75.`/Samples,
                              UI = `X97.5.`/Samples)
lin_out <- plogis(posterior_linpred(mstan))
stan_out <- data.frame(Plant=1:43, t(apply(lin_out, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  left_join(dat) %>% rename(LI = `X2.5.`,
                             LC = `X25.`,
                             M = `X50.`,
                             UC = `X75.`,
                             UI = `X97.5.`)
ggplot(stan_out %>% filter(Samples > 0) %>% mutate(Throughput = ifelse(Target <= 6, "Low", "High"))) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=2) +
  geom_point(aes(x=Plant, y=M, col=Animal, fill=Throughput), size=3, shape=21) +
  guides(col=guide_legend(title=NULL)) +
  theme_bw() + ylab("E.coli prevalence") +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_x_continuous(breaks=NULL) +
  scale_fill_manual(values=c("white", "black"))

dat2 <- prev %>% filter(Genus == "Campylobacter")
mstan2 <- stan_glmer(cbind(Isolates, Samples-Isolates) ~ Animal + (1|Plant), family='binomial',
                     data=dat2)
stan_out2 <- data.frame(Plant=1:43, t(apply(posterior_predict(mstan2), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  left_join(dat2) %>% mutate(LI = `X2.5.`/Samples,
                              LC = `X25.`/Samples,
                              M = `X50.`/Samples,
                              UC = `X75.`/Samples,
                              UI = `X97.5.`/Samples)
lin_out2 <- plogis(posterior_linpred(mstan2))
stan_out2 <- data.frame(Plant=1:43, t(apply(lin_out2, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  left_join(dat2) %>% rename(LI = `X2.5.`,
                             LC = `X25.`,
                             M = `X50.`,
                             UC = `X75.`,
                             UI = `X97.5.`)
ggplot(stan_out2 %>% filter(Samples > 0) %>% mutate(Throughput = ifelse(Target <= 6, "Low", "High"))) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=2) +
  geom_point(aes(x=Plant, y=M, col=Animal, fill=Throughput), size=3, shape=21) +
  guides(col=guide_legend(title=NULL)) +
  theme_bw() + ylab("Campylobacter prevalence") +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_x_continuous(breaks=NULL) +
  scale_fill_manual(values=c("white", "black"))

dat3 <- prev %>% filter(Genus == "Enterococci")
mstan3 <- stan_glmer(cbind(Isolates, Samples-Isolates) ~ Animal + (1|Plant), family='binomial',
                     data=dat3)
stan_out3 <- data.frame(Plant=1:43, t(apply(lin_out3, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  left_join(dat3) %>% mutate(LI = `X2.5.`/Samples,
                             LC = `X25.`/Samples,
                             M = `X50.`/Samples,
                             UC = `X75.`/Samples,
                             UI = `X97.5.`/Samples)
lin_out3 <- plogis(posterior_linpred(mstan3))
stan_out3 <- data.frame(Plant=1:43, t(apply(lin_out3, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  left_join(dat3) %>% rename(LI = `X2.5.`,
                             LC = `X25.`,
                             M = `X50.`,
                             UC = `X75.`,
                             UI = `X97.5.`)
ggplot(stan_out3 %>% filter(Samples > 0) %>% mutate(Throughput = ifelse(Target <= 6, "Low", "High"))) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=2) +
  geom_point(aes(x=Plant, y=M, col=Animal, fill=Throughput), size=3, shape=21) +
  guides(col=guide_legend(title=NULL)) +
  theme_bw() + ylab("Enterococci prevalence") +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_x_continuous(breaks=NULL) +
  scale_fill_manual(values=c("white", "black"))



# try and predict the overall prevalence directly
pred1 <- plogis(
  posterior_linpred(mstan, newdata=data.frame(Animal = levels(dat$Animal)), re.form=~0))
prev1 <- data.frame(Animal=levels(dat$Animal),
                    t(apply(pred1, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  rename(LI = `X2.5.`,
         LC = `X25.`,
         M = `X50.`,
         UC = `X75.`,
         UI = `X97.5.`)
pred2 <- plogis(
  posterior_linpred(mstan2, newdata=data.frame(Animal = levels(dat$Animal)), re.form=~0))
prev2 <- data.frame(Animal=levels(dat$Animal),
                    t(apply(pred2, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  rename(LI = `X2.5.`,
         LC = `X25.`,
         M = `X50.`,
         UC = `X75.`,
         UI = `X97.5.`)
pred3 <- plogis(
  posterior_linpred(mstan3, newdata=data.frame(Animal = levels(dat$Animal)), re.form=~0))
prev3 <- data.frame(Animal=levels(dat$Animal),
                    t(apply(pred3, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  rename(LI = `X2.5.`,
         LC = `X25.`,
         M = `X50.`,
         UC = `X75.`,
         UI = `X97.5.`)

# Table them up for the report
prev2009 <- rbind(prev1 %>% mutate(Genus='E. coli'),
      prev2 %>% mutate(Genus='Campylobacter'),
      prev3 %>% mutate(Genus='Enterococcus')) %>%
  mutate(CI = sprintf("%.1f (%.1f, %.1f)", M*100, LI*100, UI*100)) %>%
  mutate(Samples = signif(300 / LC, 2))

write.csv(prev2009, "data/2009_sample_size_calcs.csv", row.names=FALSE)










# Now look at prevalance of each type among the thingees
all <- bind_rows(stan_out, stan_out2, stan_out3)

write.csv(all, "data/fitted_prevalence.csv", row.names=FALSE)
ggplot(all %>% filter(Samples > 0) %>% mutate(Throughput = ifelse(Target <= 6, "Low", "High"))) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=2) +
  geom_point(aes(x=Plant, y=M, col=Animal, fill=Throughput), size=3, shape=21) +
  facet_wrap(~Genus, ncol=1) +
  guides(col=guide_legend(title=NULL)) +
  theme_bw() + ylab("Prevalence") +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_x_continuous(breaks=NULL) +
  scale_fill_manual(values=c("white", "black"))





d1 <- read.csv("data/tabula-Heffernan_table02.csv", skip=2) %>% gather(Cat, Percent, -X) %>%
  rename(Antimicrobial=X) %>%
  extract(Cat, into=c("Animal", "Isolates"), regex="(.*)\\.n\\.([0-9]+)", convert=TRUE) %>%
  mutate(Resistant = round(Isolates*Percent/100))
d1$Genus = c(rep("C.jejuni", 7*3), rep("C.coli", 7*4), rep("C.lari", 7*1))

d2 <- read.csv("data/tabula-Heffernan_table03.csv", skip=1) %>% gather(Cat, Percent, -X) %>%
  rename(Antimicrobial=X) %>%
  extract(Cat, into=c("Animal", "Isolates"), regex="(.*)\\.n\\.([0-9]+)", convert=TRUE) %>%
  mutate(Resistant = round(Isolates*Percent/100))
d2$Genus <- "E.coli"

d3 <- read.csv("data/tabula-Heffernan_table04.csv", skip=2, na.strings=c('', '-')) %>% gather(Cat, Percent, -X) %>%
  rename(Antimicrobial=X) %>%
  extract(Cat, into=c("Animal", "Isolates"), regex="(.*)\\.n\\.([0-9]+)", convert=TRUE) %>%
  mutate(Resistant = round(Isolates*Percent/100))
d3$Genus = rep(c("E.faecalis", "E.faecium"), each=40)
d3 

# combine the campy's
library(forcats)
resist <- rbind(d1, d2, d3) %>% filter(!Animal %in% c('Fresh.produce', 'Total')) %>%
  filter(!is.na(Resistant)) %>%
  mutate(Genus = fct_collapse(Genus, Campylobacter=c("C.jejuni", "C.coli", "C.lari"))) %>%
  mutate(Final = interaction(Animal, Genus, Antimicrobial)) %>%
  group_by(Final, Animal, Genus, Antimicrobial) %>% summarize(Isolates = sum(Isolates), Resistant = sum(Resistant)) %>%
  ungroup

maxi_mod <- stan_glm(cbind(Resistant, Isolates-Resistant) ~ Final, family='binomial',
         data=resist)

out <- cbind(resist, t(apply(posterior_predict(maxi_mod), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  mutate(LI = `2.5%`/Isolates,
         LC = `25%`/Isolates,
         M = `50%`/Isolates,
         UC = `75%`/Isolates,
         UI = `97.5%`/Isolates)

ggplot(out) +
  geom_segment(aes(x=Antimicrobial, xend=Antimicrobial, y=LI, yend=UI)) +
  geom_segment(aes(x=Antimicrobial, xend=Antimicrobial, y=LC, yend=UC, col=Animal), size=2) +
  geom_point(aes(x=Antimicrobial, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format()) +
  coord_flip() +
  facet_grid(Genus~Animal, scales = 'free_y')

resist <- resist %>% mutate(AntiGenus = interaction(Antimicrobial, Genus))

# try simr for power analysis
fit <- glmer(cbind(Resistant, Isolates-Resistant) ~ AntiGenus + (AntiGenus|Animal), family='binomial', data=resist)

# fit using a glmer - yay! Now use simr for power analysis...
library(simr)
coef(fit)["Final"][1]

# hmm, what scale are they on, and how can we figure out how much it changes by?
# I guess the idea is we can fit a separate set to the same thing, doubling up the information
# and then set the sigma accordingly?
resist_new <- rbind(cbind(resist, Study=1), cbind(resist, Study=2))
resist_new$Resistant[resist_new$Study==2] <- resist_new$Resistant[resist_new$Study == 2] + sample(1:3, nrow(resist), replace=TRUE)
resist_new$Study <- factor(resist_new$Study)
fit <- glmer(cbind(Resistant, Isolates-Resistant) ~ AntiGenus*Study + (AntiGenus*Study|Animal), family='binomial', data=resist_new)


fit <- glmer(cbind(Resistant, Isolates-Resistant) ~ (Study|Final), family='binomial', data=resist_new)
v <- VarCorr(fit)$Final
v[2,2] <- 1
v[1,2] <- v[2,1] <- 0.5
VarCorr(fit)$Final <- v
doTest(fit, compare(.~(1|Final)))

fit2 <- glm(cbind(Resistant, Isolates-Resistant) ~ Study*Final, family='binomial', data=resist_new)

powerSim(fit, compare(.~(1|Final)), nsim=10)

genius <- resist$Genus %>% unique

for (g in genius) {
  
dat <- resist %>% filter(Genus == g, !is.na(Percent))
if (dat %>% pull(Animal) %>% unique %>% length == 1) {
  mod <- stan_glm(cbind(Resistant, Isolates-Resistant) ~ Antimicrobial, family='binomial',
                  data=dat)
} else {
  mod <- stan_glm(cbind(Resistant, Isolates-Resistant) ~ Animal*Antimicrobial, family='binomial',
                  data=dat)
}

out <- cbind(dat, t(apply(posterior_predict(mod), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))) %>%
  mutate(LI = `2.5%`/Isolates,
         LC = `25%`/Isolates,
         M = `50%`/Isolates,
         UC = `75%`/Isolates,
         UI = `97.5%`/Isolates)
ggplot(out) +
  geom_segment(aes(x=Antimicrobial, xend=Antimicrobial, y=LI, yend=UI)) +
  geom_segment(aes(x=Antimicrobial, xend=Antimicrobial, y=LC, yend=UC, col=Animal), size=2) +
  geom_point(aes(x=Antimicrobial, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format()) +
  coord_flip() +
  facet_wrap(~Animal)
}









library(boot)
# and overall prevalence per animal type and genus
plant_samples %>% group_by(Animal, Genus) %>%
  summarize(Samples = sum(Samples, na.rm=TRUE), Isolates = sum(Isolates, na.rm=TRUE), Prevalence = Isolates/Samples*100, Target=sum(Target)/4)

# read in the plant data
plant_data <- read.csv("data/2009_plan.csv", na.strings = c("NA", "")) %>% fill(Genus, Animal)

par(mfrow=c(1,2))
dat <- plant_data %>% filter(Animal == "Pigs", Genus == "Campylobacter")
mod <- lm(Target ~ Throughput, data=dat)
plot(Target ~ Throughput, data=dat)
abline(mod, col='red')

dat <- plant_data %>% filter(Animal == "Very young calves", Genus == "Campylobacter")
mod <- lm(Target ~ Throughput, data=dat)
plot(Target ~ Throughput, data=dat)
abline(mod, col='red')

# read in the NMD data from April 2016 through June 2017
nmd <- read.csv("data/2016_nmd_samples.csv", na.strings=c("NA", "")) %>% fill(Animal)

