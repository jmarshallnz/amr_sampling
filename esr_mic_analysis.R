library(dplyr)
library(tidyr)
library(stringr)
library(forcats)

fix_ecoli <- function(x) {
  x <- str_replace_all(tolower(x), " ", "")
  ifelse(x == 'notdone', 'nottested', x)
}

fix_entero <- function(x) {
  x <- str_replace_all(tolower(x), " ", "")
  fct_collapse(x, speciated = c("e.faecalis", "e.faecium"))
}

# using the AMR data directly...
bobby <- read.csv("data/esr/Bobby details.csv", skip = 1, stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  select(CPH.No., Sample.Type, Date.Sampled, Quarter, NMD.Code, Campy.PCR, E.coli=Confirmed.E..coli.Result, Entero=Confirmed.Enterococcus.Result) %>%
  mutate(Animal="Calves", Campy.PCR=str_replace_all(Campy.PCR, " ", ""),
         E.coli=fix_ecoli(E.coli),
         Entero=fix_entero(Entero))

porcine <- read.csv("data/esr/Porcine details.csv", skip = 1, stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  select(CPH.No., Sample.Type, Date.Sampled, Quarter, NMD.Code, Campy.PCR, E.coli=Confirmed.E..coli.Result, Entero=Confirmed.Enterococcus.Result) %>%
  rename() %>%
  mutate(Animal="Pigs", Campy.PCR=str_replace_all(Campy.PCR, " ", ""),
         E.coli=fix_ecoli(E.coli),
         Entero=fix_entero(Entero))

poultry <- read.csv("data/esr/Poultry details.csv", skip = 1, stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  select(CPH.No., Sample.Type, Date.Sampled, Quarter, NMD.Code, Campy.PCR = Confirmed..PCR., E.coli=Confirmed.E..coli.Result, Entero=Confirmed.Enterococcus.Result) %>%
  mutate(Animal="Poultry", Campy.PCR=str_replace_all(Campy.PCR, " ", ""),
         E.coli=fix_ecoli(E.coli),
         Entero=fix_entero(Entero))

animals <- bind_rows(bobby, porcine, poultry)

# figure out prevalences. E.coli first
animals %>% filter(Sample.Type == "E. coli Petrifilm") %>%
  group_by(Animal) %>%
  summarize(TotalSamples=n(), Tested=sum(E.coli != 'nottested'), Positive=sum(E.coli == 'e.coli'))

# Campy
animals %>% filter(Sample.Type %in% c("Campy in Bolton", "Campy swab")) %>%
  mutate(Campy.PCR = fct_collapse(Campy.PCR, Campylobacter = c("C.coli", "C.jejuni", "C.jejuni/C.coli"))) %>%
  group_by(Animal) %>%
  summarize(TotalSamples=n(), Tested=sum(Campy.PCR != 'nottested'), Positive=sum(Campy.PCR == 'Campylobacter'))

# Enterococcus
animals %>% filter(Sample.Type == "Entero in MRD") %>%
  group_by(Animal) %>%
  summarize(TotalSamples=n(), Tested=sum(Entero != 'nottested'), Positive=sum(Entero == 'speciated'))

# OK, write out the file for modelling
all <- animals %>% rename(Plant = NMD.Code) %>%
  mutate(Campy.PCR = fct_collapse(Campy.PCR, Campylobacter = c("C.coli", "C.jejuni", "C.jejuni/C.coli")),
         Plant = str_replace_all(Plant, " ", ""))

ec <- all %>% filter(Sample.Type == "E. coli Petrifilm") %>%
  group_by(Plant, Animal, Quarter) %>%
  summarize(Total=n(), Samples=sum(E.coli != 'nottested'), Isolates=sum(E.coli == 'e.coli')) %>%
  mutate(Genus = "E.coli")

# Campy
ca <- all %>% filter(Sample.Type %in% c("Campy in Bolton", "Campy swab")) %>%
  group_by(Plant, Animal, Quarter) %>%
  summarize(Total=n(), Samples=sum(Campy.PCR != 'nottested'), Isolates=sum(Campy.PCR == 'Campylobacter')) %>%
  mutate(Genus = "Campylobacter")

# Enterococcus
en <- all %>% filter(Sample.Type == "Entero in MRD") %>%
  group_by(Plant, Animal, Quarter) %>%
  summarize(Total=n(), Samples=sum(Entero != 'nottested'), Isolates=sum(Entero == 'speciated')) %>%
  mutate(Genus = "Enterococci")

ps_out <- rbind(ec, ca, en)
write.csv(ps_out, "data/plant_prev_for_model.csv", row.names=FALSE)

mb_esc <- read.csv("data/esr/mbrothesc.csv")[,1:19] %>%
  gather(Antimicrobial, MIC, -CPH.No.) %>% mutate(Species='E.coli')
mb_cam <- read.csv("data/esr/mbrothcampy.csv") %>%
  gather(Antimicrobial, MIC, -CPH.No.) %>% mutate(Species='Campylobacter')
mb_ent <- read.csv("data/esr/mbrothent.csv")[1:13] %>%
  gather(Antimicrobial, MIC, -CPH.No.) %>% mutate(Species='Enterococci')

# resistance table
resist <- read.csv("data/esr/resistance.csv")

# bind everything up, and get rid of produce, clean up NMD code
all <- rbind(mb_esc, mb_cam, mb_ent) %>% left_join(animals) %>%
  mutate(Species = ifelse(Species == "Campylobacter", Campy.PCR, Species)) %>%
  left_join(resist) %>% filter(Animal != "Produce") %>%
  mutate(NMD.Code = str_replace_all(NMD.Code, ' ', '')) %>%
  filter(Species %in% c("C.jejuni", "C.coli", "E.coli", "Enterococci"),
         !Antimicrobial %in% c("BAC", "TYL")) %>%
  mutate(Species = fct_collapse(Species, Campylobacter=c("C.jejuni", "C.coli")))

# Decide if resistant or susceptible. This gives the same results as in the final
# report - yay!
desc <- all %>% group_by(Antimicrobial, Animal, Species, NMD.Code) %>%
  summarize(Susceptible = sum(MIC <= Susceptible),
            Resistant = sum(MIC >= Resistant),
            Total = n(),
            Intermediate=Total - Susceptible - Resistant) %>%
  mutate(SusPerc=Susceptible/Total*100,
         IntPerc=Intermediate/Total*100,
         ResPerc=Resistant/Total*100)

# TODO: Can we just regress based on log(MIC) instead of arbitrary class?
#       That allows linear model then instead of logistic regression...


# fit a model to all this shit. We want to allow a random effect per group I guess?
dat_glm <- desc %>% mutate(AMCSpeciesAnimal=interaction(Antimicrobial,Animal,Species)) %>%
  mutate(Susceptible= Susceptible+Intermediate, PlantEffect=interaction(Antimicrobial, Species, NMD.Code)) %>%
  select(Antimicrobial, Animal, Species, AMCSpeciesAnimal, PlantEffect, NMD=NMD.Code, Susceptible, Resistant)

#glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|PlantEffect), family='binomial', data=dat_glm)

fit <- stan_glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|PlantEffect), family='binomial', data=dat_glm)

out <- data.frame(dat_glm, t(apply(posterior_predict(fit), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  mutate(LI = `2.5%`/(Susceptible+Resistant),
         LC = `25%`/(Susceptible+Resistant),
         M = `50%`/(Susceptible+Resistant),
         UC = `75%`/(Susceptible+Resistant),
         UI = `97.5%`/(Susceptible+Resistant)) %>%
  mutate(Plant = factor(paste(Animal, NMD)))
write.csv(out, "data/fitted_resistance.csv", row.names=FALSE)

out2 <- data.frame(dat_glm, t(apply(posterior_predict(fit, re.form = ~0), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  mutate(LI = `2.5%`/(Susceptible+Resistant),
         LC = `25%`/(Susceptible+Resistant),
         M = `50%`/(Susceptible+Resistant),
         UC = `75%`/(Susceptible+Resistant),
         UI = `97.5%`/(Susceptible+Resistant)) %>%
  mutate(Plant = factor(paste(Animal, NMD)))

# This is for model predictions, rather than posterior sampled predictions
# (i.e. posterior sample prediction includes further randomness from the bernoulli/binomial
# outcome)
pos <- plogis(posterior_linpred(fit, re.form=~0))
out <- data.frame(dat_glm,t(apply(pos, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  rename(LI = `2.5%`,
         LC = `25%`,
         M = `50%`,
         UC = `75%`,
         UI = `97.5%`) %>%
  mutate(Plant = factor(paste(Animal, NMD))) %>%
  select(Antimicrobial, Species, Animal, LI, LC, M, UC, UI) %>%
  unique()

write.csv(out, "data/fitted_resistance_FE.csv", row.names=FALSE)
ggplot(out) +
  geom_segment(aes(x=Antimicrobial, xend=Antimicrobial, y=LI, yend=UI)) +
  geom_segment(aes(x=Antimicrobial, xend=Antimicrobial, y=LC, yend=UC), size=2) +
  geom_point(aes(x=Antimicrobial, y=M), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format(), limits=c(0,1)) +
  scale_x_discrete("Antimicrobial") +
  coord_flip() +
  facet_grid(Animal~Species, scales = 'free_y')



resist <- read.csv("data/fitted_resistance.csv") %>%
  left_join(read.csv("data/esr/resistance.csv") %>% select(Antimicrobial, Name))

ggplot(resist %>% filter(Species == "Campylobacter")) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=1) +
  geom_point(aes(x=Plant, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Campylobacter Resistance") +
  scale_y_continuous(labels=scales::percent_format(), limits=c(0,1)) +
  scale_x_discrete("Plant", breaks=NULL) +
  facet_wrap(~Name, scales = 'free_x')

# OK, what we want to do is some sort of an estimate of
# sample size required to detect a (say) 5% change in rate of
# one or more of the things. So we could maybe fit a glmer with
# that information in place and use simr to work out the power?

ec <- dat_glm %>% filter(Species == "E.coli")
dim(ec)

names(ec)
m <- glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|PlantEffect),
      family='binomial', data=ec)

fit <- stan_glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|NMD), family='binomial', data=ec)
fit2 <- stan_glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|PlantEffect), family='binomial', data=ec)

out <- data.frame(ec, t(apply(posterior_predict(fit), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  mutate(LI = `2.5%`/(Susceptible+Resistant),
         LC = `25%`/(Susceptible+Resistant),
         M = `50%`/(Susceptible+Resistant),
         UC = `75%`/(Susceptible+Resistant),
         UI = `97.5%`/(Susceptible+Resistant)) %>%
  mutate(Plant = factor(paste(Animal, NMD)))

ggplot(out) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=1) +
  geom_point(aes(x=Plant, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format(), limits=c(0,1)) +
  facet_wrap(~Antimicrobial, scales = 'free_x')

pos <- plogis(posterior_linpred(fit2))
out <- data.frame(ec,t(apply(pos, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  rename(LI = `2.5%`,
         LC = `25%`,
         M = `50%`,
         UC = `75%`,
         UI = `97.5%`) %>%
  mutate(Plant = factor(paste(Animal, NMD)))

ggplot(out) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=1) +
  geom_point(aes(x=Plant, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format(), limits=c(0,1)) +
  facet_wrap(~Antimicrobial, scales = 'free_x')

# try bootMer?
fit2 <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial', data=ec)
min(predict(fit2, type='response'))
out <- data.frame(ec, M=predict(fit2, type='response')) %>%
  mutate(Plant = factor(paste(Animal, NMD)))

my_pred <- function(x) {
  predict(x, type='response')
}
b <- bootMer(fit2, my_pred, nsim=100)

out <- data.frame(ec, M=predict(fit2, type='response'),
                  t(apply(b$t, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  rename(LI = `2.5%`,
         LC = `25%`,
         M2 = `50%`,
         UC = `75%`,
         UI = `97.5%`) %>%
  mutate(Plant = factor(paste(Animal, NMD)))
ggplot(out) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=1) +
  geom_point(aes(x=Plant, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format()) +
  facet_wrap(~Antimicrobial, scales = 'free_x')

library(merTools)
p1 <- predictInterval(fit2, level=c(0.5), type='probability')
p2 <- predictInterval(fit2, level=c(0.95), type='probability')

out <- data.frame(ec, M=predict(fit2, type='response'),
           p1 %>% rename(LC=lwr, UC=upr),
           p2 %>% rename(LI=lwr, UI=upr)) %>%
  mutate(Plant = factor(paste(Animal, NMD)))
ggplot(out) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=1) +
  geom_point(aes(x=Plant, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format()) +
  facet_wrap(~Antimicrobial, scales = 'free_x')


# fuck, this is quite a tricky situation, no? I mean we can do simple
# stuff just by using a simple 'rule' for power calcs for a single thing
# but we have 300 odd combinations here to test.

# what is the best thing to use in this case?
# I guess the PRIMARY thing we're interested in is whether or not
# the Human-important stuff is now important (it wasn't last time)
# If we can set up a data generating procedure that uses the noise
# in the pilot data along with an actual effect size?
# trouble is we're playing entirely with random effects on the logit
# scale at that...

# I guess we could be interested in detecting k genes being different
# but what pattern of those genes?

# We could try just having one larger maybe?? But how to simulate the
# data? I guess this is fucking hard... Maybe just go for Bonferroni??
d_ec <- ec %>% filter(Antimicrobial == "AMP") %>% ungroup %>% dplyr::select(Animal, NMD, Resistant, Susceptible)

summary(glmer(cbind(Resistant, Susceptible) ~ (1|Animal) + (1 | NMD), family='binomial', data=d_ec))

# Some sample size shit...
alpha=0.05
power=0.80
pe=0.05
pc=0.01
m=12
ICC=0.02
AR=1

# rearrange for power instead...
n <- ((qnorm(1 - alpha/2) + qnorm(power))^2 *
          (pe * (1 - pe) + pc * (1 - pc)) * (1 + (m - 1) * ICC))/(m * (pe - pc)^2)
nTemp <- 0
while (abs(n - nTemp) > 1) {
  nTemp <- n
  n <- ((qt((1 - alpha/2), df = (2 * (nTemp - 
                                          1))) + qt(power, df = (2 * (nTemp - 1))))^2 * 
            (pe * (1 - pe) + pc * (1 - pc)) * (1 + (m - 
                                                      1) * ICC))/(m * (pe - pc)^2)
}
AR = 300 * 2/m / n - 1
nE = ceiling((1/2) * n * (1 + (1/AR))* m) 
nC = ceiling((1/2) * n * (1 + AR) * m)
c(nE, nC)

# So, sample size calculation for comparison of mean etc.
library(sjstats)
library(dplyr)

library(lme4)
dt <- expand.grid(Animal = unique(dat_glm$Animal), Species = unique(dat_glm$Species))
for (i in 1:nrow(dt)) {
  md <- dat_glm %>% filter(Species == dt$Species[i], Animal == dt$Animal[i])
  fit <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial',
               data=md)
  dt$ICC[i] <- icc(fit)["NMD"]
  dt$CV[i] <- md %>% filter(Antimicrobial=="CIP") %>%
    group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
    pull(s) %>% cv
  dt$MCS[i] <- md %>% filter(Antimicrobial=="CIP") %>%
    group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
    pull(s) %>% mean
}

# TODO: Donald wants some subsampling tests. Idea is to subsample the plants,
#       Then increase the samples per plant, so total samples is the same
#       Then recompute CV/MCS, assuming the same ICC.
#       And see what that does to the design effect.
#       That should give some indication of what happens to power.
write.csv(dt, "data/sample_summary.csv", row.names=FALSE)

set.seed(3)
dt_s <- expand.grid(Animal = unique(dat_glm$Animal), Species = unique(dat_glm$Species), Iteration=1:100)
subsample_size <- 0.5
for (i in 1:nrow(dt_s)) {
  md <- dat_glm %>% filter(Species == dt_s$Species[i], Animal == dt_s$Animal[i])
  
  # subsample the plants - we want half of them
  total_samples <- md %>% filter(Antimicrobial == "CIP") %>% summarize(s=sum(Susceptible + Resistant)) %>% pull(s)
  plants <- md %>% filter(Antimicrobial == "CIP") %>% pull(NMD)
  
  # repeat this a bunch of times
  wch_plants <- sample(plants, round(length(plants)*subsample_size))
  subplants <- md %>% filter(NMD %in% wch_plants)
  total_sub <- subplants %>% filter(Antimicrobial == "CIP") %>% ungroup %>% mutate(Total=Susceptible + Resistant) %>% select(NMD, Total)
  # reallocate based on total_samples/total_sub to increase n...
  # first increase the totals
  total_sub$Total <- rmultinom(1, total_samples, total_sub$Total)
  # and then allocate resistant/susceptible
  prev_sub <- subplants %>% mutate(Prev=Resistant / (Susceptible + Resistant)) %>%
    left_join(total_sub, by="NMD") %>% mutate(Resistant = rbinom(n(), Total, Prev), Susceptible = Total - Resistant)

  if (1) {
    dt_s$ICC[i] <- dt %>% filter(Species == dt_s$Species[i], Animal == dt_s$Animal[i]) %>% pull(ICC)
  } else {
  tryCatch({
    fit <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial',
                 data=prev_sub)
    dt_s$ICC[i] <- icc(fit)["NMD"]
  },
  error = function(e) { cat("caught err\n"); dt_s$ICC[i] <- NA},
  warning = function(e) { cat("caught warn\n"); dt_s$ICC[i] <- NA}
  )
  }

  dt_s$CV[i] <- prev_sub %>% filter(Antimicrobial=="CIP") %>%
    group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
    pull(s) %>% cv
  dt_s$MCS[i] <- prev_sub %>% filter(Antimicrobial=="CIP") %>%
    group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
    pull(s) %>% mean
  cat("Done", i, "of", nrow(dt_s), "\n")
}

# now compute the design effect in each case
dt_s %>% mutate(Deff = ((CV^2+1)*MCS-1)*ICC + 1) %>% group_by(Species, Animal) %>% summarize(Deff_L = quantile(Deff, 0.1),
                                                                                           Deff_M = quantile(Deff, 0.5),
                                                                                           Deff_U = quantile(Deff, 0.9)) %>%
  left_join(dt %>% mutate(Deff = ((CV^2+1)*MCS-1)*ICC + 1) %>% select(Species, Animal, Deff), by=c("Species", "Animal"))

write.csv(dt_s, "data/design_effect_subsample.csv", row.names=FALSE)

# Some test computations
dt %>% mutate(MCS2 = MCS*2, Deff = ((CV^2+1)*MCS-1)*ICC + 1, Deff2 = ((CV^2+1)*MCS2-1)*ICC + 1, Incr=Deff2/Deff)

dt_s %>% mutate(Deff = ((CV^2+1)*MCS-1)*ICC + 1) %>% group_by(Species, Animal) %>% summarize(Deff_L = quantile(Deff, 0.1),
                                                                                             Deff_M = quantile(Deff, 0.5),
                                                                                             Deff_U = quantile(Deff, 0.9)) %>%
  left_join(dt %>% mutate(MCS2 = MCS*2, CV2=CV/2, Deff = ((CV^2+1)*MCS-1)*ICC + 1, Deff2 = ((CV2^2+1)*MCS2-1)*ICC + 1, Incr=Deff2/Deff) %>% select(Species, Animal, Deff, Deff2), by=c("Species", "Animal"))

# Hmm, better idea: Reduce dataset one plant at a time for each one. Recompute
# MCS and CV without altering sample size. Note down N. Then compute Deff/N as
# that is the effect on power.

dts <- expand.grid(Animal = unique(dat_glm$Animal), Species = unique(dat_glm$Species))
out <- list()
for (i in 1:nrow(dts)) {
  md <- dat_glm %>% filter(Species == dts$Species[i], Animal == dts$Animal[i])
  fit <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial',
               data=md)
  if (0) {
    dts$ICC[i] <- icc(fit)["NMD"]
  } else {
    dts$ICC[i] <- 0.02
  }
  
  # OK, now go through and drop plants one by one
  plants <- md %>% filter(Antimicrobial=="CIP") %>% mutate(Total = Resistant+Susceptible) %>% arrange(desc(Total))

  for (j in 2:nrow(plants)) {
    CV = cv(plants$Total[1:j])
    MCS = mean(plants$Total[1:j])
    N = sum(plants$Total[1:j])
    out[[length(out)+1]] <- data.frame(Animal = dts$Animal[i], 
                                       Species = dts$Species[i],
                                       ICC = dts$ICC[i],
                                       CV = CV,
                                       MCS = MCS,
                                       N = N,
                                       P = j)
  }
}
d_out <- do.call(rbind, out) %>% mutate(Deff = ((CV^2+1)*MCS-1)*ICC + 1, N_on_Deff = N/Deff) %>%
  mutate(Deff_on_N_I = (CV^2+1)*ICC/P, Deff_on_N_N = (1-ICC)/300, Deff_on_NT = Deff_on_N_I + Deff_on_N_N)


d_out %>% left_join(
  d_out %>% group_by(Animal, Species) %>% top_n(1, N) %>% select(Animal, Species, NT=N, DeffT=Deff,NT_on_DeffT=N_on_Deff)
) %>% mutate(ExtraN = NT_on_DeffT / N_on_Deff * N/NT * 300)

# NOPE - still not right. How about:
# Deff/N needs to be maintained. This is
# ((sd(N_i)^2/mean(N_i)^2 + 1) * mean(N_i) - 1)*ICC + 1) / N
# = (sd(N_i)^2/(N/I) + N/I - 1)*ICC + 1) / N
# = (sd(N_i)^2*I/N^2 + 1/I - 1/N)*ICC + 1/N)
# Deff ~ A * MCS - ICC + 1
# So Deff/N ~ A(I) * 1/I - ICC/N + 1/N
df_out <- d_out %>% left_join(
  d_out %>% group_by(Animal, Species) %>% top_n(1, N) %>% select(Animal, Species, Deff_on_NTT = Deff_on_NT)
) %>% mutate(Deff_on_NNr = pmax(0, Deff_on_NTT - Deff_on_N_I), Nreqd = (1-ICC)/Deff_on_NNr)

library(ggplot2)
ggplot(df_out %>% filter(P >= 5)) + geom_line(aes(x=P, y=Nreqd, col=Species)) + facet_wrap(~Animal, scales = 'free_x') +
  scale_y_continuous(limits=c(0,1500))

ggplot(df_out %>% filter(Animal=="Calves", Species=="E.coli", P >= 5)) + geom_line(aes(x=P, y=Nreqd)) +
  xlab("Number of plants") + ylab("Number of isolates") + scale_y_continuous(limits=c(0,1300), expand = c(0,1))

# TODO: Could repeat with, say ICC = 0.02?

# Sample locations

In the case of a clustered sampling design, the number of samples in total, and the
number of clusters determines the power of the design. In the case of estimating a mean or prevalence,
the variance of the estimator (which specifies detectable differences or lengths of confidence intervals)
takes the form

var(p) prop (CV_I^2 + 1)*ICC/I + (1-ICC)/N

where CV(I) is the coefficient of variation of cluster sizes, ICC is the intracluster correlation, 
I the number of clusters, and N the number of samples. When the ICC is negligible, the second term dominates
which is equivalent to a non-clustered design, and as ICC rises, the former dominates. Thus, in
general, sampling more clusters with fewer samples per cluster generally gives more power in the
case where intracluster correlation is non-trivial. Reducing the number of clusters increases
the first term, requiring the second term to reduce (by increasing N) to compensate in order
to maintain power. To illustrate this effect for this project, we consider the case with the most
plants: E. coli on very young calves, sequentially removing the plant with the least throughput
(and thus least samples) and compute the corresponding number of samples n in order to maintain
power. The figure below shows the effect. Removing the very low throughput plants (plant ##-24) has
little effect on the number of isolates, but as more are removed, the number of isolates increases.
Over 400 isolates are required if only 10 plants are sampled, increasing to 1250 should only 5
plants be sampled. In general the recommendation is to sample as many plants as feasible in
order to obtain samples, with the exception of those that would have few isolates taken due
to very low throughput.


# TODO: It seems poultry has a high ICC?? Need to look into this further,
# and adapt the design effect as needed (way fewer clusters there that
# are larger, right?)
md <- dat_glm %>% filter(Species == "Campylobacter", Animal == "Poultry")
fit <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial',
             data=md)

#cv=0.8
#mean=33

# both pigs + calves do not.
fit2 <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial', data=dat_glm %>% filter(Species == "C.coli", Animal == "Poultry"))
icc(fit2)

# working out cv/mean cluster size
dat_glm %>% filter(Species == "Campylobacter", Animal == "Poultry", Antimicrobial=="CIP") %>%
  group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
  pull(s) %>% cv
dat_glm %>% filter(Species == "Campylobacter", Animal == "Poultry", Antimicrobial=="CIP") %>%
  group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
  pull(s) %>% mean

# using same data as we had normally. This assumes equal sample sizes...
crtpwr.2prop(alpha=0.05, power=NA, m=24, n=12.5, icc=0.02, cv=0.56, p1=0.45, p2=0.55)

crtpwr.2prop(alpha=0.05, power=NA, m=24, n=12.5, icc=0.02, cv=0.56, p1=0.01, p2=0.06)
