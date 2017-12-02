library(dplyr)
library(tidyr)
library(stringr)

# using the AMR data directly...
bobby <- read.csv("data/esr/Bobby details.csv", skip = 1, stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  select(CPH.No., Date.Sampled, Quarter, NMD.Code, Campy.PCR) %>%
  mutate(Animal="Calves", Campy.PCR=str_replace_all(Campy.PCR, " ", ""))

porcine <- read.csv("data/esr/Porcine details.csv", skip = 1, stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  select(CPH.No., Date.Sampled, Quarter, NMD.Code, Campy.PCR) %>%
  mutate(Animal="Pigs", Campy.PCR=str_replace_all(Campy.PCR, " ", ""))

poultry <- read.csv("data/esr/Poultry details.csv", skip = 1, stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  select(CPH.No., Date.Sampled, Quarter, NMD.Code, Campy.PCR = Confirmed..PCR.) %>%
  mutate(Animal="Poultry", Campy.PCR=str_replace_all(Campy.PCR, " ", ""))

produce <- read.csv("data/esr/Produce details.csv", stringsAsFactors = FALSE) %>%
  extract(PHL.Sample.Number, into='CPH.No.', regex="CPH[0]*([0-9A-B]+)", convert=FALSE) %>%
  mutate(Quarter = NA, NMD.Code = NA, Campy.PCR = NA) %>%
  select(CPH.No., Date.Sampled = Date.Received, Quarter, NMD.Code, Campy.PCR) %>%
  mutate(Animal="Produce", Campy.PCR=str_replace_all(Campy.PCR, " ", ""))
  
animals <- bind_rows(bobby, porcine, poultry, produce)

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
         !Antimicrobial %in% c("BAC", "TYL"))

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

glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|PlantEffect), family='binomial', data=dat_glm)

fit <- stan_glmer(cbind(Resistant, Susceptible) ~ AMCSpeciesAnimal + (1|PlantEffect), family='binomial', data=dat_glm)

out <- data.frame(dat_glm, t(apply(posterior_predict(fit), 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))), check.names=FALSE) %>%
  mutate(LI = `2.5%`/(Susceptible+Resistant),
         LC = `25%`/(Susceptible+Resistant),
         M = `50%`/(Susceptible+Resistant),
         UC = `75%`/(Susceptible+Resistant),
         UI = `97.5%`/(Susceptible+Resistant)) %>%
  mutate(Plant = factor(paste(Animal, NMD)))

ggplot(out %>% filter(Species == "E.coli")) +
  geom_segment(aes(x=Plant, xend=Plant, y=LI, yend=UI)) +
  geom_segment(aes(x=Plant, xend=Plant, y=LC, yend=UC, col=Animal), size=1) +
  geom_point(aes(x=Plant, y=M, col=Animal), size=3, shape=21, fill='white') +
  guides(col='none') +
  theme_bw() + ylab("Resistance") +
  scale_y_continuous(labels=scales::percent_format()) +
  facet_wrap(~Antimicrobial, scales = 'free_x')

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


library(lme4)
fit2 <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial', data=ec %>% filter(Animal == "Calves"))
icc(fit2)

# TODO: It seems poultry has a high ICC?? Need to look into this further,
# and adapt the design effect as needed (way fewer clusters there that
# are larger, right?)

#cv=0.8
#mean=33

# both pigs + calves do not.
fit2 <- glmer(cbind(Resistant, Susceptible) ~ (1|AMCSpeciesAnimal) + (1|NMD), family='binomial', data=dat_glm %>% filter(Species == "C.coli", Animal == "Poultry"))
icc(fit2)

# working out cv/mean cluster size
dat_glm %>% filter(Species == "C.jejuni", Animal == "Poultry", Antimicrobial=="CIP") %>%
  group_by(NMD) %>% dplyr::summarize(s = sum(Resistant+Susceptible)) %>%
  pull(s) %>% cv

# using same data as we had normally. This assumes equal sample sizes...
crtpwr.2prop(alpha=0.05, power=NA, m=24, n=12.5, icc=0.02, cv=0.56, p1=0.45, p2=0.55)

crtpwr.2prop(alpha=0.05, power=NA, m=24, n=12.5, icc=0.02, cv=0.56, p1=0.01, p2=0.06)
