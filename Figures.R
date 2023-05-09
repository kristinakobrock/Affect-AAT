# title: "Figures"
# author: "Kristina Kobrock"
# date: "2023-02-02"

# creating new shiny figures

# libraries
library(aida)
library(ggplot2)
library(brms)
library(tidyverse)
library(tidybayes)
library(bayestestR)

# theme for plotting
theme_set(theme_aida())

# global color scheme / non-optimized
# Green RYB: 5EB045, Barn red: 811F0E, Charcoal: 424B54, maximum yellow red: FFBA49, copper crayola: DE8F6E
#green <- "#5EB045"
#red <- "#811F0E"
gray <- "#424B54"
yellow <- "#FFBA49"
copper <- "#DE8F6E"
#black <- "#000000"
# ["#399283", "#6ce9d3", "#074d65", "#a7dcf9", "#88075f", "#da82b1"]
green1 <- "#399283"
red1 <- "#88075f"
green2 <- "#6ce9d3"
red2 <- "#da82b1"
blue1 <- "#208eb7"
blue2 <- "#a7dcf9"
orange <- "#f4a95c"
black <- "#0000ff"
#project_colors = c(red1, green1, red2, green2, blue1, blue2, gray, yellow, copper, black)
#group.colors <- c(A = "#333BFF", B = "#CC6600", C ="#9633FF", D = "#E2FF33", E = "#E3DB71")
condition.colors <- c(Incongruent = orange, Congruent = blue1)
stimulus_valence.colors <- c(Positive = green2, Negative = red2)
emotional_priming.colors <- c(Positive = green1, Negative = green2)

# setting theme colors globally
#scale_colour_discrete <- function(...) {
#  scale_colour_manual(..., values = project_colors)
#}
#scale_fill_discrete <- function(...) {
#  scale_fill_manual(..., values = project_colors)
#} 

backtransform_RT <- function(list, df=df_aat_log) {
  round(exp(list * attr(df$log_RT_s, 'scaled:scale') + attr(df$log_RT_s, 'scaled:center')),2)
}


# setting working directory and loading data
setwd("/Users/kkobrock/projects/Affect-AAT/Affect-AAT")
load("affect-aat.Rdata")
load("panas.Rdata")
load("panas_item.Rdata")
load("arousal.Rdata")

# mutate for better plot labels
df_aat_log <- df_aat_log %>% 
  mutate(condition = ifelse(condition == "congruent", "Congruent", "Incongruent"),
         stimulus_valence = ifelse(stimulus_valence == "positive", "Positive", "Negative"),
         emotional_priming = ifelse(mood_induction == "positive", "Positive", "Negative"))

# get model
#fit_stimulusModelRE <- brm(file = "CC_model_stimulusvalence")
fit_model <- brm(file = "stimulus_model")

library(bayestestR)
describe_posterior(fit_model, rope_range=rope_range(fit_model))

# get posterior draws
posterior_draws <- fit_model %>%
  # just taking the population level effects here
  spread_draws(b_Intercept, b_conditionincongruent,
               b_mood_inductionpositive, b_stimulus_valencepositive,
               `b_conditionincongruent:mood_inductionpositive`,
               `b_conditionincongruent:stimulus_valencepositive`,
               `b_mood_inductionpositive:stimulus_valencepositive`,
               `b_conditionincongruent:mood_inductionpositive:stimulus_valencepositive`) %>% 
  mutate(
    congruent_negative_neg = b_Intercept,
    incongruent_negative_neg = b_Intercept + b_conditionincongruent,
    congruent_positive_neg = b_Intercept + b_mood_inductionpositive,
    congruent_negative_pos = b_Intercept + b_stimulus_valencepositive,
    incongruent_positive_neg = b_Intercept + b_conditionincongruent + 
      b_mood_inductionpositive + `b_conditionincongruent:mood_inductionpositive`,
    incongruent_negative_pos = b_Intercept + b_conditionincongruent + 
      b_stimulus_valencepositive + `b_conditionincongruent:stimulus_valencepositive`,
    congruent_positive_pos = b_Intercept + b_mood_inductionpositive + 
      b_stimulus_valencepositive + `b_mood_inductionpositive:stimulus_valencepositive`,
    incongruent_positive_pos = b_Intercept + b_conditionincongruent +
      b_mood_inductionpositive + b_stimulus_valencepositive +
      `b_conditionincongruent:mood_inductionpositive` + `b_conditionincongruent:stimulus_valencepositive` +
      `b_mood_inductionpositive:stimulus_valencepositive` + `b_conditionincongruent:mood_inductionpositive:stimulus_valencepositive`
  ) %>% 
  dplyr::select(congruent_negative_neg, congruent_negative_pos,
                congruent_positive_neg, congruent_positive_pos,
                incongruent_negative_neg, incongruent_negative_pos,
                incongruent_positive_neg, incongruent_positive_pos) %>%
  gather(key = "parameter", value = "posterior")

posteriors <- posterior_draws %>% 
  mutate(condition = ifelse(parameter == "congruent_negative_neg" | parameter == "congruent_positive_neg" |
                              parameter == "congruent_negative_pos" | parameter == "congruent_positive_pos",
                            "Congruent", "Incongruent"),
         emotional_priming = ifelse(parameter == "congruent_negative_neg" | parameter == "incongruent_negative_neg" |
                                   parameter == "congruent_negative_pos" | parameter == "incongruent_negative_pos",
                                 "Negative", "Positive"),
        stimulus_valence = ifelse(parameter == "congruent_negative_neg" | parameter == "congruent_positive_neg" |
                                    parameter == "incongruent_negative_neg" | parameter == "incongruent_positive_neg",
                                  "Negative", "Positive")
  ) %>% 
  select(condition, emotional_priming, stimulus_valence, posterior) 

#-----------------------------------
# Figure 3: AAB: condition main effect
hypothesis(fit_model, "conditionincongruent > 0")

# compute CrIs as HDIs
posterior_CrIs_AAB <- posteriors %>% 
  group_by(condition) %>% 
  summarise(mean_posterior = mean(posterior),
            `95lowerCrI` = HDInterval::hdi(posterior, credMass = 0.95)[1],
            `95higherCrI` = HDInterval::hdi(posterior, credMass = 0.95)[2])

#plot
plot_AAB <- ggplot() +
  geom_jitter(aes(x = condition, y = log_RT_s), data = df_aat_log, alpha = 0.5, color=yellow) +
  geom_errorbar(aes(x = condition, ymin=`95lowerCrI`, ymax=`95higherCrI`, group=condition), data = posterior_CrIs_AAB, color=gray) + 
  geom_line(aes(x = condition, y = mean_posterior, group = NA), data = posterior_CrIs_AAB, linetype=2, color=gray) + 
  geom_point(aes(x = condition, y = mean_posterior), data = posterior_CrIs_AAB, size=3, color=copper) +  
  ylim(-2, 2) +
  ylab("log reaction times (scaled)")
plot_AAB

plot_AABms <- ggplot() +
  geom_jitter(aes(x = condition, y = reaction_times, color=condition), data = df_aat_log, alpha = 0.5) +
  geom_errorbar(aes(x = condition, ymin=backtransform_RT(`95lowerCrI`), 
                    ymax=backtransform_RT(`95higherCrI`), group=condition), 
                data = posterior_CrIs_AAB, linewidth=0.8, color=gray) +
  geom_line(aes(x = condition, y = backtransform_RT(mean_posterior), group = NA), 
            data = posterior_CrIs_AAB, linetype=2, linewidth=0.8, color=gray) +
  geom_point(aes(x = condition, y = backtransform_RT(mean_posterior)), 
             data = posterior_CrIs_AAB, size=5, color=gray) +
  ylim(700, 1000) +
  ylab("Reaction Times (ms)") +
  xlab("Condition") +
  theme(text = element_text(size = 14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none") +
  # custom colors
  scale_colour_manual(values=condition.colors) 

plot_AABms

ggsave(filename = "Figures/AAB_main.jpg",
       plot = plot_AABms,
       device = "jpg",
       bg = "white",
       width = 240, 
       height = 140,
       units = "mm",
       dpi = 300)

#-----------------------------------
# Figure 4: stimulus valence main effect
hypothesis(fit_model, "stimulus_valencepositive < 0")
hypothesis(fit_model, "Intercept > 0")

# compute CrIs as HDIs
posterior_CrIs_stimulus <- posteriors %>% 
  group_by(stimulus_valence) %>% 
  summarise(mean_posterior = mean(posterior),
            `95lowerCrI` = HDInterval::hdi(posterior, credMass = 0.95)[1],
            `95higherCrI` = HDInterval::hdi(posterior, credMass = 0.95)[2])

# plot
plot_stimulus <- ggplot() +
  geom_jitter(aes(x = stimulus_valence, y = log_RT_s, color=stimulus_valence), 
              data = df_aat_log, alpha = 0.15, width=0.3) +
  geom_errorbar(aes(x = stimulus_valence, ymin=`95lowerCrI`, ymax=`95higherCrI`, group=stimulus_valence), 
                color=gray, data = posterior_CrIs_stimulus, width=0.8, linewidth=0.8) +
  geom_line(aes(x = stimulus_valence, y = mean_posterior, group = NA), 
            data = posterior_CrIs_stimulus, color=gray, linetype=2, linewidth=0.8) +
  geom_point(aes(x = stimulus_valence, y = mean_posterior, color=stimulus_valence), 
             data = posterior_CrIs_stimulus, size=5) +
  ylim(-1, 1) +
  theme(legend.position = "bottom") +
  labs(y = "log reaction times (scaled)", 
       x = "",
       color = "stimulus valence")
plot_stimulus

plot_stimulus_ms <- ggplot() +
  geom_jitter(aes(x = stimulus_valence, y = reaction_times, color=stimulus_valence), 
              data = df_aat_log, alpha = 0.2, width=0.3) +
  geom_errorbar(aes(x = stimulus_valence, ymin=backtransform_RT(`95lowerCrI`), ymax=backtransform_RT(`95higherCrI`), 
                    group=stimulus_valence), color=gray, data = posterior_CrIs_stimulus, width=0.8, linewidth=0.8) +
  geom_line(aes(x = stimulus_valence, y = backtransform_RT(mean_posterior), group = NA), 
            data = posterior_CrIs_stimulus, color=gray, linetype=2, linewidth=0.8) +
  geom_point(aes(x = stimulus_valence, y = backtransform_RT(mean_posterior), color=stimulus_valence),
             data = posterior_CrIs_stimulus, size=5) +
  theme(legend.position = "right") +
  ylim(700, 1000) +
  labs(y = "reaction times (ms)", 
       x = "",
       color = "stimulus valence") +
  theme(text = element_text(size = 30),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))
plot_stimulus_ms

#-----------------------------------
# Figure 5: PANAS (arousal)

# get model
model_arousal <- brm(file = "panas_arousal_model") # need to set adapt_delta above 0.8

describe_posterior(model_arousal, "mean", rope_range = rope_range(model_arousal))
hypothesis(model_arousal, "mood_inductionpositive < 0")
hypothesis(model_arousal, "Intercept > 0")

# get posterior draws
posterior_draws_arousal <- model_arousal %>%
  # just taking the population level effects here
  spread_draws(b_Intercept, b_mood_inductionpositive) %>% 
  mutate(
    MI_negative = b_Intercept,
    MI_positive = b_Intercept + b_mood_inductionpositive
  ) %>% 
  dplyr::select(MI_negative, MI_positive) %>%
  gather(key = "parameter", value = "posterior")

posteriors_arousal <- posterior_draws_arousal %>% 
  mutate(mood_induction = ifelse(parameter == "MI_negative",
                            "negative", "positive")
  ) %>% 
  select(mood_induction, posterior) 

# compute CrIs as HDIs
posterior_CrIs_arousal <- posteriors_arousal %>% 
  group_by(mood_induction) %>% 
  summarise(mean_posterior = mean(posterior),
            `95lowerCrI` = HDInterval::hdi(posterior, credMass = 0.95)[1],
            `95higherCrI` = HDInterval::hdi(posterior, credMass = 0.95)[2])

# sumstats
df_arousal_agg <- df_arousal_exc %>%
  group_by(mood_induction) %>% 
  summarize(mean_arousal = mean(arousal_score))

# plot
plot_panas_arousal <- df_arousal_exc %>% 
  ggplot(aes(x = mood_induction, y = arousal_score, color = mood_induction)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(alpha = 0.3, size = 2) +#, position = position_jitter(width = 0.1)) +
  geom_line(aes(group = subject_id), alpha = 0.1, color = gray) +
  geom_line(aes(x = mood_induction, y = mean_posterior, group = NA), 
             data = posterior_CrIs_arousal, color = gray, linetype = 2, linewidth = 1) +
  geom_errorbar(aes(x = mood_induction, y = mean_posterior, ymin = `95lowerCrI`, ymax = `95higherCrI`),
                data = posterior_CrIs_arousal, width = 0.2, color = gray) +
  geom_point(aes(x = mood_induction, y = mean_posterior),
             data = posterior_CrIs_arousal, size = 5) +
  ylab("arousal increase") +
  xlab("mood induction") +
  ylim(-21, 23) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 30),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))
plot_panas_arousal

#-----------------------------------
# Figure 6: PANAS (affect)

# mutate such that negative items are on negative scale
df_affect_plot <- df_affect_excs %>% 
  mutate(item_score = ifelse(item_affect == "NA", -item_score, item_score),
         panas_nr = ifelse(panas_nr == 1, "(1) before", "(2) after"))

# sumstats
panas_affect_agg <- df_affect_plot %>% 
  group_by(item_affect, mood_induction, panas_nr) %>% 
  summarize(item_mean = mean(item_score))
panas_affect_agg

# plot
plot_panas_affect <- df_affect_plot %>% 
  ggplot(aes(x = panas_nr, y = item_score, color = mood_induction), alpha = 0.1) +
  facet_grid(cols = vars(mood_induction)) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, color = mood_induction), alpha = 0.001) +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.2)) +
  #geom_line(aes(group = interaction(subject_id, item_affect)), color = gray, alpha = 0.2) +
  geom_point(aes(x = panas_nr, y = item_mean), data = panas_affect_agg, size = 4) +
  geom_line(aes(x = panas_nr, y = item_mean, group = item_affect), data = panas_affect_agg, 
            linetype = 1, linewidth = 1) +
  scale_y_continuous(limits = c(-5, 6), 
                     breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6),
                     labels = c("extremely (-5)", "quite a bit (-4)","moderately (-3)",
                              "a little (-2)","very slightly or not at all (-1)",
                              "negative affect:", "very slightly or not at all (1)", "a little (2)", 
                              "moderately (3)", "quite a bit (4)", "extremely (5)", "positive affect:")) +
  ylab("") +
  xlab("emotional priming") +
  theme(legend.position = "none")
plot_panas_affect

# plot summed affect scores
df_affect_scores <- df_affect_excs %>% 
  mutate(item_affect = ifelse(item_affect == "PA", "positive affect", "negative affect"),
         panas_nr = ifelse(panas_nr == "1", "(1) before", "(2) after")) %>% 
  group_by(subject_id, mood_induction, item_affect, panas_nr) %>% 
  summarize(affect_score = sum(item_score))

#sumstat
df_scores_agg <- df_affect_scores %>% 
  group_by(mood_induction, item_affect, panas_nr) %>% 
  summarize(mean = mean(affect_score))

# plot
plot_panas_affect_scores <- df_affect_scores %>% 
  ggplot(aes(x = panas_nr, y = affect_score, color = item_affect), alpha = 0.1) +
  facet_grid(item_affect ~ mood_induction) +
  #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, color = mood_induction), alpha = 0.001) +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.2)) +
  #geom_line(aes(group = interaction(subject_id, item_affect)), color = gray, alpha = 0.2) +
  geom_point(aes(x = panas_nr, y = mean), data = df_scores_agg, size = 4) +
  geom_line(aes(x = panas_nr, y = mean, group = item_affect), data = df_scores_agg, 
            linetype = 1, linewidth = 1) +
  ylab("affect score") +
  xlab("emotional priming") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 30),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))
plot_panas_affect_scores

#-----------------------------------
# Figure 7: Interaction: emotional priming and stimulus valence

# hypotheses
hypothesis(fit_model, "conditionincongruent > 0")
hypothesis(fit_model, "stimulus_valencepositive < 0")
hypothesis(fit_model, "mood_inductionpositive:stimulus_valencepositive < 0")

# compute CrIs as HDIs
posterior_CrIs_int <- posteriors %>% 
  group_by(condition, mood_induction, stimulus_valence) %>% 
  summarise(mean_posterior = mean(posterior),
            `95lowerCrI` = HDInterval::hdi(posterior, credMass = 0.95)[1],
            `95higherCrI` = HDInterval::hdi(posterior, credMass = 0.95)[2])

df_aat_log1 <- df_aat_log %>% 
  mutate(mood_induction = ifelse(mood_induction == "p", "positive", "negative"))

# plot
plot_3way_ms <- ggplot() +
  # plotting condition, mood_induction and stimulus valence
  geom_point(aes(x = stimulus_valence, y = reaction_times, color = stimulus_valence),
         data = df_aat_log1, alpha = 0.05, size = 1, position = position_jitter(width = 0.3)) +
  facet_grid(mood_induction ~ condition) +
  #theme(strip.background = element_rect(fill = copper)) +
  #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, color = mood_induction), 
  #          data = df_aat_log1, alpha = 0.001) +
  #geom_point(alpha = 0.05, size = 1,
             # adding a little bit of jitter to make the points visible
  #           position = position_jitter(width = 0.3)) +
  # adding posterior means and a connecting line in order to make the differences between locations and conditions 
  # visible
  geom_point(data = posterior_CrIs_int,
             aes(x = stimulus_valence, y = backtransform_RT(mean_posterior)), color = gray, size = 4) +
  geom_line(data = posterior_CrIs_int,
            aes(x = stimulus_valence, y = backtransform_RT(mean_posterior), group = mood_induction), color = gray) +
  geom_errorbar(data = posterior_CrIs_int,
                aes(x = stimulus_valence, y = backtransform_RT(mean_posterior), 
                    ymin = backtransform_RT(`95lowerCrI`), 
                    ymax = backtransform_RT(`95higherCrI`), 
                    group = condition), color=gray, width = 0.4) +
  ylim(750, 950) + 
  # axes labels and titles
  labs(
    #title = "3-way interaction",
    y = "reaction times (ms)",
    x = "stimulus valence"
  ) +
  guides(fill = "none", color = "none") +
  theme(text = element_text(size = 30),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))

plot_3way_ms

#-----------------------------------
#-----------------------------------
# More plots (might go into the appendix)

# Item-wise PANAS plot with lines 
plot_panas_item <- df_affect_excs %>% 
  ggplot(aes(x = panas_nr, y = item_score, color = item_affect)) +
  geom_point(aes(x = panas_nr, y = item_score, group = subject_id), alpha = 0.1) +
  geom_line(aes(group = subject_id), alpha = 0.1) +
  facet_grid(mood_induction ~ panas_item)
plot_panas_item

# Item-wise PANAS plot with means
df_item_agg <- df_affect_excs %>% 
  group_by(mood_induction, panas_nr, panas_item, item_affect) %>% 
  summarize(mean = mean(item_score))

plot_panas_item_mean <- df_affect_excs %>% 
  ggplot(aes(x = panas_nr, y = item_score, color = item_affect)) +
  geom_point(aes(x = panas_nr, y = mean), data = df_item_agg) +
  geom_line(aes(x = panas_nr, y = mean, group = panas_item), data = df_item_agg, color = gray) +
  facet_grid(mood_induction ~ panas_item) +
  scale_y_continuous(breaks=c(1,2,3,4,5), limits = c(1,5)) +
  theme(text = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 8),
        legend.title = element_text(size=10),
        strip.text.y = element_text(size = 10)
        )
plot_panas_item_mean

#-----------------------------------
# plot difference distributions

