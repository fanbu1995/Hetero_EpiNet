# Jan 27, 2024: analyze "posterior predictive checks"

library(dplyr)
library(ggplot2)
library(tidyr)

pathHeader = "~/Downloads/Hetero-EpiNet-postCheck/postSimulate"
numSequences = c(0:49)

#i = 0

allSumm = NULL

for(i in numSequences){
  thePath = paste0(pathHeader,i)
  events = readr::read_csv(file.path(thePath, "dat.csv"))
  
  # count event types
  h1Summ = events %>% filter(time <= 35) %>%
    summarize(internalInfection = sum(event %in% c(9,10)),
              externalInfection = sum(event %in% c(109, 110)),
              activation = sum(event %in% c(3:5)),
              deletion = sum(event %in% c(6:8)),
              tmax = max(time)) %>%
    mutate(rep = i,
           period = 0)
  
  h2Summ = events %>% filter(time > 35)
  if(nrow(h2Summ) == 0){
    h2Summ = NULL
  }else{
    h2Summ = h2Summ %>%
      summarize(internalInfection = sum(event %in% c(9,10)),
                externalInfection = sum(event %in% c(109, 110)),
                activation = sum(event %in% c(3:5)),
                deletion = sum(event %in% c(6:8)),
                tmax = max(time)-35) %>%
      mutate(rep = i,
             period = 1)
  }
  
  allSumm = bind_rows(allSumm, h1Summ, h2Summ)
}
## look at total infection case counts
allSummInfec = allSumm %>% group_by(rep) %>%
  summarize(internalInfection = sum(internalInfection),
            externalInfection = sum(externalInfection)) %>%
  ungroup()

## turn into long format
allSumm_long = allSumm %>% 
  pivot_longer(internalInfection:deletion)

## plot distribution of infection cases and compare with truth
ggplot(allSummInfec, aes(x=internalInfection)) + 
  geom_density(fill = "#FBC511", alpha = 0.5) +
  geom_vline(xintercept = 18, linewidth = 1.5)+
  labs(y = "", x="internal infection") +
  theme_bw(base_size = 15)

ggplot(allSummInfec, aes(x=externalInfection))  + 
  geom_density(fill = "#69AED5", alpha = 0.5) +
  geom_vline(xintercept = 16, linewidth = 1.5)+
  labs(y = "", x= "external infection") +
  theme_bw(base_size = 15)


