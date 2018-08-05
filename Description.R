library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(caret)
#library(pls) # for mvr
library(naivebayes)
library(ggfortify)
library(plotly)

setwd("~/PycharmProjects/BioHack18")

load("~/PycharmProjects/BioHack18/.RData")
#save.image("~/PycharmProjects/BioHack18/.RData")

patients <- read.csv(file = "RECORDS.txt", header = F, sep = "/",
                     stringsAsFactors = F)

controls <- read.csv(file = "database/CONTROLS", header = F, sep = "/",
                     stringsAsFactors = F)

(Nsubsribes <- patients %>% 
  group_by(V1) %>% 
  summarise(N = n()) %>% 
  as.data.frame() %>% 
  group_by(N) %>% 
  summarise(n()))

controls %>% 
  group_by(V1) %>% 
  summarise(N = n()) %>% 
  as.data.frame() %>% 
  group_by(N) %>% 
  summarise(n())

patients %>% 
  group_by(V1) %>% 
  summarise(N = n()) %>% 
  as.data.frame() %>% 
  group_by(N) %>% 
  summarise(NN = n()) %>% 
  as.data.frame() %>% 

  ggplot(aes(N, NN, fill = N))+ 
  geom_bar(stat = "Identity")+
  geom_text(aes(label = NN), nudge_y = 10, size = 8)+
  scale_x_continuous(name = "Количество наблюдений", breaks = c(1:7))+
  scale_y_continuous(name = "Количество пациентов")+
  scale_fill_continuous(low = 'green', high = 'red')+
  guides(fill = F)+
  theme_classic()+
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 23),
        axis.title = element_text(size = 23),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 23),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 16))

fil <- read.csv(file = "Unfiltered.txt", sep = "\t", header = F)

unfil <- read.csv(file = "Filtered.txt", sep = "\t", header = F)

df <- cbind(fil, unfil) %>% 
  rownames_to_column() %>% 
  mutate(rowname = as.numeric(rowname)*2)

df %>% 
  slice(200:550) %>% 
  gather(-rowname, key = key, value = value) %>% 
  
  ggplot(aes(rowname, value))+
  geom_line()+
  facet_wrap(~key, nrow = 2, labeller = as_labeller(c(`V1` = "До фильтрации",
                                                         `V1.1` = "После фильтрации")))+
  
  scale_x_continuous("Время, мс")+
  scale_y_continuous("Напряжение", breaks = NULL)+
  
  theme_classic()+
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 23),
        axis.title = element_text(size = 23),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 23),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 16))
  

X <- read.csv(file = "totalX.txt", header = F, sep = "\t")
y <- read.csv(file = "label.txt", header = F, sep = "\t") %>% 
  transmute(y = V1)

Xb <- read.csv(file = "bootX.txt", header = F, sep = "\t")
yb <- read.csv(file = "boot_lab.txt", header = F, sep = "\t") %>% 
  transmute(y = V1)

X9b <- read.csv(file = "9bootX.txt", header = F, sep = "\t")
y9b <- read.csv(file = "9boot_lab.txt", header = F, sep = "\t") %>% 
  transmute(y = V1)

Xmb <- read.csv(file = "max_totalX.txt", header = F, sep = "\t")
ymb <- read.csv(file = "max_label.txt", header = F, sep = "\t") %>% 
  transmute(y = V1)

data <- cbind(X, y) %>% 
  mutate(y = factor(y))

datab <- cbind(Xb, yb) %>% 
  mutate(y = factor(y))

data9b <- cbind(X9b, y9b) %>% 
  mutate(y = factor(y))

datamb <- cbind(Xmb, ymb) %>% 
  mutate(y = factor(y))

models <-  list()

control <- trainControl(method = "cv", selectionFunction = "oneSE",
                        index = createMultiFolds(data$y,
                                                 k = 5,
                                                 times = 10))
controlb <- trainControl(method = "cv", selectionFunction = "oneSE",
                         index = createMultiFolds(datab$y,
                                                  k = 5,
                                                  times = 10))
control9b <- trainControl(method = "cv", selectionFunction = "oneSE",
                          index = createMultiFolds(data9b$y,
                                                   k = 5,
                                                   times = 10))
controlmb <- trainControl(method = "cv", selectionFunction = "oneSE",
                          index = createMultiFolds(datamb$y,
                                                   k = 5,
                                                   times = 10))

models$lda <- train(y ~ .,
                    data = data,
                    method = "lda", 
                    tuneLength = 5,
                    trControl = control,
                    preProcess = c("nzv", "center", "scale"))

models$lda1 <- train(y ~ .,
                     data = datab,
                     method = "lda", 
                     tuneLength = 5,
                     trControl = control,
                     preProcess = c("nzv", "center", "scale"))

models$nb1 <- train(y ~ .,
                   data = datab,
                   method = "naive_bayes", 
                   tuneLength = 8,
                   trControl = control,
                   preProcess = c("nzv", "center", "scale"))

models$nb <- train(y ~ .,
                    data = data,
                    method = "naive_bayes", 
                    tuneLength = 8,
                    trControl = control,
                    preProcess = c("nzv", "center", "scale"))

models$lda9b <- train(y ~ .,
                    data = data9b,
                    method = "lda", 
                    tuneLength = 7,
                    trControl = control,
                    preProcess = c("nzv", "center", "scale"))

models$nb9b <- train(y ~ .,
                    data = data9b,
                    method = "naive_bayes", 
                    tuneLength = 7,
                    trControl = control,
                    preProcess = c("nzv", "center", "scale"))

models$ldamb <- train(y ~ .,
                      data = datamb,
                      method = "lda", 
                      tuneLength = 7,
                      trControl = control,
                      preProcess = c("nzv", "center", "scale"))

models$nbmb <- train(y ~ .,
                     data = datamb,
                     method = "naive_bayes", 
                     tuneLength = 7,
                     trControl = control,
                     preProcess = c("nzv", "center", "scale"))

models$qdamb <- train(y ~ .,
                     data = datamb,
                     method = "qda", 
                     tuneLength = 7,
                     trControl = controlmb,
                     preProcess = c("nzv", "center", "scale"))

models$qda9b <- train(y ~ .,
                     data = data9b,
                     method = "qda", 
                     tuneLength = 7,
                     trControl = control9b,
                     preProcess = c("nzv", "center", "scale"))

models$qdab1 <- train(y ~ .,
                     data = datab,
                     method = "qda", 
                     tuneLength = 7,
                     trControl = controlb,
                     preProcess = c("nzv", "center", "scale"))

models$qda <- train(y ~ .,
                     data = data,
                     method = "qda", 
                     tuneLength = 7,
                     trControl = control,
                     preProcess = c("nzv", "center", "scale"))

models$svm <- train(y ~ .,
                    data = data,
                    method = "svmLinear", 
                    tuneLength = 7,
                    trControl = control,
                    preProcess = c("nzv", "center", "scale"))

models$svm1 <- train(y ~ .,
                    data = datab,
                    method = "svmLinear", 
                    tuneLength = 7,
                    trControl = controlb,
                    preProcess = c("nzv", "center", "scale"))

models$svm9b <- train(y ~ .,
                    data = data9b,
                    method = "svmLinear", 
                    tuneLength = 7,
                    trControl = control9b,
                    preProcess = c("nzv", "center", "scale"))

models$svmmb <- train(y ~ .,
                    data = datamb,
                    method = "svmLinear", 
                    tuneLength = 7,
                    trControl = controlmb,
                    preProcess = c("nzv", "center", "scale"))

bwplot(resamples(models), metric = "Accuracy")

nrow(unique(ymb))

bwplot(resamples(list("10 индивидов" = models$lda,
                 "18 индивидов" = models$lda1,
                 "18 индивидов (9)" = models$lda9b,
                 "268 индивидов" = models$ldamb)), metric = "Accuracy",
       main = list(label="Линейный дискриминантный анализ",cex=1.5),
       scales = list(y=list(cex = 1.5), x=list(cex=1.5)),
       strip=strip.custom(factor.levels=c("Точность", "Каппа")))

bwplot(resamples(list("10 индивидов" = models$nb,
                      "18 индивидов," = models$nb1,
                      "18 индивидов (9)" = models$nb9b,
                      "268 индивидов" = models$nbmb)), metric = "Accuracy",
       main = list(label="Наивный Байесовский классификатор",cex=1.5),
       scales = list(y=list(cex = 1.5), x=list(cex=1.5)),
       strip=strip.custom(factor.levels=c("Точность", "Каппа")))

bwplot(resamples(list("10 индивидов" = models$qda,
                      "18 индивидов" = models$qdab1,
                      "18 индивидов (9)" = models$qda9b,
                      "268 индивидов" = models$qdamb)), metric = "Accuracy",
       main = list(label="Квадратичный дискриминантный анализ",cex=1.5),
       scales = list(y=list(cex = 1.5), x=list(cex=1.5)),
       strip=strip.custom(factor.levels=c("Точность", "Каппа")))

bwplot(resamples(list("10 индивидов" = models$svm,
                      "18 индивидов" = models$svm1,
                      "18 индивидов (9)" = models$svm9b)), metric = "Accuracy",
       main = list(label="Метод опорных векторов",cex=1.5),
       scales = list(y=list(cex = 1.5), x=list(cex=1.5)),
       strip=strip.custom(factor.levels=c("Точность", "Каппа")))

bwplot(resamples(list("Линейный ДА" = models$ldamb,
                      "Наивный Байес" = models$nbmb,
                      "Квадратичный ДА" = models$qdamb)), metric = "Accuracy",
       main = list(label="Качество модели для полной выборки",cex=1.5),
       scales = list(y=list(cex = 1.5), x=list(cex=1.5)),
       strip=strip.custom(factor.levels=c("Точность", "Каппа")))

plot(varImp(models$qdamb), 4)

pca <- prcomp(datamb[,-10])

data9bn <- data9b %>% 
  transmute("R.y.I" = V2,
            "Q.x" = V3,
            "Q.y.I" = V4,
            "S.x" = V5,
            "S.y.I" = V6,
            "R.y.II" = V7,
            "Q.y.II" = V8,
            "S.y.II" = V9)

pcas <- prcomp(data9bn, scale. = T)

autoplot(pcas, data = data9b, colour = 'y',
         loadings = T, loadings.colour = 'blue', loadings.label = T,
         loadings.label.size = 10)+
  guides(colour=F)+
  scale_x_continuous(name = "Главная компонента 1 (28.52%)")+
  scale_y_continuous(name = "Главная компонента 2 (21.62%)")+
  theme_classic()+
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 23),
        axis.title = element_text(size = 23),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 23),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 16))+
  ggtitle("Выборка 3 в координатах главных компонент")


  
  