rm(list = ls())

library(tidyverse)
library(nleqslv)
library(MASS)


# simulation study :

   # exp
exp_sim = function(n, t){
   df1 = 0  # data store
   df2 = matrix(ncol = 3, nrow = length(n))  # D max store
   R = 10000
   
   for(i in 1:length(n)){
      M = matrix(rexp(R*n[i], rate = 1/t), nrow = R, byrow = T)
      Z = sqrt(n[i])*(apply(M, 1, mean) - t)/t
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   df1 = df1[,-1]
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                       names_to = 'sample_size',
                       values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>%
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 30, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(1.2,0.8,1.2,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black',
                                      size = 20)) -> p
   
   rownames(df2) = paste(1:length(n))
   colnames(df2) = c('n','D_max','P value')
   return(list(df2,p))
}

exp_sim(c(10,30,100,200), 20)

setwd('D:/CLT Project')
ggsave(height = 9, width = 16,
       filename = 'Exp_sim(5).png')



   # beta
beta_sim = function(n, a, b){
   df1 = 0
   df2 = matrix(ncol = 3, nrow = length(n))
   R = 10000
   
   m = a/(a+b)
   v = (a*b)/(((a+b)^2)*(a+b+1))
   
   for(i in 1:length(n)){
      M = matrix(rbeta(R*n[i], a, b), nrow = R, byrow = T)
      Z = sqrt(n[i])*(apply(M, 1, mean) - m)/sqrt(v)
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   df1 = df1[,-1]
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                       names_to = 'sample_size',
                       values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>%
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 30, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(1.2,0.8,1.2,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black', 
                                      size = 20)) -> p
   
   colnames(df2) = c('n','D_max','P value')
   rownames(df2) = paste(1:length(n))
   return(list(p,df2))
}

beta_sim(c(10,30,100,200), 1, 5)

setwd('D:/CLT Project')
ggsave(height = 9, width = 16,
       filename = 'Beta_sim(1,5).png')



   # lognormal unction:

logn_sim = function(n, p1, p2){
   R = 10000
   df1 = 0
   df2 = matrix(ncol = 3, nrow = length(n))
   set.seed(1)
   
   m = exp(p1 + (p2^2)/2)
   v = exp(2*p1 + p2^2)*(exp(p2^2) - 1)
   
   for(i in 1:length(n)){
      M = matrix(rlnorm(R*n[i], meanlog = p1,sdlog = p2), 
                 nrow = R, byrow = T)
      Z = sqrt(n[i])*(apply(M, 1, mean) - m)/sqrt(v)
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   df1 = df1[,-1]; df
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                       names_to = 'sample_size',
                       values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>%
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 40, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(.~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(1.2,0.8,1.2,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black', 
                                      size = 20)) -> p
   colnames(df2) = c('n','D_max','P value')
   rownames(df2) = paste(1:length(n))
   return(list(p,df2))
}

logn_sim(c(10,30,100,200), 5, 1)

ggsave(height = 9, width = 16,
       filename = 'log_sim(0,0.5).png')




   # pareto function:

par_sim = function(n, a, x0){    # t = theta, pareto parameter
   R = 10000
   df1 = 0
   df2 = matrix(ncol = 3, nrow = length(n))
   set.seed(1)
   
   m = a*x0/(a-1)
   v = a*(x0^2)/((a-2)*(a-1)^2)
   
   for(i in 1:length(n)){
      M = matrix(rPareto(R*n[i], alpha = a, t = x0), 
                 nrow = R, byrow = T)
      Z = sqrt(n[i])*(apply(M, 1, mean) - m)/sqrt(v)
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   df1 = df1[,-1]
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                       names_to = 'sample_size',
                       values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>% 
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 40, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(1.2,0.8,1.2,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black', 
                                      size = 20)) -> p
   colnames(df2) = c('n','D_max','P value')
   rownames(df2) = paste(1:length(n))
   return(list(p,df2))
}

par_sim(c(10,30,100,200), 20, 5)

ggsave(height = 9, width = 16,
       filename = 'pareto_sim(20,5).png')



   # mixture distribution:

d1 = function(N, a, m1, s1, m2, s2){
   x = runif(N)
   ifelse(x <= a, rnorm(N,m1,s1), rnorm(N,m2,s2)) %>% 
      return()
}

mix_sim = function(n, a, m1, s1, m2, s2){
   R = 10000
   df1 = 0
   df2 = matrix(ncol = 3, nrow = length(n))
   set.seed(1)
   
   m = a*m1 + (1-a)*m2
   v = a*(s1^2 + m1^2) + (1-a)*(s2^2 + m2^2) - (a*m1 + (1-a)*m2)^2
   
   for(i in 1:length(n)){
      M = matrix(d1(R*n[i], a, m1, s1, m2, s2), 
                 nrow = R, byrow = T)
      Z = sqrt(n[i])*(apply(M, 1, mean) - m)/sqrt(v)
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   df1 = df1[,-1]
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                       names_to = 'sample_size',
                       values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>% 
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 40, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(2,0.8,1.8,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black', 
                                      size = 20)) -> p
   colnames(df2) = c('n','D_max','P value')
   rownames(df2) = paste(1:length(n))
   return(list(p,df2))
}

mix_sim(c(20,30,100,200), 0.2, 0, 1, 15, 1)

ggsave(height = 9, width = 16,
       filename = '0.2_0_1_15_1.png')


   # Application:
bw = read.csv("D:/Datasets/US_Baby_weights/us_new_born_2million.csv")
View(bw)

d = bw %>% filter(state == 'FL') %>% 
   select(weight_pounds) %>% na.omit()
colnames(d) = c('wt')

sum(is.na(d))
summary(d)
dim(d)
head(d,10)
pm = mean(d$wt); pm

ggplot(d, aes(wt)) + geom_histogram(bins = 40,
                                    colour = 'black',
                                    fill = 7,
                                    aes(y = ..density..)) +
   theme_minimal() + 
   labs(title = 'Distribution of baby weights in the state',
        x = 'Baby weights (lbs)', y = 'Frequency density') +
   theme(plot.margin = unit(c(1,2,1,2), 'cm'),
         plot.title = element_text(face = 'bold',
                                   size = 28, hjust = 0.5),
         axis.title = element_text(face = 'bold', size = 20),
         axis.text = element_text(face = 'bold', size = 16),
         axis.title.x = element_text(vjust = -3),
         axis.title.y = element_text(vjust = 5)) + 
   scale_x_continuous(n.breaks = 8) + 
   geom_vline(xintercept = pm, colour = 'red', lty = 2,
              lwd = 1.5)

setwd("D:/CLT Project")
ggsave(filename = "population_dist.png",
       height = 9, width = 16)




R = 10000
n = 10
s1 = 0

for(i in 1:R)
   s1[i] = mean(sample(d$wt, replace = T, size = 10))

m1 = mean(s1); m1
s1 = as.data.frame(s1)
sd(s1$s1)

ggplot(s1, aes(s1)) + geom_histogram(bins = 40,
                                    colour = 'black',
                                    fill = 'lightblue',
                                    aes(y = ..density..)) +
   theme_minimal() + 
   labs(title = 'n = 10',
        x = 'Sample means (lbs)', y = 'Frequency density') +
   theme(plot.title = element_text(face = 'bold',
                                   size = 20, hjust = 0.5),
         plot.margin = unit(c(1,2,1,2), 'cm'),
         axis.title = element_text(face = 'bold', size = 20),
         axis.text = element_text(face = 'bold', size = 16),
         axis.title.x = element_text(vjust = -3),
         axis.title.y = element_text(vjust = 5)) + 
   scale_x_continuous(n.breaks = 8) +
   geom_vline(xintercept = c(pm,m1), colour = c('red','darkgreen'),
              lty = c(2,1))
   

setwd("D:/CLT Project")
ggsave(filename = "dist1.png",
       height = 9, width = 16)







s2 = 0
n = c(20,30,50,100,200,500)
R = 10000

for(i in 1:length(n))
   for(j in 1:R)
      s2[(i-1)*R+j] = mean(sample(d$wt, size = n[i], replace = T))


s3 = data.frame('x' = s2, 'n' = rep(n, each = R))
s3$n = as.factor(s3$n)
View(s3)

a = aggregate(x ~ n, mean, data = s3); a
colnames(a) = c('Sample size', 'means')
a


name = c('20' = 'n = 20',
         '30' = 'n = 30',
         '50' = 'n = 50',
         '100' = 'n = 100',
         '200' = 'n = 200',
         '500' = 'n = 500')


ggplot(s3, aes(x, fill = n)) +
   geom_histogram(aes(y = ..density..),
                  bins = 40, colour = 'black') + 
   facet_wrap(.~n, 
     labeller = labeller(n = as_labeller(name, label_context))) +
   labs(x = 'Sample means', y = 'Frequency density') +
   theme_bw() +
   theme(plot.margin = unit(c(2,0.8,2,0.8), "cm"),
         axis.title = element_text(face = 'bold', size = 20),
         axis.title.x = element_text(vjust = -3),
         axis.title.y = element_text(vjust = 5),
         axis.text = element_text(face = 'bold', size = 14),
         legend.position = "none",
         strip.text = element_text(face = 'bold',
                                   colour = 'black', size = 20))+
   scale_y_continuous(n.breaks = 8)

setwd("D:/CLT Project")
ggsave(filename = "dist2.png",
       height = 9, width = 16)


#==================================================================
# Cauchy test :
set.seed(0)

X = rcauchy(161,0.2,1)
Xm = median(X); Xm
k = 80

pdf = function(x){
   k_fact = factorial(2*k + 1)/(factorial(k)^2)
   ((0.5 + (1/pi)*atan(x))^k)*
      ((0.5 - (1/pi)*atan(x))^k)/(pi*(1+x^2))*k_fact
}


f = function(c) (c(integrate(pdf, c, Inf)$value - 0.05))


c_val = nleqslv(0, f)$x; c_val


# Power curve :

f = function(mu)   # power function
{
   k_fact = factorial(2*k + 1)/(factorial(k)^2)
   
   f = function(x)
   {
      ((0.5 + (1/pi)*atan(x-mu))^k)*
         ((0.5 - (1/pi)*atan(x-mu))^k)/(pi*(1+(x-mu)^2))*k_fact
   }
   I = integrate(f, 0.205, Inf)$value
   return(I)
}

mu = seq(0,0.7, 0.001)
p1 = 0

for(i in 1:length(mu)) (p1[i] = f(mu[i]))

p2 = pnorm(8.077799*mu - 1.644854)

df = data.frame('mu' = rep(mu,2),
                'power' = c(p1,p2),
                'Type' = rep(c('Exact distribution',
                               'Asymptotic distribution'),
                             each = length(mu)))

# View(df)

my_col = c('red', 'blue')
df %>% ggplot(aes(x = mu, y = power, colour = Type)) +
   geom_line(aes(colour = Type, linetype = Type), 
             lwd = 2, alpha = 0.4) +
   theme_minimal() +
   labs(x = 'Parameter values',
        y = 'Power values')+
   theme(axis.title = element_text(face = 'bold',
                                   size = 18),
         plot.margin = unit(c(0.7,0.7,0.7,1.2), 'cm'),
         legend.title = element_text(face = 'bold',
                                     size = 16),
         legend.text = element_text(face = 'bold', 
                                    size = 14),
         axis.text = element_text(face = 'bold',
                                  size = 14),
         axis.title.x = element_text(vjust = -3),
         axis.title.y = element_text(vjust = 5)) +
   scale_x_continuous(n.breaks = 12) +
   scale_y_continuous(n.breaks = 14, minor_breaks = NULL) +
   geom_hline(yintercept = 0, lty = 2, lwd = 1) +
   geom_vline(xintercept = 0, lty = 2, lwd = 1) +
   scale_color_manual(values = my_col)

setwd("D:/CLT Project")
ggsave(filename = 'power curves.png',
       height = 9, width = 16)


#========================================================================
# Confidence interval :-

   # exact approach
CI_func = function(n,x)
{
   a = 0.05
   
   fl = function(l){
      f = function(t){
         (t^(x-1))*((1-t)^(n-x))
      }
         (integrate(f, 0, l)$value)/beta(x, n-x+1) - a/2
   }
   c1 = nleqslv(x/n, fl)$x
   
   fu = function(u){
      f = function(t){
         (t^x)*((1-t)^(n-x-1))
      }
         (integrate(f, 0, u)$value)/beta(x+1, n-x) - (1-a/2)
   }
   c2 = nleqslv(x/n, fu)$x
   
   return(c(c1,c2))
}

coeff_func = function(n, condition){
   R = 1000; CI = 0
   count1 = 0; count2 = 0
   
   set.seed(1)
   for(i in 1:R){
      x = rbinom(1, n, 0.5)
      CI = CI_func(n, x)
      if(CI[1] < 0.5 & CI[2] > 0.5) (count1 = count1 + 1)
   }
   
   set.seed(1)
   for(i in 1:R){
      x = rbinom(1, n, 0.5)
      c1 = (x/n) - qnorm(0.975)*sqrt(x*(n-x)/(n^3))
      c2 = (x/n) + qnorm(0.975)*sqrt(x*(n-x)/(n^3))
      if(c1 < 0.5 & c2 > 0.5) (count2 = count2 + 1)
   }
   
   if(condition == 'exact'){
      return(count1/R)
   }else if(condition == 'asym'){
      return(count2/R)
   }
}

N = seq(20,300,10)
cov_exact = 0
cov_asym = 0

for(i in 1:length(N)){
   cov_exact[i] = coeff_func(N[i], 'exact')
   cov_asym[i] = coeff_func(N[i], 'asym')
}

df = tibble(N, 'Exact' = cov_exact, 
      'Asymptotic' = cov_asym)


my_col = c('red','blue')
df %>% pivot_longer(Exact:Asymptotic,
                    names_to = 'Approach type',
                    values_to = 'coeff') %>% 
   ggplot(aes(x = N, y = coeff, color = `Approach type`,
              linetype = `Approach type`)) +
   geom_line(lwd = 1) + geom_point(color = 'black')+
   theme_bw() +
   labs(x = 'Values of n', y = 'Empirical coverage') +
   theme(axis.title = element_text(face = 'bold', size = 16),
         axis.text = element_text(face = 'bold', size = 12),
         plot.margin = unit(c(1.5,0.7,3,1), 'cm'),
         axis.title.x = element_text(vjust = -3),
         axis.title.y = element_text(vjust = 5),
         legend.title = element_text(face = 'bold', size = 16),
         legend.text = element_text(face = 'bold', size = 14)) +
   scale_x_continuous(n.breaks = 12) +
   scale_y_continuous(n.breaks = 10) +
   scale_color_manual(values = my_col) +
   geom_hline(yintercept = 0.95, lwd = 1, color = 'grey',
              lty = 5)

setwd("D:/CLT Project")
ggsave(filename = "CI_plot(new).png",
       height = 9, width = 16)


#====================================================================
n = 30:150
R = 10000
p = 0.8

cov1 = 0
cov2 = 0

set.seed(1)
for(i in 1:length(n))
{
   X = rbinom(R, n[i], p)
   lcl1 = X/n[i] - qnorm(0.975)*sqrt(X*(n[i] - X)/n[i]^3)
   ucl1 = X/n[i] + qnorm(0.975)*sqrt(X*(n[i] - X)/n[i]^3)
   
   lcl2 = (sin(asin(sqrt(X/n[i])) - qnorm(0.975)/(2*sqrt(n[i]))))^2
   ucl2 = (sin(asin(sqrt(X/n[i])) + qnorm(0.975)/(2*sqrt(n[i]))))^2
   
   cov1[i] = length(p[p > lcl1 & p < ucl1])/R
   cov2[i] = length(p[p > lcl2 & p < ucl2])/R
}

matplot(n, cbind(cov1, cov2), type = 'l', 
        col = c('red','blue'), lwd = 3)
abline(h = 0.95)

#============================================================================
# correlated
x = rexp(200, 10)

cov_func = function(k){
	M = matrix(ncol = k, nrow = k)
	diag(M) = rep(1,k)
	
	for(i in 1:(k-1)){
		r = round(rev(sort(rexp(k-i, 15))),3)
		M[i,(i+1):k] = r
		M[(i+1):k,i] = r
	}
	return(M)
}

corr_sim = function(n){   
   R = 10000
   df1 = 0
   df2 = matrix(ncol = 3, nrow = length(n))
   
   for(i in 1:length(n)){
      cov_mat = cov_func(n[i])
      M = mvrnorm(10000, rep(0,n[i]), cov_mat)
      
      Z = sqrt(n[i])*apply(M, 1, mean)/sqrt(1+
                        (sum(cov_mat)-n[i])/n[i])
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   df1 = df1[,-1]
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                        names_to = 'sample_size',
                        values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>% 
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 40, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(1.2,0.8,1.2,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black', 
                                      size = 20)) -> p
   colnames(df2) = c('n','D_max','P value')
   rownames(df2) = paste(1:length(n))
   return(list(p,df2))
}

setwd('D:/CLT Project')
corr_sim(c(10,30,100,200))

ggsave(height = 9, width = 16,
       filename = 'corr.png')
#====================================================================
# non-identical case

nonid_func = function(n){
   df1 = 0  # data store
   df2 = matrix(ncol = 3, nrow = length(n))  # D max store
   R = 10000
   set.seed(1)
   
   for(i in 1:length(n)){
      
      m = runif(n[i],0,10)
      sd = runif(n[i],1,3)
      
      M = matrix(ncol = n[i], nrow = R)
      for(j in 1:R) (M[j,] = rnorm(n[i], m, sd))
      
      Z = (apply(M, 1, mean) - mean(m))/sqrt(sum(sd^2)/(n[i]^2))
      k = ks.test(Z, 'pnorm')
      df2[i,] = c(n[i], round(k$statistic,5), 
                  round(k$p.value, 5))
      df1 = df1 %>% bind_cols(Z)
   }
   
   df1 = df1[,-1]
   colnames(df1) = paste('n =', n)
   
   df1 %>% pivot_longer(everything(), 
                        names_to = 'sample_size',
                        values_to = 'x') %>% 
      mutate(sample_size = factor(sample_size,
                                  levels = paste('n =',n))) %>%
      ggplot(aes(x = x, fill = sample_size)) + 
      geom_histogram(aes(y = ..density..),
                     bins = 30, color = 'black') +
      stat_function(fun = dnorm, args = list(0,1),
                    lty = 2, lwd = 1) +
      facet_grid(~sample_size) +
      labs(x = '', y = 'Frequency density') +
      theme_bw() +
      theme(plot.margin = unit(c(1.2,0.8,1.2,0.8), "cm"),
            axis.title.y = element_text(vjust = 5),
            axis.text = element_text(face = 'bold', size = 14),
            legend.position = "none",
            strip.text = element_text(face = 'bold',
                                      colour = 'black',
                                      size = 20)) -> p
   
   rownames(df2) = paste(1:length(n))
   colnames(df2) = c('n','D_max','P value')
   return(list(df2,p))
}

nonid_func(c(10,30,100,200))

setwd('D:/CLT Project')
ggsave(height = 9, width = 16,
       filename = 'nonid.png')








