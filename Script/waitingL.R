# Date: Jul 6, 2015 $
# Author: tjiang $
# Purpose: produce waitting length #

e = 1e-9
k = 8
p = 0.15

t1 = log(e)
t2 = log(1- (1-p)^k)

w1 = 1/p
w2 = k*(1-p)^k
w3 = 1 - (1-p)^k

wl = t1/t2 * (w1 - w2/w3)

print(c(k,wl))
