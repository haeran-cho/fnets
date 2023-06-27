library(fnets)
## Alessi, Barigozzi, and Capasso method for restricted models
x <- fnets::restricted
abc <- factor.number(x, fm.restricted = TRUE)
print(abc)
plot(abc)

## Eigenvalue ratio method
er <- factor.number(x, method = "er", fm.restricted = TRUE)
print(er)
plot(er)

## Hallin and Liška method for unrestricted models
x <- fnets::unrestricted
hl <- factor.number(x, fm.restricted = FALSE)
print(hl)
plot(hl)
