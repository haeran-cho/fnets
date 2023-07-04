library(fnets)
## Alessi, Barigozzi, and Capasso method for restricted models
abc <- factor.number(data.restricted, fm.restricted = TRUE)
print(abc)
plot(abc)

## Eigenvalue ratio method
er <- factor.number(data.restricted, method = "er", fm.restricted = TRUE)
print(er)
plot(er)

## Hallin and Liška method for unrestricted models
hl <- factor.number(data.unrestricted, fm.restricted = FALSE)
print(hl)
plot(hl)
