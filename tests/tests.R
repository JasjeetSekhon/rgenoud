suppressMessages(library(rgenoud))

#maximize the sin function
sin1 <- genoud(sin, nvars=1, max=TRUE, unif.seed=912821, int.seed=93058, print.level=0)
sin1$value <- signif(sin1$value,6)
sin1$par <- signif(sin1$par,6)
sin1$gradients <- signif(sin1$gradients,6)
print(sin1)

#minimize the sin function
sin2 <- genoud(sin, nvars=1, max=FALSE, unif.seed=912821, int.seed=93058, print.level=0)
sin2$value <- signif(sin2$value,6)
sin2$par <- signif(sin2$par,6)
sin2$gradients <- signif(sin2$gradients,6)
print(sin2)

