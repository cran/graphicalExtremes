## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
library(graphicalExtremes)
library(dplyr)
library(ggplot2)
library(igraph)

## ----include=FALSE, eval=FALSE------------------------------------------------
#  rightAlign <- function(txt){
#      paste0('<div style="text-align: right">', txt, '</div>')
#  }
#  now <- Sys.time()
#  knitr::knit_hooks$set(timeit = function(before, options, envir) {
#      if (before) {
#          now <<- Sys.time()
#      } else {
#          rightAlign(sprintf("\n_Code execution time: %.3f seconds._\n", Sys.time() - now))
#      }
#  })

## ----include=FALSE------------------------------------------------------------
par0 <- par()
knitr::knit_hooks$set(nomargin = function(before) {
    if(before){
        par0 <<- par(mar=c(0,0,0,0))
    } else{
        par(par0)
    }
})
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = '#>',
    error = TRUE,
    fig.align = 'center',
    timeit = TRUE
)
theme_set(
    theme_bw()
    + theme(
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = "white"),
        plot.caption = element_text( size = 7.5, hjust = 0, margin = margin(t = 15)),
        text = element_text(size = 11),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size = 0.25)
    )
)

## ----fig.width=6, fig.height=4------------------------------------------------
# Plot the physical locations of measuring stations, indicating flow volume by width
plotDanube(useStationVolume = TRUE, useConnectionVolume = TRUE, xyRatio = 1.5)

## ----nomargin=TRUE------------------------------------------------------------
# Plot only the flow graph
danube_flow <- getDanubeFlowGraph()
plotDanubeIGraph(graph = danube_flow)

## ----fig.show='hold', out.width='50%', fig.align=''---------------------------
# Plot some pairwise scatter plots
plotData <- as_tibble(danube$data_clustered)
ggplot(plotData) + geom_point(aes(x = X1, y = X2))
ggplot(plotData) + geom_point(aes(x = X22, y = X28))

## -----------------------------------------------------------------------------
X <- danube$data_clustered
Y <- data2mpareto(data = X, p = .8)

## -----------------------------------------------------------------------------
theta <- .5
Ysim <- rmpareto(n = 100, model = "logistic", d = 2, par = theta)

## -----------------------------------------------------------------------------
G <- cbind(c(0, 1.5), c(1.5, 0))
Ysim <- rmpareto(n = 100, model = "HR", par = G)

## -----------------------------------------------------------------------------
chi_hat <- emp_chi(data = X, p = .8)

## -----------------------------------------------------------------------------
chi_hat <- emp_chi(data = Y)

## -----------------------------------------------------------------------------
Gamma_1_hat <- emp_vario(data = X, k = 1, p = 0.8)

## -----------------------------------------------------------------------------
Gamma_hat <- emp_vario(data = X, p = 0.8)
Gamma_hat <- emp_vario(data = Y)

## ----out.width='33%', fig.align='', fig.show='hold', nomargin=TRUE------------
# A tree graph, a decomposable graph, and a non-decomposable (cycle) graph
plot(graph_from_literal(1--2--3, 2--4))
plot(graph_from_literal(1--2--3--1, 2--4))
plot(graph_from_literal(1--2--3--4--1))

## -----------------------------------------------------------------------------
set.seed(42)

my_model <- generate_random_model(d = 4, graph_type = "tree")
Ysim <- rmpareto_tree(
    n = 100,
    model = "HR",
    tree = my_model$graph,
    par = my_model$Gamma
)

theta_vec <- c(.2, .8, .3)
Ysim <- rmpareto_tree(
    n = 100,
    model = "logistic",
    tree = my_model$graph,
    par = theta_vec
)

## -----------------------------------------------------------------------------
# Utility function to plot fitted parameters against true/empirical ones
plot_fitted_params <- function(G0, G1, xlab = 'True', ylab = 'Fitted'){
    return(
        ggplot()
        + geom_point(aes(
            x = G0[upper.tri(G0)],
            y = G1[upper.tri(G1)]
        ))
        + geom_abline(slope = 1, intercept = 0)
        + xlab(xlab)
        + ylab(ylab)
    )
}

## -----------------------------------------------------------------------------
# Get flow graph from package
flow_graph <- getDanubeFlowGraph()

# Fit parameter matrix with this graph structure
flow_Gamma <- fmpareto_graph_HR(data = Y, graph = flow_graph)

# Compute likelihood/ICs and plot parameters
flow_loglik <- loglik_HR(
    data = Y,
    Gamma = flow_Gamma,
    graph = flow_graph
)
plot_fitted_params(chi_hat, Gamma2chi(flow_Gamma), xlab = 'Empirical')

## ----nomargin=TRUE, fig.align='', fig.show='hold', out.width='50%'------------
# Fit tree grpah to the data
emst_fit <- emst(data = Y, method = "vario")

# Compute likelihood/ICs, and plot fitted graph, parameters
loglik_emst <- loglik_HR(
    data = Y,
    Gamma = emst_fit$Gamma,
    graph = emst_fit$graph
)
plotDanubeIGraph(graph = emst_fit$graph)
plot_fitted_params(chi_hat, Gamma2chi(emst_fit$Gamma), xlab = 'Empirical')

## ----nomargin=TRUE------------------------------------------------------------
Gamma <- cbind(
    c(0, 1.5, 1.5, 2),
    c(1.5, 0, 2, 1.5),
    c(1.5, 2, 0, 1.5),
    c(2, 1.5, 1.5, 0)
)
Gamma2Sigma(Gamma, k = 1)
Gamma2Theta(Gamma)
Gamma2graph(Gamma)

## ----nomargin=TRUE------------------------------------------------------------
set.seed(42)
my_tree <- generate_random_tree(d = 4)
Gamma_vec <- c(.5, 1.4, .8)
Gamma_comp <- complete_Gamma(Gamma = Gamma_vec, graph = my_tree)
print(Gamma_comp)
plot(Gamma2graph(Gamma_comp))

## ----nomargin=TRUE------------------------------------------------------------
G <- rbind(
    c(0, 5, 7, 6, NA),
    c(5, 0, 14, 15, NA),
    c(7, 14, 0, 5, 5),
    c(6, 15, 5, 0, 6),
    c(NA, NA, 5, 6, 0)
)
Gamma_comp <- complete_Gamma(Gamma = G)
plot(Gamma2graph(Gamma_comp))

## ----nomargin=TRUE------------------------------------------------------------
set.seed(42)
G <- rbind(
    c(0, 5, 7, 6, 6),
    c(5, 0, 14, 15, 13),
    c(7, 14, 0, 5, 5),
    c(6, 15, 5, 0, 6),
    c(6, 13, 5, 6, 0)
)
my_graph <- generate_random_connected_graph(d = 5, m = 5)
Gamma_comp <- complete_Gamma(Gamma = G, graph = my_graph)
plot(Gamma2graph(Gamma_comp))

## ----nomargin=TRUE------------------------------------------------------------
set.seed(42)
d <- 10
my_model <- generate_random_model(d = d, graph_type = "decomposable")
plot(my_model$graph)
Ysim <- rmpareto(n = 100, d = d, model = "HR", par = my_model$Gamma)
my_fit <- fmpareto_graph_HR(data = Ysim, graph = my_model$graph, p = NULL)
plot_fitted_params(Gamma2chi(my_model$Gamma), Gamma2chi(my_fit))

## ----nomargin=TRUE------------------------------------------------------------
set.seed(1)
d <- 20
my_model <- generate_random_model(d = d, graph_type = "general")
plot(my_model$graph)
Ysim <- rmpareto(n = 100, d = d, model = "HR", par = my_model$Gamma)
my_fit <- fmpareto_graph_HR(data = Ysim, graph = my_model$graph, p = NULL, method = "vario")
plot_fitted_params(Gamma2chi(my_model$Gamma), Gamma2chi(my_fit))

## ----message=FALSE, warning=FALSE---------------------------------------------
# Run eglearn for a suitable list of penalization parameters
rholist <- seq(0, 0.1, length.out = 11)
eglearn_fit <- eglearn(
    data = Y,
    rholist = rholist,
    complete_Gamma = TRUE
)

# Compute the corresponding likelihoods/ICs
logliks_eglearn <- sapply(
    seq_along(rholist),
    FUN = function(j) {
        Gamma <- eglearn_fit$Gamma[[j]]
        if(is.null(Gamma)) return(c(NA, NA, NA))
        loglik_HR(
            data = Y,
            Gamma = Gamma,
            graph = eglearn_fit$graph[[j]]
        )
    }
)

# Plot the BIC vs rho/number of edges
ggplot(mapping = aes(x = rholist, y = logliks_eglearn['bic', ])) +
    geom_line() +
    geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
    geom_hline(aes(yintercept = flow_loglik['bic']), lty = "dashed") +
    geom_hline(aes(yintercept = loglik_emst['bic']), lty = "dotted") +
    xlab("rho") +
    ylab("BIC") +
    scale_x_continuous(
        breaks = rholist,
        labels = round(rholist, 3),
        sec.axis = sec_axis(
            trans = ~., breaks = rholist,
            labels = sapply(eglearn_fit$graph, igraph::gsize),
            name = "Number of edges"
        )
    )

## ----nomargin=TRUE, fig.align='', fig.show='hold', out.width='50%'------------
# Compare fitted chi to empirical one
best_index <- which.min(logliks_eglearn['bic', ])
best_Gamma <- eglearn_fit$Gamma[[best_index]]
best_graph <- eglearn_fit$graph[[best_index]]
plotDanubeIGraph(graph = best_graph)
plot_fitted_params(chi_hat, Gamma2chi(best_Gamma), xlab='Empirical')

