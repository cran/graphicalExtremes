## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = '#>',
    error = TRUE
)

## ----setup, include=FALSE-----------------------------------------------------
library(graphicalExtremes)
library(dplyr)
library(ggplot2)
theme_set(theme_bw() +
    theme(
        plot.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = "white"),
        plot.caption = element_text(
            size = 7.5, hjust = 0,
            margin = margin(t = 15)
        ),
        text = element_text(size = 11),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size = 0.25)
    )
)

## ---- fig.show='hold', out.width="50%"----------------------------------------
danube <- graphicalExtremes::danube

X <- danube$data_clustered

ggplot(X %>% as_tibble()) +
    geom_point(aes(x = X1, y = X2))
ggplot(X %>% as_tibble()) +
    geom_point(aes(x = X22, y = X28))

## ---- fig.align='center'------------------------------------------------------
danube_flow <- igraph::graph_from_edgelist(danube$flow_edges)
danube_plot_coords <- as.matrix(danube$info[, c("PlotCoordX", "PlotCoordY")])
plot(danube_flow, layout = danube_plot_coords, edge.arrow.size = .3)

## ---- fig.align='center'------------------------------------------------------
Y <- data2mpareto(data = X, p = .8)

## ---- fig.align='center'------------------------------------------------------
G <- cbind(c(0, 1.5), c(1.5, 0))
Ysim <- rmpareto(n = 100, model = "HR", d = 2, par = G)

theta <- .5
Ysim <- rmpareto(n = 100, model = "logistic", d = 2, par = theta)

## -----------------------------------------------------------------------------
chi_hat <- emp_chi(data = X, p = .8)

## -----------------------------------------------------------------------------
chi_hat <- emp_chi(data = Y)

## -----------------------------------------------------------------------------
Gamma_1_hat <- emp_vario(data = X, k = 1, p = 0.8)

## -----------------------------------------------------------------------------
Gamma_hat <- emp_vario(data = X, p = 0.8)
Gamma_hat <- emp_vario(data = Y)

## ---- echo = FALSE, fig.align='center', out.width="50%"-----------------------
# All defaults
knitr::include_graphics("images/graphs.png")

## -----------------------------------------------------------------------------
set.seed(42)

my_model <- generate_random_model(d = 4, graph_type = "tree")
Ysim <- rmpareto_tree(
    n = 100, model = "HR",
    tree = my_model$graph,
    par = my_model$Gamma
)

theta_vec <- c(.2, .8, .3)
Ysim <- rmpareto_tree(
    n = 100, model = "logistic",
    tree = my_model$graph,
    par = theta_vec
)

## ---- fig.align='center'------------------------------------------------------
danube_tree <- igraph::as.undirected(danube_flow)
danube_flow_fit <- fmpareto_graph_HR(data = X, graph = danube_tree, p = .8)

## ---- fig.align='center'------------------------------------------------------
danube_emst_fit <- emst(data = X, p = .8, method = "vario")
plot(danube_emst_fit$graph, layout = danube_plot_coords)

## ---- fig.align='center'------------------------------------------------------
ggplot() +
    geom_point(aes(x = c(Gamma2chi(danube_flow_fit)), y = c(chi_hat))) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Fitted") +
    ylab("Empirical") +
    coord_cartesian(xlim = c(0.4, 1), ylim = c(0.4, 1))

## ---- eval=FALSE--------------------------------------------------------------
#  Gamma <- cbind(
#      c(0, 1.5, 1.5, 2),
#      c(1.5, 0, 2, 1.5),
#      c(1.5, 2, 0, 1.5),
#      c(2, 1.5, 1.5, 0)
#  )
#  Gamma2Sigma(Gamma, k = 1)
#  Gamma2Theta(Gamma)
#  my_graph <- Gamma2graph(Gamma)

## ---- fig.align='center', echo=FALSE------------------------------------------
Gamma <- cbind(
    c(0, 1.5, 1.5, 2),
    c(1.5, 0, 2, 1.5),
    c(1.5, 2, 0, 1.5),
    c(2, 1.5, 1.5, 0)
)
Gamma2Sigma(Gamma, k = 1)
round(Gamma2Theta(Gamma), 3)
my_graph <- Gamma2graph(Gamma)

## ---- fig.align='center'------------------------------------------------------
set.seed(42)
my_tree <- generate_random_model(d = 4, graph_type = "tree")$graph
Gamma_vec <- c(.5, 1.4, .8)
complete_Gamma(Gamma = Gamma_vec, graph = my_tree)
plot(my_tree)

## ---- fig.align='center'------------------------------------------------------
G <- rbind(
    c(0, 5, 7, 6, NA),
    c(5, 0, 14, 15, NA),
    c(7, 14, 0, 5, 5),
    c(6, 15, 5, 0, 6),
    c(NA, NA, 5, 6, 0)
)
complete_Gamma(Gamma = G)
my_graph <- Gamma2graph(complete_Gamma(G))

## ---- fig.align='center'------------------------------------------------------
set.seed(42)
G <- rbind(
    c(0, 5, 7, 6, 6),
    c(5, 0, 14, 15, 13),
    c(7, 14, 0, 5, 5),
    c(6, 15, 5, 0, 6),
    c(6, 13, 5, 6, 0)
)
my_graph <- generate_random_model(d = 5, graph_type = "general", m = 5)$graph
complete_Gamma(Gamma = G, graph = my_graph)
plot(my_graph)

## ---- fig.align='center'------------------------------------------------------
set.seed(42)
d <- 10
my_model <- generate_random_model(d = d, graph_type = "decomposable")
plot(my_model$graph)
Ysim <- rmpareto(n = 100, d = d, model = "HR", par = my_model$Gamma)
my_fit <- fmpareto_graph_HR(data = Ysim, graph = my_model$graph, p = NULL)
ggplot() +
    geom_point(aes(x = c(Gamma2chi(my_fit)), y = c(Gamma2chi(my_model$Gamma)))) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Fitted") +
    ylab("True") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

## ---- fig.align='center'------------------------------------------------------
set.seed(42)
d <- 20
my_model <- generate_random_model(d = d, graph_type = "general")
plot(my_model$graph)
Ysim <- rmpareto(n = 100, d = d, model = "HR", par = my_model$Gamma)
my_fit <- fmpareto_graph_HR(data = Ysim, graph = my_model$graph, p = NULL, method = "vario")
ggplot() +
    geom_point(aes(x = c(Gamma2chi(my_fit)), y = c(Gamma2chi(my_model$Gamma)))) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Fitted") +
    ylab("True") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

## ---- fig.align='center', message=FALSE, warning=FALSE------------------------
set.seed(42)
# (Ysim from previous example is used)
rholist <- seq(1e-5, 0.04, length.out = 10)
eglasso_fit <- eglearn(
    data = Ysim,
    rholist = rholist,
    complete_Gamma = TRUE
)

eglasso_loglik <- sapply(1:length(rholist), FUN = function(j) {
    loglik_HR(
        data = Ysim, p = .8, Gamma = eglasso_fit$Gamma[[j]],
        graph = eglasso_fit$graph[[j]]
    )
})

emst_fit <- emst(Ysim)

danube_loglik_tree <- loglik_HR(
    data = Ysim, p = .8, Gamma = emst_fit$Gamma,
    graph = emst_fit$graph
)

ggplot(mapping = aes(x = rholist, y = eglasso_loglik[3, ])) +
    geom_line() +
    geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
    geom_hline(aes(yintercept = danube_loglik_tree[3]), lty = "dashed") +
    xlab("rho") +
    ylab("BIC") +
    scale_x_continuous(
        breaks = rholist,
        labels = round(rholist, 3),
        sec.axis = sec_axis(
            trans = ~., breaks = rholist,
            labels = sapply(eglasso_fit$graph, igraph::gsize),
            name = "Number of edges"
        )
    )

## ---- fig.align='center', message=FALSE, warning=FALSE------------------------
best_Gamma <- eglasso_fit$Gamma[[which.min(eglasso_loglik[3, ])]]
chi_hat <- emp_chi(Ysim)

ggplot() +
    geom_point(aes(x = c(Gamma2chi(best_Gamma)), y = c(chi_hat))) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Fitted") +
    ylab("Empirical") +
    coord_cartesian(
        ylim = c(0, 1),
        xlim = c(0, 1)
    )

