## ----echo = TRUE, message=FALSE, warning=FALSE--------------------------------
library(graphicalExtremes)
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
        plot.caption=element_text( size=7.5, hjust=0, margin=margin(t=15)),
        text = element_text(size = 11),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size = 0.25)
    )
)


## ----fig.width=6, fig.height=3------------------------------------------------
plotFlights(plotConnections = FALSE, map = "world", xyRatio = 2)

## ----fig.width=6, fig.height=4------------------------------------------------
# Get IATAs from Texas cluster
IATAs <- getFlightDelayData('IATAs', 'tcCluster')

# Plot airports + connections (indicating number of flights by thickness)
plotFlights(
    IATAs,
    useAirportNFlights = TRUE,
    useConnectionNFlights = TRUE
)

## -----------------------------------------------------------------------------
# Check whether all dates from the train-test-split are available
# (due to size restrictions, the CRAN version of this package does not contain all dates)
allDatesAvailable <- tryCatch({
    getFlightDelayData('dates', dateFilter = c('tcAll'))
    TRUE
}, error = function(...) FALSE)
cat('All dates avilable:', allDatesAvailable, '\n')

# Create train and test data sets
if(allDatesAvailable){
    # Use train-test-split and threshold p from article
    matEst <- drop(getFlightDelayData('delays', 'tcCluster', 'tcTrain'))
    matVal <- drop(getFlightDelayData('delays', 'tcCluster', 'tcTest'))
    p <- 0.95
} else {
    # Take all available dates that do not contain NAs
    mat <- drop(getFlightDelayData('delays', 'tcCluster', 'all'))
    rowHasNA <- apply(is.na(mat), 1, any)
    mat <- mat[!rowHasNA, ]
    
    # Split dates in half
    splitInd <- floor(nrow(mat)/2)
    matEst <- mat[1:splitInd,]
    matVal <- mat[(splitInd+1):nrow(mat),]
    
    # Use a lower threshold to compensate for the smaller dataset
    p <- 0.8
}

## -----------------------------------------------------------------------------
# Compute the empirical extremal correlation matrix
emp_chi_mat <- emp_chi(matEst, p = p)

# Utility function to plot fitted parameters against true/empirical ones
plot_fitted_params <- function(G0, G1, xlab = 'Empirical', ylab = 'Fitted'){
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
# Compute undirected flights per connection between selected airports
nYears <- dim(flights$flightCounts)[3]
flightsPerConnection <- apply(flights$flightCounts[IATAs,IATAs,], c(1, 2), sum)
flightsPerConnectionUD <- flightsPerConnection + t(flightsPerConnection)

# Make flight graph from adjacency matrix
A <- 1 * (flightsPerConnectionUD > nYears * 12)
flight_graph <- graph_from_adjacency_matrix(A, diag = FALSE, mode = "undirected")

# Plot flight graph
plotFlights(IATAs, graph = flight_graph, clipMap = 1.3, xyRatio = 1)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Fit the model
model_fit <- fmpareto_graph_HR(
    data = matEst,
    graph = flight_graph,
    p = p,
    method = "vario"
)

# Compute likelihood/ICs
flights_loglik_graph <- loglik_HR(
    data = matVal,
    p = p,
    graph = flight_graph,
    Gamma = model_fit
)
cat("Flight graph test-loglikelihood =", round(flights_loglik_graph['loglik'], 2), "\n")

# Plot fitted parameters
plot_fitted_params(emp_chi_mat, Gamma2chi(model_fit))

## ----collapse=TRUE, fig.align='', fig.show='hold', out.width='50%'------------
# Fit the model
flights_emst_fit <- emst(data = matEst, p = p, method = "vario")

# Compute likelihood/ICs
flights_loglik_tree <- loglik_HR(
    data = matVal,
    p = p,
    Gamma = flights_emst_fit$Gamma,
    graph = flights_emst_fit$graph
)
cat("Tree test-loglikelihood =", round(flights_loglik_tree['likelihood'], 2), "\n")

# Plot fitted graph, parameters
plotFlights(
    IATAs,
    graph = flights_emst_fit$graph,
    xyRatio = 1,
    clipMap = 1.3
)
plot_fitted_params(emp_chi_mat, Gamma2chi(flights_emst_fit$Gamma))

## ----message=FALSE, warning=FALSE---------------------------------------------
# Fit the model
rholist <- seq(0, 0.50, length.out = 11)
flights_eglearn_fit <- eglearn(matEst, p = p, rholist = rholist, complete_Gamma = TRUE)

# Compute likelihood/ICs
flights_loglik <- sapply(seq_along(rholist), FUN = function(j) {
    loglik_HR(
        data = matVal,
        p = p,
        Gamma = flights_eglearn_fit$Gamma[[j]],
        graph = flights_eglearn_fit$graph[[j]]
    )
})

## -----------------------------------------------------------------------------
ggplot(
    mapping = aes(x = rholist, y = flights_loglik['loglik', ])) +
    geom_line() +
    geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
    geom_hline(aes(yintercept = flights_loglik_graph['loglik']), lty = "dashed") +
    geom_hline(aes(yintercept = flights_loglik_tree['loglik']), lty = "dotted") +
    xlab("rho") +
    ylab("Log-likelihood") +
    scale_x_continuous(
        breaks = rholist,
        labels = round(rholist, 3),
        sec.axis = sec_axis(
            trans = ~., breaks = rholist,
            labels = sapply(
                flights_eglearn_fit$graph,
                igraph::gsize
            ),
            name = "Number of edges"
        )
)

best_index <- which.max(flights_loglik['loglik',])
cat('Best rho =', rholist[best_index], '\n')
cat('Corresponding test-loglikelihood =', flights_loglik['loglik', best_index])

## ----fig.align='', fig.show='hold', out.width='50%'---------------------------
best_graph <- flights_eglearn_fit$graph[[best_index]]
best_Gamma <- flights_eglearn_fit$Gamma[[best_index]]
plotFlights(IATAs, graph = best_graph, clipMap = 1.3, xyRatio = 1)
plot_fitted_params(emp_chi_mat, Gamma2chi(best_Gamma))

