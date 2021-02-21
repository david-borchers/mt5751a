#' @name tigerch
#' 
#' @title Simulated capture-recapture tiger survey dataset.
#' 
#' @docType data
#' @description Simulated capture history from a simulated population of 81 tigers in the Sundarbans
#' region of India. Simulation is based on a real camera trap survey of the tigers. (Object 
#' \link{tigersecrch} contains the same data, but with spatial capture information.)
#' @usage tigerch
#' @format A capture history matrix of 56 animals (rows) on 5 occasions (columns), with 0 representing
#' not captured and 1 representing captured.
#' @source Simulated
#' @examples
#'  data(tigerch)
NULL


#' @name tigersecrch
#' 
#' @title Simulated spatial capture-recapture tiger survey dataset.
#' 
#' @docType data
#' @description Simulated spatial capture history from a simulated population of 81 tigers in the Sundarbans
#' region of India. Simulation is based on a real camera trap survey of the tigers. (Object 
#' \link{tigerch} contains the same capture history data, but with spatial capture information removed.)
#' @usage tigersecrch
#' @format An \code{secr} \code{capthist} object containing 56 animals on 5 occasions.
#' @seealso \code{\link{tigermask}} 
#' @source Simulated
#' @examples
#'  data(tigersecrch)
NULL




#' @name tigermask
#' 
#' @title The \code{secr} \code{mask} object for simualted tiger survey.
#' 
#' @docType data
#' @description The \code{secr} \code{mask} object for the \code{tigersecrech} object.
#' @usage tigermask
#' @format An \code{secr} \code{mask} object defining the survey region.
#' 
#' @source Simulated
#' @examples
#'  data(tigermask)
NULL


#' @name tigerpop
#' 
#' @title A simulated tiger population.
#' 
#' @docType data
#' @description The \code{secr} \code{popn} object defining the population that was 
#' surveyed to get the capture history \code{tigersecrch}.
#' @usage tigerpop
#' @format An \code{secr} \code{mask} object defining the survey region.
#' 
#' @source Simulated
#' @examples
#'  data(tigerpop)
NULL


#' @name edwards.eberhardt
#' 
#' @title Rabbit capture-recapture data
#' 
#' @docType data
#' @description A capture-recapture data set on rabbits derived from Edwards and Eberhardt (1967) that accompanies 
#' programme MARK as an example analysis using the closed population models.
#' @usage edwards.eberhardt
#' @format A capture history matrix of 76 rabbits (rows) on 18 occasions (columns), with 0 representing
#' not captured and 1 representing captured.
#' @source Edwards, W.R. and L.L. Eberhardt 1967. Estimating cottontail abundance from live trapping data. J. Wildl. Manage. 31:87-96.
#' @examples
#'  data(edwards.eberhardt)
NULL


#' @name bowhead_LT
#' 
#' @title Line transect survey data from bowhead whale survey.
#' 
#' @docType data
#' @description Line transect survey data from the bowhead whale survey reported in Rekdal et al. (2015).
#' @usage bowhead_LT
#' @format A data frame of line transect data suitable for analysis with package \code{Distance}.
#' @source Rekdal, S.L., Hansen, R.G., Borchers, D., Bachmann, L., Laidre, K.L., Wiig, O., Nielsen, N.H., 
#' Fossette, S., Tervo, O. and Heide-Jorgensen, M.P., 2015. Trends in bowhead whales in West Greenland: 
#' Aerial surveys vs. genetic capture-recapture analyses. Marine Mammal Science, 31(1), pp.133-154. 
#' @examples
#'  data(bowhead_LT)
NULL



#' @name bowhead_PS
#' 
#' @title Plot sampling survey data from Bowhead whale survey.
#' 
#' @docType data
#' @description Data from the bowhead whale survey reported in Rekdal et al. (2015), but formatted as 
#' a plot sampling survey.
#' @usage bowhead_PS
#' @format A data frame with one row per plot and columns \code{stratum}, \code{transect} (plot), \code{A} 
#' (stratum area), \code{a} (plot area), and \code{n} (count of whale groups in plot).
#' @source Rekdal, S.L., Hansen, R.G., Borchers, D., Bachmann, L., Laidre, K.L., Wiig, O., Nielsen, N.H., 
#' Fossette, S., Tervo, O. and Heide-Jorgensen, M.P., 2015. Trends in bowhead whales in West Greenland: 
#' Aerial surveys vs. genetic capture-recapture analyses. Marine Mammal Science, 31(1), pp.133-154. 
#' @examples
#'  data(bowhead_PS)
NULL


#' @name plotsample1
#' 
#' @title Simulated plot sample data.
#' 
#' @docType data
#' @description Simulated plot sample data from 20 plots in a 100m by 100m survey region.
#' @usage plotsample1
#' @format A data frame with one row per plot and columns \code{stratum}, \code{plot}, \code{A} 
#' (stratum area), \code{a} (plot area), and \code{n} (count of whale groups in plot).
#' @source Simulated 
#' @examples
#'  data(plotsample1)
NULL



#' @name plotsample2
#' 
#' @title Simulated plot sample data.
#' 
#' @docType data
#' @description Simulated plot sample data from 20 plots in a 100m by 100m survey region.
#' @usage plotsample2
#' @format A data frame with one row per plot and columns \code{stratum}, \code{plot}, \code{A} 
#' (stratum area), \code{a} (plot area), and \code{n} (count of whale groups in plot).
#' @source Simulated 
#' @examples
#'  data(plotsample2)
NULL

