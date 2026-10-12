#' Traits of 220 populations of the Ranunculus auricomus plant species complex (Karbstein et al. 2021)
#'
#' A dataset containing reproductive and cytogenetic traits,
#' and environmental characteristics of 220
#' plant populations comprising 80 taxa. 
#' It has been published by Karbstein et al. 2021 
#' in Molecular Ecology (https://onlinelibrary.wiley.com/doi/10.1111/mec.15919).
#'
#'
#' @format A data frame with 220 rows and 12 variables:
#' \itemize{
#'   \item population_id: the official population identifier
#'   \item sexual_seed_perc: the percent of sexual to all formed seeds averaged
#'   across all flowers within the population
#'   \item sexual_seed_prop: the ratio of sexual to all formed seeds averaged
#'   across all flowers within the population
#'   \item altitude: meters above sea level
#'   \item latitude_gps_n: latitude coordinate of the measurement location
#'   (Northern Hemisphere)
#'   \item longitude_gps_e: longitude coordinate of the measurement location
#'   (Eastern Hemisphere)
#'   \item solar_radiation: amount of incoming solar radiation at the location
#'   (kJ m⁻² day⁻¹; see BioClim database)
#'   \item habitat: description of the environment where the measurement was
#'   taken (e.g., forest, meadow, or swamp)
#'   \item ploidy: the number of entire sets of chromosomes in a cell
#'   \item number_flower_leaves: number of flower or floral leaves of a plant
#'   \item annual_mean_temperature: annual mean temperature at the location
#'   (see 'Bio1' in BioClim database)
#'   \item isothermality: ratio of mean diurnal range divided by temperature
#'   annual range (see 'Bio3' in BioClim database)
#' }
#' @docType data
#' @usage data(plants)
#' @keywords datasets
"plants"
