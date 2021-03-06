% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Creation_v0.6.R
\name{GenerateLandcover}
\alias{GenerateLandcover}
\title{GenerateLandcover creates a ascii land cover dataset with a corresponding patch number data set.}
\usage{
GenerateLandcover(
  NumberSeeds,
  LCVPath,
  PatchPath,
  runno = NULL,
  MatrixHabs = 10,
  GenExtent = 250,
  FullExtent = 1000,
  BufferFullExtent = 0,
  cover = 0.1,
  DistributionPatches = "Random",
  DistributionSeeds = "Uniform",
  NumberEachMatrixHabs = 1:200,
  MatrixHabsEqualSize = FALSE,
  MatrixPatchDistribution = "Random",
  MarixSeedDistribution = "Uniform"
)
}
\arguments{
\item{NumberSeeds}{integer -}

\item{LCVPath}{character - string of the output path of the land cover ascii}

\item{PatchPath}{character - string of the output path of the patch ascii}

\item{runno}{integer - defaults to NULL}

\item{MatrixHabs}{integer -number of habitats in the matrix}

\item{GenExtent}{integer - default 250, this is the extent that the land cover is generated at. Can be slow with large values.}

\item{FullExtent}{integer final extent of the landcover}

\item{BufferFullExtent}{integer - buffer around the edge of the land cover}

\item{cover}{floating point between 0 and 1.}

\item{DistributionPatches}{character - how seed points are distributed to grow the patches. Options are "Random" or "Poisson", "Clustered" or "Thomas", or "Regular", "Strauss" or"Dispersed"
and defaults to random.}

\item{DistributionSeeds}{character - a distribution from which to pick the size of patches from. Opetions are
"Gaussian" or"Normal", "Pareto", or the default "Uniform". Uniform gives patches of almost any size and the maximum variation.}

\item{NumberEachMatrixHabs}{list of numbers - default is 1:200}

\item{MatrixHabsEqualSize}{boolean - true or false}

\item{MatrixPatchDistribution}{character - how seed points are distributed to grow the patches. Options are "Random" or "Poisson", "Clustered" or "Thomas", or "Regular", "Strauss" or"Dispersed"
and defaults to random.}

\item{MarixSeedDistribution}{character - a distribution from which to pick the size of patches from. Opetions are
"Gaussian" or"Normal", "Pareto", or the default "Uniform". Uniform gives patches of almost any size and the maximum variation.}
}
\value{
list of arrays - the list contains two arrays, the "Patches" and "InterestLCV" (land cover).
}
\description{
This function uses other functions to create a land cover with the focal habitat haveing a number
of patches.
}
\examples{
   library(raster)
   library(LcvGen)

   NewLandcover <- GenerateLandcover(5, "Lcv5.asc", "Patches5.asc")
   plot(raster(NewLandcover$LCV))
   plot(raster(NewLandcover$Patches))

   NewLandcover <- GenerateLandcover(4, "Lcv5.asc", "Patches5.asc",
                                     DistributionPatches = "Dispersed",
                                     NumberEachMatrixHabs = 4,
                                     MatrixPatchDistribution = "Dispersed")

}
